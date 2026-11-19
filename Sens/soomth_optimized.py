import os
import time
import numpy as np
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from scipy.interpolate import griddata
from scipy.signal import fftconvolve
from concurrent.futures import ProcessPoolExecutor, as_completed


# ===================== User parameters =====================
TRAVEL_FILE = "../../Input/array_p.txt"
FORWARD_KERNEL_TEMPLATE = "../../Output/Kernel/Kernel_vel_s{Ns}_r{Nr}.dat"
RESAMPLING_KERNEL_TEMPLATE = "../../Output/Kernel/Kernel_interpolated_s{Ns}_r{Nr}.dat"
TARGET_KERNEL_TEMPLATE = "../../Output/Kernel/Kernel_soomth_s{Ns}_r{Nr}.dat"
OUTPUT_KERNEL_TEMPLATE = "../../Output/Kernel/Kernel_interpolated_s{Ns}_r{Nr}.dat"

SMOOTH_TEMPLATE = "../../Output/Kernel/Smooth_s{Ns}_r{Nr}.dat"
SENSITIVITY_TEMPLATE = "../../Output/Kernel/Sensitivity_s{Ns}_r{Nr}.dat"

CONFIG_FILE = "../../Input/parameters_inversion3d.config"
MESH_SUMMARY_FILE = "../../Output/Mesh/mesh_summary.txt"
ARRIVAL_FILE = "../../Input/arrival_p.txt"

OUTPUT_FMT = "%.15e %.15e %.15e"

# 并行进程数。None 表示自动使用配置文件中的 num_cores。
N_WORKERS = None

# 只处理前 N 条数据。None 表示全部处理。
MAX_TASKS = None

# 密网格加密倍数。
NUM_DENSE = 20

# 二维剖面基础网格节点数，独立于三维反演网格的 nx_inv、nz_inv。
PROFILE_NX = 5
PROFILE_NZ = 5

# 高斯平滑截断半径。原程序使用 10.0，但数值上通常 3.0 已足够。
THD_GAUSSIAN = 10.0
# ===========================================================

@dataclass
class SmoothingParameters:
    wavelength: float
    num_cores: int


def read_smoothing_parameters(filename):
    """兼容 17 行旧配置和含正演网格数的 18 行配置。"""
    with Path(filename).open(encoding="utf-8") as file:
        rows = [line.split("c:", 1)[0].split() for line in file]
    rows = [row for row in rows if row]
    if len(rows) not in (17, 18):
        raise ValueError(f"{filename}: 反演配置格式不正确（需要 17 或 18 行）")

    offset = len(rows) - 17  # 18 行版本在震源数前增加了一行正演网格数
    parameters = SmoothingParameters(
        wavelength=float(rows[7 + offset][0]),
        num_cores=int(rows[16 + offset][0]),
    )
    if parameters.wavelength <= 0 or parameters.num_cores < 1:
        raise ValueError(f"{filename}: 主波长或核数无效")
    return parameters


def read_mesh_summary(filename):
    """返回 (震源编号, 台站编号) 对应的 (模型长度, 模型深度)。"""
    rows = read_array(filename, 13)
    sizes = {}
    for row in rows:
        ns, nr = row[:2]
        length, depth = row[3:5]  # model_length、depth；不是 line_length
        if not np.isfinite(row[:5]).all() or ns != int(ns) or nr != int(nr) or length <= 0 or depth <= 0:
            raise ValueError(f"{filename}: 测线编号或模型范围无效: {row[:5]}")
        pair = (int(ns), int(nr))
        if pair in sizes:
            raise ValueError(f"{filename}: 重复的震源-台站编号 {pair}")
        sizes[pair] = (float(length), float(depth))
    return sizes


def read_array(filename, min_cols):
    filename = Path(filename)
    if not filename.exists():
        raise FileNotFoundError(f"cannot open this file -- {filename}")

    data = np.loadtxt(filename)
    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] < min_cols:
        raise ValueError(f"{filename} must contain at least {min_cols} columns")

    return data


def make_regular_grid(x_min, x_max, z_min, z_max, dx, dz):
    num_x = round((x_max - x_min) / dx)
    num_z = round((z_max - z_min) / dz)

    x = x_min + np.arange(num_x) * dx
    z = z_min + np.arange(num_z) * dz
    xx, zz = np.meshgrid(x, z, indexing="ij")

    grid = np.zeros((num_x * num_z, 3), dtype=float)
    grid[:, 0] = xx.ravel()
    grid[:, 1] = zz.ravel()
    return grid, num_x, num_z


def interpolate_kernel_to_target(array_ker_f, array_ker_t):
    source_points = array_ker_f[:, 0:2]
    source_values = array_ker_f[:, 2]
    target_points = array_ker_t[:, 0:2]

    linear_values = griddata(source_points, source_values, target_points, method="linear")

    bad = np.isnan(linear_values)
    if np.any(bad):
        linear_values[bad] = griddata(
            source_points, source_values, target_points[bad], method="nearest"
        )

    array_resampled = array_ker_t.copy()
    array_resampled[:, 2] = linear_values
    return array_resampled


def get_resampled_kernel(id_ns, id_nr, x_min, x_max, z_min, z_max, dx, dz):
    file_ker_f = FORWARD_KERNEL_TEMPLATE.format(Ns=id_ns, Nr=id_nr)
    array_ker_f = read_array(file_ker_f, min_cols=3)[:, 0:3].copy()

    # array_ker_f[:, 1] *= -1.0

    array_ker_t, _, _ = make_regular_grid(x_min, x_max, z_min, z_max, dx, dz)
    return interpolate_kernel_to_target(array_ker_f, array_ker_t)


@lru_cache(maxsize=16)
def _gaussian_kernel_norm(nx1, nz1, dx1, dz1, dd, thd):
    radius_x = int(np.ceil(thd * dd / dx1))
    radius_z = int(np.ceil(thd * dd / dz1))

    rx = np.arange(-radius_x, radius_x + 1, dtype=float)
    rz = np.arange(-radius_z, radius_z + 1, dtype=float)
    dist2 = (rx[:, None] * dx1)**2 + (rz[None, :] * dz1)**2

    fa = 4.0 / (np.pi * dd * dd)
    kernel = fa * np.exp(-4.0 * dist2 / (dd * dd))
    norm = fftconvolve(np.ones((nx1, nz1), dtype=float), kernel, mode="same")
    return kernel, norm


def gaussian_smooth_griddata(kernel_grid1, nx1, nz1, dx1, dz1, kernel_grid2, dd, thd):
    values1 = kernel_grid1[:, 2].reshape(nx1, nz1)
    kernel, ff2 = _gaussian_kernel_norm(nx1, nz1, dx1, dz1, dd, thd)

    ff1 = fftconvolve(values1, kernel, mode="same")
    smooth_values1 = np.divide(ff1, ff2, out=values1.copy(), where=ff2 != 0.0)

    # 当前 build_tasks 中 dx=dx1、dz=dz1，因此两个规则网格完全一致。
    # 这种情况下不再做第二次 griddata，这是最主要的提速点之一。
    if kernel_grid1.shape == kernel_grid2.shape and np.array_equal(kernel_grid1[:, :2], kernel_grid2[:, :2]):
        out = kernel_grid2.copy()
        out[:, 2] = smooth_values1.ravel()
        return out

    smooth_values2 = griddata(
        points=kernel_grid1[:, :2],
        values=smooth_values1.ravel(),
        xi=kernel_grid2[:, :2],
        method="linear",
        fill_value=0.0,
    )

    out = kernel_grid2.copy()
    out[:, 2] = smooth_values2
    return out


def process_one_task(task):
    (
        task_id,
        id_ns,
        id_nr,
        x_min,
        x_max,
        z_min,
        z_max,
        dx1,
        dz1,
        nx1,
        nz1,
        dx,
        dz,
        nx,
        nz,
        len_gaussian,
        thd_gaussian,
    ) = task

    output_file = TARGET_KERNEL_TEMPLATE.format(Ns=id_ns, Nr=id_nr)

    array_ker_f_resampled = get_resampled_kernel(
        id_ns, id_nr, x_min, x_max, z_min, z_max, dx1, dz1
    )

    array_ker_f_smooth, _, _ = make_regular_grid(x_min, x_max, z_min, z_max, dx, dz)
    array_ker_f_smooth = gaussian_smooth_griddata(
        array_ker_f_resampled,
        nx1,
        nz1,
        dx1,
        dz1,
        array_ker_f_smooth,
        len_gaussian,
        thd_gaussian,
    )

    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(output_file, array_ker_f_smooth, fmt=OUTPUT_FMT)

    return task_id, id_ns, id_nr, output_file


def build_tasks(parameters, array_data, mesh_sizes):
    len_gaussian = parameters.wavelength / 1.0
    thd_gaussian = THD_GAUSSIAN

    if PROFILE_NX < 2 or PROFILE_NZ < 2 or NUM_DENSE < 1:
        raise ValueError("PROFILE_NX、PROFILE_NZ 必须至少为 2，NUM_DENSE 必须至少为 1")

    n_task = array_data.shape[0] if MAX_TASKS is None else min(MAX_TASKS, array_data.shape[0])

    tasks = []
    for i in range(n_task):
        id_ns = int(array_data[i, 0])
        id_nr = int(array_data[i, 1])
        pair = (id_ns, id_nr)
        if pair not in mesh_sizes:
            raise ValueError(f"mesh_summary.txt 缺少震源-台站对 {pair}")

        model_length, depth = mesh_sizes[pair]
        x_min, x_max = 0.0, model_length
        z_min, z_max = -depth, 0.0
        dx0 = model_length / (PROFILE_NX - 1)
        dz0 = depth / (PROFILE_NZ - 1)
        dx1 = dx0 / NUM_DENSE
        dz1 = dz0 / NUM_DENSE
        dx, dz = dx1, dz1
        nx1 = (PROFILE_NX - 1) * NUM_DENSE
        nz1 = (PROFILE_NZ - 1) * NUM_DENSE
        nx, nz = nx1, nz1
        tasks.append(
            (
                i,
                id_ns,
                id_nr,
                x_min,
                x_max,
                z_min,
                z_max,
                dx1,
                dz1,
                nx1,
                nz1,
                dx,
                dz,
                nx,
                nz,
                len_gaussian,
                thd_gaussian,
            )
        )
    return tasks


def main():
    t0 = time.time()

    parameters = read_smoothing_parameters(CONFIG_FILE)
    mesh_sizes = read_mesh_summary(MESH_SUMMARY_FILE)
    array_data = read_array(ARRIVAL_FILE, 3)
    tasks = build_tasks(parameters, array_data, mesh_sizes)

    n_workers = parameters.num_cores if N_WORKERS is None else N_WORKERS
    n_workers = max(1, min(int(n_workers), os.cpu_count() or 1, len(tasks)))

    print(f"Total tasks - smooth : {len(tasks)}")
    print(f"Workers     - smooth : {n_workers}")

    if n_workers == 1:
        for task in tasks:
            task_id, id_ns, id_nr, output_file = process_one_task(task)
            # print(f"[{task_id + 1}/{len(tasks)}] Ns={id_ns}, Nr={id_nr}, output={output_file}")
    else:
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            future_to_task = {executor.submit(process_one_task, task): task for task in tasks}
            done = 0
            for future in as_completed(future_to_task):
                done += 1
                task = future_to_task[future]
                try:
                    task_id, id_ns, id_nr, output_file = future.result()
                    # print(f"[{done}/{len(tasks)}] Ns={id_ns}, Nr={id_nr}, output={output_file}")
                except Exception as exc:
                    print(f"[{done}/{len(tasks)}] task_id={task[0]} failed: {exc}")
                    raise

    # print(f"Elapsed time: {time.time() - t0:.2f} s")


if __name__ == "__main__":
    main()
