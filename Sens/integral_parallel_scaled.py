import numpy as np
import time
import os
from dataclasses import dataclass
from pathlib import Path
from scipy.signal import convolve2d
from typing import Union, List, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed

# 避免多进程 + BLAS/OpenMP 线程同时开启导致过度并行。
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

# ======================== 用户参数 ========================
CONFIG_FILE = "../../Input/parameters_inversion.config"
ERROR_FILE = "../../Output/Error/Error.txt"

KERNEL_SOOMTH = "../../Output/Kernel/Kernel_soomth_s{Ns}_r{Nr}.dat"
KERNEL_INTEGRAL = "../../Output/Kernel/Kernel_integral_s{Ns}_r{Nr}.dat"

SENSITIVITY_FILE = "../../Output/Model/Inver/Sensitivity.dat"
SENSITIVITY_IN_FILE = "../../Output/Model/Inver/Sensitivity_in.dat"

OUTPUT_FMT = "%.15e %.15e %.15e"

# ================= 敏感核幅值均衡参数 =================
SCALE_KERNELS = True
SCALE_REFERENCE = "median"   # "median" or "mean"
SCALE_POWER = 1.0            # 1.0=完全均衡，0.5=较温和
MIN_SCALE_FACTOR = 0.01
MAX_SCALE_FACTOR = 100.0

SAVE_SCALED_KERNELS = True
KERNEL_INTEGRAL_SCALED = "../../Output/Kernel/Kernel_integral_scaled_s{Ns}_r{Nr}.dat"
KERNEL_SCALE_INFO = "../../Output/Kernel/Kernel_scale_info.dat"
# =======================================================

# None：使用 parameters_inversion.config 中的 num_cores
# 也可以直接写成 4、8、16 等。
N_WORKERS = None

# None：处理全部数据；调试时可设为 10 等。
MAX_TASKS = None

# =========================================================


@dataclass
class InversionConfig:
    x1: float
    x2: float
    z1: float
    z2: float
    nx: int
    nz: int
    num_sources: int
    num_receivers: int
    num_iterations: int
    max_perturbation: float
    damping_factor: float
    smooth_factor: float
    num_cores: int
    wavelength: float
    min_window_cc: float
    min_waveform_cc: float
    integration_interval: int
    inversion_method: int


def read_inversion_config(config_file: Union[str, Path]) -> InversionConfig:
    path = Path(config_file)
    if not path.is_file():
        raise FileNotFoundError(f"配置文件不存在: {path}")

    rows: List[Tuple[int, List[float]]] = []

    with path.open("r", encoding="utf-8") as file:
        for line_number, line in enumerate(file, start=1):
            value_text = line.split("c:", 1)[0].strip()
            if not value_text:
                continue

            try:
                rows.append((line_number, [float(x) for x in value_text.split()]))
            except ValueError as exc:
                raise ValueError(
                    f"第 {line_number} 行包含无效的数字格式: {line.rstrip()}"
                ) from exc

    expected_counts = [2, 2, 2] + [1] * 12

    if len(rows) != len(expected_counts):
        raise ValueError(
            f"期望找到 {len(expected_counts)} 行参数，但实际解析了 {len(rows)} 行"
        )

    for (line_number, values), expected in zip(rows, expected_counts):
        if len(values) != expected:
            raise ValueError(
                f"第 {line_number} 行需要 {expected} 个参数，但实际找到 {len(values)} 个"
            )

    values = [row_values for _, row_values in rows]

    return InversionConfig(
        x1=values[0][0],
        x2=values[0][1],
        z1=values[1][0],
        z2=values[1][1],
        nx=int(values[2][0]),
        nz=int(values[2][1]),
        num_sources=int(values[3][0]),
        num_receivers=int(values[4][0]),
        num_iterations=int(values[5][0]),
        max_perturbation=values[6][0],
        damping_factor=values[7][0],
        smooth_factor=values[8][0],
        num_cores=int(values[9][0]),
        wavelength=values[10][0],
        min_window_cc=values[11][0],
        min_waveform_cc=values[12][0],
        integration_interval=int(values[13][0]),
        inversion_method=int(values[14][0]),
    )


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


def integrate_array(array_data, num_x, num_z, dx, dz):
    """利用二维卷积完成局部积分。"""
    values = array_data[:, 2].reshape(num_x, num_z)

    wx = np.array([[1.0, 2.0, 1.0]], dtype=np.float64) / 6.0 * dx
    wz = np.array([[1.0], [2.0], [1.0]], dtype=np.float64) / 6.0 * dz
    kernel = wz @ wx

    integral_values = convolve2d(
        values,
        kernel,
        mode="same",
        boundary="fill",
        fillvalue=0.0,
    )

    array_integral = array_data.copy()
    array_integral[:, 2] = integral_values.ravel()

    return array_integral


def project2invergrid(array_ker_inver1, dx1, dz1, array_integral, num_x1, num_z1):
    """利用双线性权重将积分网格投影到反演网格。"""
    array_data2 = array_ker_inver1.copy()
    array_data2[:, 2] = 0.0

    x_min = array_ker_inver1[0, 0]
    z_min = array_ker_inver1[0, 1]

    xx = array_integral[:, 0]
    zz = array_integral[:, 1]
    value = array_integral[:, 2]

    fx = (xx - x_min) / dx1
    fz = (zz - z_min) / dz1

    ix = np.floor(fx).astype(np.int32)
    iz = np.floor(fz).astype(np.int32)

    tx = fx - ix
    tz = fz - iz

    weights = [
        (0, 0, (1.0 - tx) * (1.0 - tz)),
        (1, 0, tx * (1.0 - tz)),
        (0, 1, (1.0 - tx) * tz),
        (1, 1, tx * tz),
    ]

    for dx_idx, dz_idx, weight in weights:
        ix_new = ix + dx_idx
        iz_new = iz + dz_idx

        mask = (
            (ix_new >= 0)
            & (ix_new < num_x1)
            & (iz_new >= 0)
            & (iz_new < num_z1)
        )

        target = ix_new[mask] * num_z1 + iz_new[mask]

        np.add.at(
            array_data2[:, 2],
            target,
            value[mask] * weight[mask],
        )

    return array_data2


# ===================== 并行进程全局数据 =====================

_G_ARRAY_KER_INVER1 = None
_G_DX1 = None
_G_DZ1 = None
_G_NUM_X1 = None
_G_NUM_Z1 = None


def init_worker(x_min1, z_min1, dx1, dz1, num_x1, num_z1):
    """
    每个子进程启动时只执行一次。
    在进程内部创建反演网格，避免每个任务重复传输大数组。
    """
    global _G_ARRAY_KER_INVER1
    global _G_DX1, _G_DZ1, _G_NUM_X1, _G_NUM_Z1

    _G_DX1 = dx1
    _G_DZ1 = dz1
    _G_NUM_X1 = num_x1
    _G_NUM_Z1 = num_z1

    i_idx, j_idx = np.meshgrid(
        np.arange(num_x1),
        np.arange(num_z1),
        indexing="ij",
    )

    _G_ARRAY_KER_INVER1 = np.zeros(
        (num_x1 * num_z1, 3),
        dtype=np.float64,
    )

    _G_ARRAY_KER_INVER1[:, 0] = x_min1 + dx1 * i_idx.ravel()
    _G_ARRAY_KER_INVER1[:, 1] = z_min1 + dz1 * j_idx.ravel()


def process_one_kernel(task):
    """
    单个震源-接收器核的完整处理：
        读取平滑核 -> 积分 -> 投影 -> 保存
    """
    task_id, id_ns, id_nr = task

    file_ker_smooth = KERNEL_SOOMTH.format(Ns=id_ns, Nr=id_nr)
    file_ker_integral = KERNEL_INTEGRAL.format(Ns=id_ns, Nr=id_nr)

    array_ker_smooth = read_array(file_ker_smooth, 3)

    x_min = np.min(array_ker_smooth[:, 0])
    x_max = np.max(array_ker_smooth[:, 0])
    z_min = np.min(array_ker_smooth[:, 1])
    z_max = np.max(array_ker_smooth[:, 1])

    dz2 = array_ker_smooth[1, 1] - array_ker_smooth[0, 1]
    if dz2 == 0.0:
        raise ValueError(f"{file_ker_smooth}: dz2 = 0")

    num_z2 = round((z_max - z_min) / dz2 + 1)

    if num_z2 >= array_ker_smooth.shape[0]:
        raise ValueError(
            f"{file_ker_smooth}: 无法根据数据确定 x 方向网格间距"
        )

    dx2 = array_ker_smooth[num_z2, 0] - array_ker_smooth[0, 0]
    if dx2 == 0.0:
        raise ValueError(f"{file_ker_smooth}: dx2 = 0")

    num_x2 = round((x_max - x_min) / dx2 + 1)

    expected = num_x2 * num_z2
    if expected != array_ker_smooth.shape[0]:
        raise ValueError(
            f"{file_ker_smooth}: 网格尺寸不匹配，"
            f"num_x2*num_z2={expected}, 实际点数={array_ker_smooth.shape[0]}"
        )

    array_integral = integrate_array(
        array_ker_smooth,
        num_x2,
        num_z2,
        dx2,
        dz2,
    )

    array_ker_inver2 = project2invergrid(
        _G_ARRAY_KER_INVER1,
        _G_DX1,
        _G_DZ1,
        array_integral,
        _G_NUM_X1,
        _G_NUM_Z1,
    )

    Path(file_ker_integral).parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(file_ker_integral, array_ker_inver2, fmt=OUTPUT_FMT)

    return task_id, id_ns, id_nr, file_ker_integral


def build_tasks(array_data):
    n_task = array_data.shape[0]
    if MAX_TASKS is not None:
        n_task = min(n_task, MAX_TASKS)

    tasks = []

    for i in range(n_task):
        id_ns = int(array_data[i, 0])
        id_nr = int(array_data[i, 1])
        tasks.append((i, id_ns, id_nr))

    return tasks


def get_kernel_max_abs(file_ker_integral):
    """返回单个积分敏感核第三列的最大绝对值。"""
    data = read_array(file_ker_integral, 3)
    return float(np.max(np.abs(data[:, 2])))


def build_kernel_scale_factors(array_data):
    """
    统计每个震源-台站对的 max(|K|)，并计算放缩系数。

    M_i   = max(|K_i|)
    M_ref = 所有非零 M_i 的中位数或平均值

    c_i = (M_ref / M_i) ** SCALE_POWER

    M_i < M_ref -> c_i > 1，放大
    M_i > M_ref -> c_i < 1，缩小
    """
    num_data = array_data.shape[0]
    if MAX_TASKS is not None:
        num_data = min(num_data, MAX_TASKS)

    info = np.zeros((num_data, 3), dtype=np.float64)

    for i in range(num_data):
        id_ns = int(array_data[i, 0])
        id_nr = int(array_data[i, 1])
        file_ker_integral = KERNEL_INTEGRAL.format(Ns=id_ns, Nr=id_nr)

        info[i, 0] = id_ns
        info[i, 1] = id_nr
        info[i, 2] = get_kernel_max_abs(file_ker_integral)

    positive = info[:, 2] > 0.0
    if not np.any(positive):
        raise ValueError("所有积分敏感核最大绝对值均为0，无法计算放缩系数。")

    if SCALE_REFERENCE == "median":
        reference = float(np.median(info[positive, 2]))
    elif SCALE_REFERENCE == "mean":
        reference = float(np.mean(info[positive, 2]))
    else:
        raise ValueError('SCALE_REFERENCE 只能是 "median" 或 "mean"。')

    scale_factors = np.ones(num_data, dtype=np.float64)

    for i in range(num_data):
        max_abs = info[i, 2]

        if max_abs > 0.0:
            scale = (reference / max_abs) ** SCALE_POWER
            scale = float(np.clip(scale, MIN_SCALE_FACTOR, MAX_SCALE_FACTOR))
        else:
            scale = 1.0

        scale_factors[i] = scale

    scaled_max_abs = info[:, 2] * scale_factors
    scale_info = np.column_stack(
        (info[:, 0], info[:, 1], info[:, 2], scale_factors, scaled_max_abs)
    )

    Path(KERNEL_SCALE_INFO).parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        KERNEL_SCALE_INFO,
        scale_info,
        fmt="%d %d %.15e %.15e %.15e",
    )

    print(f"Scale reference max(|K|): {reference:.15e}")
    # print(f"Scale info file         : {KERNEL_SCALE_INFO}")

    return scale_factors


def sum_integral_kernels(array_data, num_x1, num_z1):
    """
    所有积分核完成后：
      1. 统计每个核的 max(|K|)
      2. 计算放缩系数
      3. 小核放大、大核缩小
      4. 使用放缩后的核进行累加
    """
    num_data = array_data.shape[0]
    if MAX_TASKS is not None:
        num_data = min(num_data, MAX_TASKS)

    if SCALE_KERNELS:
        scale_factors = build_kernel_scale_factors(array_data)
    else:
        scale_factors = np.ones(num_data, dtype=np.float64)

    num_sens_sum = num_x1 * num_z1
    array_sens_sum = np.zeros((num_sens_sum, 3), dtype=np.float64)
    last_array = None

    Path(SENSITIVITY_IN_FILE).parent.mkdir(parents=True, exist_ok=True)

    with open(SENSITIVITY_IN_FILE, "w", encoding="utf-8") as file:
        for i in range(num_data):
            id_ns = int(array_data[i, 0])
            id_nr = int(array_data[i, 1])

            file_ker_integral = KERNEL_INTEGRAL.format(Ns=id_ns, Nr=id_nr)
            array_ker_integral = read_array(file_ker_integral, 3)

            if array_ker_integral.shape[0] != num_sens_sum:
                raise ValueError(
                    f"{file_ker_integral}: 点数 {array_ker_integral.shape[0]} "
                    f"与反演网格点数 {num_sens_sum} 不一致"
                )

            scale_factor = scale_factors[i]

            array_ker_scaled = array_ker_integral.copy()
            array_ker_scaled[:, 2] *= scale_factor

            array_sens_sum[:, 2] += array_ker_scaled[:, 2]

            np.savetxt(file, array_ker_scaled[:, 2], fmt="%.15e")

            if SAVE_SCALED_KERNELS:
                file_scaled = KERNEL_INTEGRAL_SCALED.format(Ns=id_ns, Nr=id_nr)
                Path(file_scaled).parent.mkdir(parents=True, exist_ok=True)
                np.savetxt(file_scaled, array_ker_scaled, fmt=OUTPUT_FMT)

            last_array = array_ker_scaled

    if last_array is None:
        raise ValueError("没有可汇总的积分核文件")

    array_sens_sum[:, 0] = last_array[:, 0]
    array_sens_sum[:, 1] = last_array[:, 1]

    Path(SENSITIVITY_FILE).parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(SENSITIVITY_FILE, array_sens_sum, fmt=OUTPUT_FMT)

def main():
    t0 = time.time()

    parameters = read_inversion_config(CONFIG_FILE)
    array_data = read_array(ERROR_FILE, 4)
    tasks = build_tasks(array_data)

    num_x1 = parameters.nx
    num_z1 = parameters.nz

    x_min1 = parameters.x1
    x_max1 = parameters.x2
    z_min1 = parameters.z1
    z_max1 = parameters.z2

    dx1 = (x_max1 - x_min1) / (num_x1 - 1)
    dz1 = (z_max1 - z_min1) / (num_z1 - 1)

    n_workers = parameters.num_cores if N_WORKERS is None else N_WORKERS
    n_workers = max(1, min(int(n_workers), os.cpu_count() or 1, len(tasks)))

    # print(f"Total tasks : {len(tasks)}")
    # print(f"Workers     : {n_workers}")

    if n_workers == 1:
        init_worker(
            x_min1,
            z_min1,
            dx1,
            dz1,
            num_x1,
            num_z1,
        )

        for task in tasks:
            process_one_kernel(task)

    else:
        with ProcessPoolExecutor(
            max_workers=n_workers,
            initializer=init_worker,
            initargs=(
                x_min1,
                z_min1,
                dx1,
                dz1,
                num_x1,
                num_z1,
            ),
        ) as executor:

            future_to_task = {
                executor.submit(process_one_kernel, task): task
                for task in tasks
            }

            done = 0

            for future in as_completed(future_to_task):
                done += 1
                task = future_to_task[future]

                try:
                    task_id, id_ns, id_nr, output_file = future.result()
                    # print(
                    #     f"[{done}/{len(tasks)}] "
                    #     f"Ns={id_ns}, Nr={id_nr}"
                    # )
                except Exception as exc:
                    print(
                        f"[{done}/{len(tasks)}] "
                        f"task_id={task[0]} failed: {exc}"
                    )
                    raise

    # 汇总必须等所有积分核生成完以后再执行。
    sum_integral_kernels(
        array_data,
        num_x1,
        num_z1,
    )

    # print(f"Elapsed time: {time.time() - t0:.2f} s")


if __name__ == "__main__":
    main()
