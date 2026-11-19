#!/usr/bin/env python3
"""二维核局部积分 → 沿测线映射到三维 → 投影到三维反演网格。"""

from itertools import product
from pathlib import Path

import numpy as np
from scipy.signal import convolve2d

from read_file import read_inversion_config3d


CONFIG_FILE = Path("../../Input/parameters_inversion3d.config")
ERROR_FILE = Path("../../Output/Error/Error.txt")
LINE_FILE = Path("../../Output/Model/profile2d/line.txt")
SMOOTH_TEMPLATE = "../../Output/Kernel/Kernel_soomth_s{ns}_r{nr}.dat"
INTEGRAL_2D_TEMPLATE = "../../Output/Kernel/Kernel_integral2d_s{ns}_r{nr}.dat"
SCALED_2D_TEMPLATE = "../../Output/Kernel/Kernel_integral2d_scaled_s{ns}_r{nr}.dat"
PROFILE_3D_TEMPLATE = "../../Output/Kernel/Kernel_profile3d_s{ns}_r{nr}.dat"
INVERSION_3D_TEMPLATE = "../../Output/Kernel/Kernel_inversion3d_s{ns}_r{nr}.dat"
SCALE_INFO_FILE = Path("../../Output/Kernel/Kernel_scale_info_3d.dat")
OUTPUT_FMT = "%.15e"

# 按各条积分核的 max(|K|) 做跨测线幅值均衡；不改变单条核的空间形状。
SCALE_REFERENCE = "median"  # median 或 mean
SCALE_POWER = 1.0
MIN_SCALE_FACTOR = 0.01
MAX_SCALE_FACTOR = 100.0


def read_array(filename, columns):
    data = np.loadtxt(filename, ndmin=2)
    if data.shape[0] == 0 or data.shape[1] < columns or not np.isfinite(data).all():
        raise ValueError(f"{filename}: 需要至少 {columns} 列有限数值")
    return data


def read_lines(filename):
    """九列：ns nr xs ys zs xe ye ze line_length（单位 m）。"""
    rows = read_array(filename, 9)
    if rows.shape[1] != 9:
        raise ValueError(f"{filename}: line.txt 必须恰好为九列")
    lines = {}
    for row in rows:
        ns, nr = row[:2]
        if ns < 1 or nr < 1 or ns != int(ns) or nr != int(nr) or row[8] <= 0:
            raise ValueError(f"{filename}: 无效编号或测线长度：{row}")
        pair = (int(ns), int(nr))
        if pair in lines:
            raise ValueError(f"{filename}: 重复的震源–台站对 {pair}")
        lines[pair] = (row[2:5], row[5:8], row[8])
    return lines


def integrate_profile(kernel):
    """输入按 x 外层、z 内层排列；z 为负值，向地表递增。"""
    x, z = np.unique(kernel[:, 0]), np.unique(kernel[:, 1])
    nx, nz = len(x), len(z)
    if nx < 2 or nz < 2 or nx * nz != len(kernel):
        raise ValueError("二维核必须是至少 2×2 的完整规则网格")
    if not (np.allclose(kernel[:, 0].reshape(nx, nz), x[:, None])
            and np.allclose(kernel[:, 1].reshape(nx, nz), z[None, :])):
        raise ValueError("二维核需要按 x 外层、z 内层递增的顺序排列")

    dx, dz = x[1] - x[0], z[1] - z[0]
    if not (np.allclose(np.diff(x), dx) and np.allclose(np.diff(z), dz)):
        raise ValueError("二维核坐标必须等间距")

    weights = np.outer([1.0, 4.0, 1.0], [1.0, 4.0, 1.0]) * (dx * dz / 36.0)
    values = convolve2d(kernel[:, 2].reshape(nx, nz), weights,
                        mode="same", boundary="fill", fillvalue=0.0)
    return np.column_stack((kernel[:, :2], values.ravel()))


def profile_to_3d(integral_2d, line):
    start, end, length = line
    t = integral_2d[:, 0] / length
    xyz = start + t[:, None] * (end - start)
    xyz[:, 2] -= integral_2d[:, 1]  # 负的二维 z → 正的三维深度
    return np.column_stack((xyz, integral_2d[:, 2]))


def inversion_grid(config):
    shape = (config.nx_inv, config.ny_inv, config.nz_inv)
    lower = np.array([config.x1, config.y1, config.z1], dtype=float)
    upper = np.array([config.x2, config.y2, config.z2], dtype=float)
    step = (upper - lower) / (np.array(shape) - 1)
    axes = [np.linspace(a, b, n) for a, b, n in zip(lower, upper, shape)]
    xyz = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)
    return np.column_stack((xyz, np.zeros(len(xyz)))), lower, upper, step, shape


def project_to_inversion_grid(profile_3d, grid, lower, upper, step, shape):
    """与 project2invergrid 相同的累加方式，扩展为周围八节点权重。"""
    result = grid.copy()
    xyz, values = profile_3d[:, :3], profile_3d[:, 3]
    tolerance = 1e-10 * np.max(upper - lower)
    inside = np.all((xyz >= lower - tolerance) & (xyz <= upper + tolerance), axis=1)
    fraction = np.clip((xyz[inside] - lower) / step, 0, np.array(shape) - 1)
    base = np.floor(fraction).astype(np.int64)
    part = fraction - base

    for ox, oy, oz in product((0, 1), repeat=3):
        offsets = np.array([ox, oy, oz])
        indices = base + offsets
        valid = np.all(indices < shape, axis=1)
        weight = np.prod(np.where(offsets == 1, part, 1 - part), axis=1)
        flat = np.ravel_multi_index(indices[valid].T, shape)
        np.add.at(result[:, 3], flat, values[inside][valid] * weight[valid])

    return result, int(np.count_nonzero(~inside))


def save(filename, array):
    filename = Path(filename)
    filename.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(filename, array, fmt=OUTPUT_FMT)


def kernel_scale_factors(maxima):
    positive = maxima > 0
    if not np.any(positive):
        raise ValueError("所有积分敏感核均为零，无法计算幅值均衡系数")
    if SCALE_REFERENCE == "median":
        reference = np.median(maxima[positive])
    elif SCALE_REFERENCE == "mean":
        reference = np.mean(maxima[positive])
    else:
        raise ValueError("SCALE_REFERENCE 必须是 median 或 mean")
    factors = np.ones(len(maxima))
    factors[positive] = np.clip((reference / maxima[positive]) ** SCALE_POWER,
                                MIN_SCALE_FACTOR, MAX_SCALE_FACTOR)
    return reference, factors


def main():
    config = read_inversion_config3d(CONFIG_FILE)
    lines = read_lines(LINE_FILE)
    errors = read_array(ERROR_FILE, 4)
    grid, lower, upper, step, shape = inversion_grid(config)
    seen = set()

    pairs = []
    maxima = []
    # 第一遍：逐对读取平滑核，计算并保存未经幅值均衡的二维积分核。
    for row in errors:
        pair = (int(row[0]), int(row[1]))
        if pair in seen or pair not in lines or not np.array_equal(row[:2], pair):
            raise ValueError(f"Error.txt 中的编号重复、无效或 line.txt 缺少该测线：{pair}")
        seen.add(pair)
        pairs.append(pair)
        ns, nr = pair
        kernel = read_array(SMOOTH_TEMPLATE.format(ns=ns, nr=nr), 3)[:, :3]
        integral_2d = integrate_profile(kernel)
        maxima.append(np.max(np.abs(integral_2d[:, 2])))
        save(INTEGRAL_2D_TEMPLATE.format(ns=ns, nr=nr), integral_2d)

    reference, factors = kernel_scale_factors(np.asarray(maxima))
    info = np.column_stack((np.asarray(pairs), maxima, factors,
                            np.asarray(maxima) * factors))
    SCALE_INFO_FILE.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(SCALE_INFO_FILE, info, fmt="%d %d %.15e %.15e %.15e")
    print(f"积分核幅值参考值：{reference:.15e}")

    # 第二遍：按整条核的系数缩放，再映射三维并投影到三维反演网格。
    for (ns, nr), factor in zip(pairs, factors):
        integral_2d = read_array(INTEGRAL_2D_TEMPLATE.format(ns=ns, nr=nr), 3)
        integral_2d[:, 2] *= factor
        save(SCALED_2D_TEMPLATE.format(ns=ns, nr=nr), integral_2d)
        profile_3d = profile_to_3d(integral_2d, lines[(ns, nr)])
        projected, outside = project_to_inversion_grid(
            profile_3d, grid, lower, upper, step, shape)

        save(PROFILE_3D_TEMPLATE.format(ns=ns, nr=nr), profile_3d)
        save(INVERSION_3D_TEMPLATE.format(ns=ns, nr=nr), projected)
        print(f"s{ns}_r{nr}: 缩放系数={factor:.6g}，{outside} 个点位于反演网格外")


if __name__ == "__main__":
    main()
