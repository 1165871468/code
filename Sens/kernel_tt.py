#!/usr/bin/env python3
"""汇总三维反演网格上的走时积分敏感核。"""

import numpy as np
from pathlib import Path
from read_file import read_inversion_config3d


CONFIG_FILE = Path("../../Input/parameters_inversion3d.config")
ERROR_FILE = Path("../../Output/Error/Error.txt")
KERNEL_TEMPLATE = "../../Output/Kernel/Kernel_inversion3d_s{ns}_r{nr}.dat"
SENSITIVITY_IN_FILE = Path("../../Output/Model/Inver/Sensitivity_in_tt.dat")
SENSITIVITY_SUM_FILE = Path("../../Output/Model/Inver/Sensitivity_sum_tt.dat")
OUTPUT_FMT = "%.15e"


def inversion_coordinates(config):
    """节点顺序：x 外层、y 中层、z 内层。"""
    axes = (
        np.linspace(config.x1, config.x2, config.nx_inv),
        np.linspace(config.y1, config.y2, config.ny_inv),
        np.linspace(config.z1, config.z2, config.nz_inv),
    )
    return np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)


def read_pairs(path):
    data = np.loadtxt(path, ndmin=2)
    if data.shape[1] < 2 or not np.isfinite(data).all():
        raise ValueError(f"{path}: 需要至少两列有限数值（震源编号、台站编号）")
    ids = data[:, :2]
    if np.any(ids < 1) or np.any(ids != np.floor(ids)):
        raise ValueError(f"{path}: 震源和台站编号必须为正整数")
    pairs = [tuple(row) for row in ids.astype(int)]
    if len(pairs) != len(set(pairs)):
        raise ValueError(f"{path}: 存在重复的震源—台站对")
    return pairs


def read_kernel(path, coordinates):
    kernel = np.loadtxt(path, ndmin=2)
    if kernel.shape != (len(coordinates), 4) or not np.isfinite(kernel).all():
        raise ValueError(f"{path}: 必须为 {len(coordinates)} 行、4 列有限数值 x y z K")
    if not np.allclose(kernel[:, :3], coordinates, rtol=1e-10, atol=1e-7):
        raise ValueError(f"{path}: 三维网格坐标或行顺序与反演参数不一致")
    return kernel[:, 3]


def main():
    config = read_inversion_config3d(CONFIG_FILE)
    if min(config.nx_inv, config.ny_inv, config.nz_inv) < 2:
        raise ValueError("三维反演网格各方向至少需要 2 个节点")
    coordinates = inversion_coordinates(config)
    pairs = read_pairs(ERROR_FILE)
    total = np.zeros(len(coordinates), dtype=float)

    SENSITIVITY_IN_FILE.parent.mkdir(parents=True, exist_ok=True)
    SENSITIVITY_SUM_FILE.parent.mkdir(parents=True, exist_ok=True)
    # 每个震源—台站对占一段连续的 N 个标量，顺序与 Error.txt 一致
    with SENSITIVITY_IN_FILE.open("w", encoding="utf-8") as stream:
        for ns, nr in pairs:
            path = Path(KERNEL_TEMPLATE.format(ns=ns, nr=nr))
            values = read_kernel(path, coordinates)
            np.savetxt(stream, values, fmt=OUTPUT_FMT)
            total += values

    np.savetxt(SENSITIVITY_SUM_FILE, np.column_stack((coordinates, total)), fmt=OUTPUT_FMT)
    print(f"已汇总 {len(pairs)} 条走时积分敏感核，每条 {len(coordinates)} 个网格点")
    print(SENSITIVITY_IN_FILE)
    print(SENSITIVITY_SUM_FILE)


if __name__ == "__main__":
    main()
