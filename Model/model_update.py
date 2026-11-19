#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import time
import read_file

from pathlib import Path
from scipy.interpolate import RegularGridInterpolator

config_file = Path("../../Input/parameters_inversion3d.config")
config = read_file.read_inversion_config3d(config_file)

# print(config)
# time.sleep(1000)

# ======================== parameters ========================

file_per = Path("vel_per3d1.txt")
file_ini = Path("vel_ini3d1.txt")
file_out = Path("vel_new3d.txt")

# True:
#     perturbation = dVp / Vp
#     Vp_new = Vp_old * (1 + perturbation)
#
# False:
#     perturbation = dVp，单位 m/s
#     Vp_new = Vp_old + perturbation
relative_per = True

output_fmt = "%.6f %.6f %.6f %.6f %.6f %.6f"

# ====================================================
# def build_perturbation_grid(data):
#     """
#     将 x y z perturbation 数据整理成规则三维网格。
#     """
#     x = np.unique(data[:, 0])
#     y = np.unique(data[:, 1])
#     z = np.unique(data[:, 2])

#     nx = x.size
#     ny = y.size
#     nz = z.size

#     expected_size = nx * ny * nz

#     if data.shape[0] != expected_size:
#         raise ValueError(
#             "perturbation file is not a complete regular 3-D grid: "
#             f"rows={data.shape[0]}, expected={expected_size}"
#         )

#     # 不依赖输入文件原始排序。
#     index_x = {value: i for i, value in enumerate(x)}
#     index_y = {value: i for i, value in enumerate(y)}
#     index_z = {value: i for i, value in enumerate(z)}

#     grid = np.empty((nx, ny, nz), dtype=np.float64)

#     for row in data:
#         ix = index_x[row[0]]
#         iy = index_y[row[1]]
#         iz = index_z[row[2]]

#         grid[ix, iy, iz] = row[3]

#     return x, y, z, grid


def build_perturbation_grid(data):
    """
    将 x y z perturbation 数据整理成规则三维网格。

    data:
        第1列 x
        第2列 y
        第3列 z
        第4列 perturbation
    """

    if data.ndim != 2 or data.shape[1] < 4:
        raise ValueError(
            "data must have at least 4 columns: x y z perturbation"
        )

    x, ix = np.unique(data[:, 0], return_inverse=True)
    y, iy = np.unique(data[:, 1], return_inverse=True)
    z, iz = np.unique(data[:, 2], return_inverse=True)

    nx, ny, nz = x.size, y.size, z.size
    expected_size = nx * ny * nz

    if data.shape[0] != expected_size:
        raise ValueError(
            "perturbation file is not a complete regular 3-D grid: "
            f"rows={data.shape[0]}, expected={expected_size}"
        )

    # 将三维索引转换为唯一的一维编号
    linear_index = np.ravel_multi_index(
        (ix, iy, iz),
        (nx, ny, nz)
    )

    # 检查是否存在重复 (x, y, z)
    if np.unique(linear_index).size != expected_size:
        raise ValueError(
            "perturbation file contains duplicate or missing "
            "(x, y, z) grid points"
        )

    grid = np.empty((nx, ny, nz), dtype=np.float64)

    grid[ix, iy, iz] = data[:, 3]

    return x, y, z, grid


def check_model_range(x, y, z, model):
    """
    检查输入模型是否完全位于反演扰动网格范围内。
    """
    xmin, xmax = model[:, 0].min(), model[:, 0].max()
    ymin, ymax = model[:, 1].min(), model[:, 1].max()
    zmin, zmax = model[:, 2].min(), model[:, 2].max()

    if xmin < x[0] or xmax > x[-1]:
        raise ValueError(
            f"initial-model x range [{xmin}, {xmax}] is outside "
            f"perturbation range [{x[0]}, {x[-1]}]"
        )

    if ymin < y[0] or ymax > y[-1]:
        raise ValueError(
            f"initial-model y range [{ymin}, {ymax}] is outside "
            f"perturbation range [{y[0]}, {y[-1]}]"
        )

    if zmin < z[0] or zmax > z[-1]:
        raise ValueError(
            f"initial-model z range [{zmin}, {zmax}] is outside "
            f"perturbation range [{z[0]}, {z[-1]}]"
        )


def interpolate_perturbation(data_per, model):
    """
    将粗反演网格上的扰动插值到输入模型网格。
    """
    x, y, z, grid = build_perturbation_grid(data_per)

    check_model_range(x, y, z, model)

    print("Perturbation grid:")
    print(f"  nx = {x.size}")
    print(f"  ny = {y.size}")
    print(f"  nz = {z.size}")

    print("Perturbation range:")
    print(f"  x = {x[0]:.6f} ~ {x[-1]:.6f}")
    print(f"  y = {y[0]:.6f} ~ {y[-1]:.6f}")
    print(f"  z = {z[0]:.6f} ~ {z[-1]:.6f}")

    interpolator = RegularGridInterpolator(
        (x, y, z),
        grid,
        method="linear",
        bounds_error=True,
    )

    target_points = model[:, 0:3]

    perturbation_interp = interpolator(target_points)

    return perturbation_interp


def update_model(model, perturbation):
    """
    根据插值后的扰动更新 Vp。

    Vs 和 rho 当前保持不变。
    """
    output = model.copy()

    if relative_per:
        output[:, 3] = model[:, 3] * (1.0 + perturbation)
    else:
        output[:, 3] = model[:, 3] + perturbation

    return output


def save_model(filename, data):
    filename = Path(filename)

    filename.parent.mkdir(parents=True, exist_ok=True)

    np.savetxt(
        filename,
        data,
        fmt=output_fmt,
    )


def main():
    # 1. 读取反演扰动
    data_per = read_file.read_perturbation(file_per)

    # 2. 读取输入模型
    model_ini = read_file.read_initial_model(file_ini)

    # 3. 将粗网格扰动插值到输入模型网格
    perturbation_interp = interpolate_perturbation(
        data_per,
        model_ini,
    )

    # 4. 更新模型
    model_new = update_model(
        model_ini,
        perturbation_interp,
    )

    # 5. 输出
    save_model(
        file_out,
        model_new,
    )

    print("")
    print(f"Initial model : {file_ini}")
    print(f"Perturbation  : {file_per}")
    print(f"Output model  : {file_out}")

    print("")
    print(
        f"Perturbation min/max = "
        f"{perturbation_interp.min():.6e} / "
        f"{perturbation_interp.max():.6e}"
    )

    print(
        f"Vp old min/max       = "
        f"{model_ini[:, 3].min():.6f} / "
        f"{model_ini[:, 3].max():.6f}"
    )

    print(
        f"Vp new min/max       = "
        f"{model_new[:, 3].min():.6f} / "
        f"{model_new[:, 3].max():.6f}"
    )


if __name__ == "__main__":
    main()