#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import time
from pathlib import Path
from read_file import read_inversion_config3d
from scipy.interpolate import RegularGridInterpolator

# ======================== 用户参数 ========================

# MODEL_FILE = Path("../../Output/Model/model_2d3d/model3d_ini.txt")

config_file = Path("../../Input/parameters_inversion3d.config")
model_file = Path("../../Output/Model/Forver/vel3d.txt")
line_file = Path("../../Output/Model/profile2d/line.txt")

output_dir_3d = Path("../../Output/Model/profile3d/")
output_dir_2d = Path("../../Output/Model/profile2d/")

# 三维模型预期参数，用于检查 model3d_ini.txt
config = read_inversion_config3d(config_file)
model_x_start = config.x1
model_x_end = config.x2
model_y_start = config.y1
model_y_end = config.y2
model_z_start = config.z1
model_z_end = config.z2

# print(MODEL_X_START, MODEL_X_END, MODEL_Y_START, MODEL_Y_END, MODEL_Z_START, MODEL_Z_END)
# time.sleep(1000)

numNx = config.nx_for
numNy = config.ny_for
numNz = config.nz_for

# print(numNx, numNy, numNz)
# time.sleep(1000)

model_dx = (model_x_end - model_x_start) / (config.nx_for - 1)
model_dy = (model_y_end - model_y_start) / (config.ny_for - 1)
model_dz = (model_z_end - model_z_start) / (config.nz_for - 1)


# print(model_dx, model_dy, model_dz)
# time.sleep(1000)

# 二维剖面沿测线方向采样间隔
profile_dx = 1000.0

# 若测线总长不能被 PROFILE_DX 整除，是否强制加入测线终点
include_line_end = True

coord_rtol = 1.0e-7
coord_atol = 1.0e-6

output_fmt_3d = "%.6f %.6f %.6f %.2f %.2f %.2f"
output_fmt_2d = "%.6f %.6f %.2f %.2f %.2f"

# =========================================================

def read_regular_model(filename):
    """读取规则三维模型：x y z vp vs rho。"""
    filename = Path(filename)

    if not filename.is_file():
        raise FileNotFoundError(f"cannot open model file: {filename}")

    data = np.loadtxt(filename, dtype=np.float64, ndmin=2)

    if data.shape[1] < 6:
        raise ValueError(f"{filename} must contain at least 6 columns")

    if not np.all(np.isfinite(data[:, :6])):
        raise ValueError(f"{filename} contains NaN or Inf")

    x = np.unique(data[:, 0])
    y = np.unique(data[:, 1])
    z = np.unique(data[:, 2])

    nx, ny, nz = len(x), len(y), len(z)

    if data.shape[0] != nx*ny*nz:
        raise ValueError(
            f"model is not a complete regular grid: "
            f"rows={data.shape[0]}, nx*ny*nz={nx * ny * nz}"
        )

    ix = np.searchsorted(x, data[:, 0])
    iy = np.searchsorted(y, data[:, 1])
    iz = np.searchsorted(z, data[:, 2])

    vp = np.empty((nx, ny, nz), dtype=np.float64)
    vs = np.empty((nx, ny, nz), dtype=np.float64)
    rho = np.empty((nx, ny, nz), dtype=np.float64)

    vp[ix, iy, iz] = data[:, 3]
    vs[ix, iy, iz] = data[:, 4]
    rho[ix, iy, iz] = data[:, 5]

    return x, y, z, vp, vs, rho


def check_uniform_axis(axis, expected_start, expected_end, expected_d, name):
    """检查规则网格某一坐标轴。"""
    if axis.size < 2:
        raise ValueError(f"{name} direction must contain at least 2 nodes")

    if not np.isclose(axis[0], expected_start, rtol=coord_rtol, atol=coord_atol):
        raise ValueError(
            f"{name} start mismatch: file={axis[0]:.15e}, "
            f"expected={expected_start:.15e}"
        )

    if not np.isclose(axis[-1], expected_end, rtol=coord_rtol, atol=coord_atol):
        raise ValueError(
            f"{name} end mismatch: file={axis[-1]:.15e}, "
            f"expected={expected_end:.15e}"
        )

    spacing = np.diff(axis)

    if not np.allclose(spacing, expected_d, rtol=coord_rtol, atol=coord_atol):
        raise ValueError(
            f"{name} spacing mismatch: "
            f"file min/max=({spacing.min():.15e}, {spacing.max():.15e}), "
            f"expected={expected_d:.15e}"
        )

    expected_n = int(round((expected_end - expected_start) / expected_d)) + 1

    if axis.size != expected_n:
        raise ValueError(
            f"{name} node count mismatch: file={axis.size}, expected={expected_n}"
        )


def check_model_parameters(x, y, z):
    """检查三维模型参数是否与用户定义一致。"""
    check_uniform_axis(x, model_x_start, model_x_end, model_dx, "X")
    check_uniform_axis(y, model_y_start, model_y_end, model_dy, "Y")
    check_uniform_axis(z, model_z_start, model_z_end, model_dz, "Z")

    print("3D model parameter check: OK")
    print(f"Model grid: nx={len(x)}, ny={len(y)}, nz={len(z)}")


def read_lines(filename):
    """读取 line.txt：Ns Nr xs ys zs xe ye ze。"""
    filename = Path(filename)

    if not filename.is_file():
        raise FileNotFoundError(f"cannot open line file: {filename}")

    data = np.loadtxt(filename, dtype=np.float64, ndmin=2)

    if data.shape[1] < 8:
        raise ValueError(f"{filename} must contain at least 8 columns")

    if not np.all(np.isfinite(data[:, :8])):
        raise ValueError(f"{filename} contains NaN or Inf")

    return data[:, :8]


def build_profile_distance(line_length, spacing):
    """建立距离测线起点的二维局部坐标 x2。"""
    if spacing <= 0.0:
        raise ValueError("PROFILE_DX must be greater than 0")

    if line_length <= 0.0:
        raise ValueError("line length must be greater than 0")

    n = int(np.floor(line_length / spacing))
    distance = np.arange(n + 1, dtype=np.float64) * spacing

    if include_line_end:
        if not np.isclose(
            distance[-1],
            line_length,
            rtol=coord_rtol,
            atol=coord_atol,
        ):
            distance = np.append(distance, line_length)
        else:
            distance[-1] = line_length

    return distance


def extract_vertical_profile(line_row, z_axis, vp_interp, vs_interp, rho_interp):
    """
    提取沿测线垂直向下的二维剖面。

    返回：
        profile_3d: x1 y1 z1 vp vs rho
        profile_2d: x2 z1 vp vs rho
    """
    ns = int(round(line_row[0]))
    nr = int(round(line_row[1]))

    x_start, y_start, z_start = line_row[2:5]
    x_end, y_end, z_end = line_row[5:8]

    dx = x_end - x_start
    dy = y_end - y_start

    # 垂直剖面使用水平测线方向；
    # 深度统一采用三维模型原始 z_axis。
    line_length = float(np.hypot(dx, dy))


    if line_length <= 0.0:
        raise ValueError(
            f"Ns={ns}, Nr={nr}: horizontal line length is zero"
        )

    distance = build_profile_distance(line_length, profile_dx)

    ux = dx / line_length
    uy = dy / line_length

    x_line = x_start + distance * ux
    y_line = y_start + distance * uy

    nd = distance.size
    nz = z_axis.size

    xx = np.repeat(x_line, nz)
    yy = np.repeat(y_line, nz)
    zz = np.tile(z_axis, nd)
    dd = np.repeat(distance, nz)

    points = np.column_stack((xx, yy, zz))

    vp = vp_interp(points)
    vs = vs_interp(points)
    rho = rho_interp(points)
    # print("111")
    # time.sleep(1000)
    profile_3d = np.column_stack((xx, yy, zz, vp, vs, rho))
    profile_2d = np.column_stack((dd, zz, vp, vs, rho))

    return ns, nr, profile_3d, profile_2d


def save_profile_3d(ns, nr, profile_3d):
    """保存 x1 y1 z1 vp vs rho。"""
    output_dir_3d.mkdir(parents=True, exist_ok=True)

    filename = output_dir_3d / f"profile_s{ns}_r{nr}.txt"
    np.savetxt(filename, profile_3d, fmt=output_fmt_3d)

    return filename


def save_profile_2d(ns, nr, profile_2d):
    """保存 x2 z1 vp vs rho，其中 x2 为距测线起点的距离。"""
    output_dir_2d.mkdir(parents=True, exist_ok=True)

    filename = output_dir_2d / f"profile2d_s{ns}_r{nr}.txt"
    np.savetxt(filename, profile_2d, fmt=output_fmt_2d)

    return filename


def main():
    print("\n--------------------------------")
    
    print("Extracting 2D slices from 3D")
    # 1. 读取三维初始模型
    x, y, z, vp, vs, rho = read_regular_model(model_file)

    # 2. 检查三维模型范围和网格参数
    check_model_parameters(x, y, z)

    # 3. 建立三线性插值器
    vp_interp = RegularGridInterpolator(
        (x, y, z), vp, method="linear", bounds_error=True
    )
    vs_interp = RegularGridInterpolator(
        (x, y, z), vs, method="linear", bounds_error=True
    )
    rho_interp = RegularGridInterpolator(
        (x, y, z), rho, method="linear", bounds_error=True
    )

    # 4. 逐条测线提取二维垂直剖面
    line_data = read_lines(line_file)

    for row in line_data:
        ns, nr, profile_3d, profile_2d = extract_vertical_profile(
            row,
            z,
            vp_interp,
            vs_interp,
            rho_interp,
        )

        # file_3d = save_profile_3d(ns, nr, profile_3d)
        file_2d = save_profile_2d(ns, nr, profile_2d)

    print("--------------------------------\n")

if __name__ == "__main__":
    main()
