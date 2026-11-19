#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""根据 line.txt, 将垂直剖面的 x2 z vp vs rho 转为 x y z vp vs rho。"""

from pathlib import Path
import numpy as np
import time

# ======================== 用户参数 ========================
# 相对路径以运行命令时所在的目录为基准；所有长度使用相同单位。
line_file = Path("../../Output/Model/profile2d/line.txt")
input_dir = Path("../../Output/Model/profile2d/")
output_dir = Path("../../Output/Model/profile3d/")
input_pattern = "profile2d_s{ns}_r{nr}.txt"
output_pattern = "profile3d_s{ns}_r{nr}.txt"

# 长度检查允许 line.txt 的小数舍入误差。
coord_rtol = 1.0e-7
coord_atol = 1.0e-6
# 17位有效数字, 避免额外降低原模型参数的数值精度。
output_fmt = "%.6f"
# =========================================================

def read_lines(filename):
    """读取 Ns Nr xs ys zs xe ye ze L, 并检查编号和长度。"""
    lines = np.loadtxt(filename, dtype=np.float64, ndmin=2)
    if lines.shape[0] == 0 or lines.shape[1] != 9:
        raise ValueError(f"{filename}: expected 9 columns and at least one row")
    if not np.all(np.isfinite(lines)):
        raise ValueError(f"{filename}: contains NaN or Inf")
    if np.any(lines[:, :2] != np.rint(lines[:, :2])):
        raise ValueError("Ns and Nr must be integers")
    if len(np.unique(lines[:, :2], axis=0)) != len(lines):
        raise ValueError("Duplicate (Ns, Nr) in line file")

    lengths = np.hypot(lines[:, 5] - lines[:, 2], lines[:, 6] - lines[:, 3])
    valid = (lines[:, 8] > 0) & (lengths > 0)
    valid &= np.isclose(lines[:, 8], lengths, rtol=coord_rtol, atol=coord_atol)
    if not np.all(valid):
        rows = (np.flatnonzero(~valid) + 1).tolist()
        raise ValueError(f"Invalid horizontal line length in rows: {rows}")
    return lines


def convert_profile(data, line_row):
    """z是原三维模型的绝对z坐标, 不是相对起点深度。"""
    if data.ndim != 2 or data.shape[0] == 0 or data.shape[1] != 5:
        raise ValueError("2D model must contain exactly 5 columns: x2 z vp vs rho")
    if not np.all(np.isfinite(data)):
        raise ValueError("2D model contains NaN or Inf")

    xs, ys = line_row[2:4]
    xe, ye = line_row[5:7]
    # length = line_row[8]

    # 与提取二维切片时采用相同的长度计算方法。
    length = float(np.hypot(xe - xs, ye - ys))
    if length <= 0.0:
        raise ValueError("Horizontal line length must be greater than 0")
    # 最后一列长度用于一致性检查。
    if not np.isclose(length, line_row[8], rtol=coord_rtol, atol=coord_atol):
        raise ValueError(
            f"Line length mismatch: calculated={length}, recorded={line_row[8]}"
    )

    distance = data[:, 0]

    tolerance = coord_atol + coord_rtol * length
    if np.any(distance < -tolerance) or np.any(distance > length + tolerance):
        raise ValueError(f"x2 must lie within [0, {length}]")

    # 按每个点实际的x2转换，允许非等间距；保持输入顺序和深度。
    fraction = distance / length

    # print(data[:10, :2])
    # print(distance)
    # print(length)
    # print(fraction)
    # time.sleep(1000)
    x = xs + fraction * (xe - xs)
    y = ys + fraction * (ye - ys)
    return np.column_stack((x, y, data[:, 1:]))


def main():
    lines = read_lines(line_file)
    output_dir.mkdir(parents=True, exist_ok=True)
    converted = skipped = 0

    for row in lines:
        ns, nr = int(row[0]), int(row[1])
        source = input_dir / input_pattern.format(ns=ns, nr=nr)
        target = output_dir / output_pattern.format(ns=ns, nr=nr)

        if not source.is_file():
            print(f"Skip missing file: {source}")
            skipped += 1
            continue
        if source.resolve() == target.resolve():
            raise ValueError(f"Input and output must differ: {source}")

        data = np.loadtxt(source, dtype=np.float64, ndmin=2)

        try:
            result = convert_profile(data, row)
        except ValueError as error:
            raise ValueError(f"{source}: {error}") from error

        # 覆盖已有输出文件，写入本次转换结果
        with target.open("w", encoding="utf-8") as stream:
            np.savetxt(stream, result, fmt=output_fmt)
        converted += 1
        # print(f"Saved: {target} ({len(result)} points)")

    print(f"profile_2d_to_3d done: converted={converted}, missing={skipped}")

if __name__ == "__main__":
    main()
