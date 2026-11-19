#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
根据 integral_parallel_scaled.py 已生成的单个震源-台站积分敏感核，
直接构造共震源双差敏感核。

前置程序：
    python3 integral_parallel_scaled.py

输入：
    ../../Output/Error/Error_dd_cs.txt
        Ns Nr1 Nr2 error_dd

    ../../Output/Kernel/Kernel_integral_s{Ns}_r{Nr}.dat
        x z K

输出：
    ../../Output/Model/Inver/Sensitivity_in_dd_cs.dat
        每个共震源双差核的第三列连续写入，用于反演矩阵

    ../../Output/Model/Inver/Sensitivity_dd_cs.dat
        所有共震源双差核的累加结果

可选输出：
    ../../Output/Kernel/Kernel_dd_cs_s{Ns}_r{Nr1}_r{Nr2}.dat
"""

from pathlib import Path
import numpy as np


# ======================== 用户参数 ========================

ERROR_DD_FILE = Path("../../Output/Error/Error_dd_cs.txt")

# 推荐使用未放缩的单核构造双差：
#     Kdd = K(s,r1) - K(s,r2)
KERNEL_INTEGRAL = "../../Output/Kernel/Kernel_integral_s{Ns}_r{Nr}.dat"

# 如果确实希望直接使用 integral_parallel_scaled.py 已放缩的单核，
# 把 USE_SCALED_SINGLE_KERNELS 改为 True。
KERNEL_INTEGRAL_SCALED = "../../Output/Kernel/Kernel_integral_scaled_s{Ns}_r{Nr}.dat"
USE_SCALED_SINGLE_KERNELS = True

SENSITIVITY_IN_FILE = Path("../../Output/Model/Inver/Sensitivity_in_dd_cs.dat")
SENSITIVITY_SUM_FILE = Path("../../Output/Model/Inver/Sensitivity_dd_cs.dat")

SAVE_DD_KERNELS = False
KERNEL_DD_TEMPLATE = "../../Output/Kernel/Kernel_dd_cs_s{Ns}_r{Nr1}_r{Nr2}.dat"

OUTPUT_FMT = "%.15e %.15e %.15e"

# None：全部处理；调试时可以设为 10。
MAX_TASKS = None

# =========================================================


def read_array(filename, min_cols):
    filename = Path(filename)

    if not filename.is_file():
        raise FileNotFoundError(f"cannot open this file -- {filename}")

    data = np.loadtxt(filename, dtype=np.float64, ndmin=2)

    if data.shape[1] < min_cols:
        raise ValueError(f"{filename} must contain at least {min_cols} columns")

    if not np.all(np.isfinite(data[:, :min_cols])):
        raise ValueError(f"{filename} contains NaN or Inf")

    return data


def get_single_kernel_file(ns, nr):
    template = KERNEL_INTEGRAL_SCALED if USE_SCALED_SINGLE_KERNELS else KERNEL_INTEGRAL
    return Path(template.format(Ns=ns, Nr=nr))


def build_common_source_dd_kernel(ns, nr1, nr2):
    """
    共震源双差敏感核：
        Kdd(s,r1,r2) = K(s,r1) - K(s,r2)
    """
    file1 = get_single_kernel_file(ns, nr1)
    file2 = get_single_kernel_file(ns, nr2)

    ker1 = read_array(file1, 3)
    ker2 = read_array(file2, 3)

    if ker1.shape != ker2.shape:
        raise ValueError(
            f"kernel shape mismatch:\n"
            f"  {file1}: {ker1.shape}\n"
            f"  {file2}: {ker2.shape}"
        )

    if not np.allclose(
        ker1[:, :2],
        ker2[:, :2],
        rtol=1.0e-10,
        atol=1.0e-12,
    ):
        raise ValueError(
            f"kernel coordinates are inconsistent:\n"
            f"  {file1}\n"
            f"  {file2}"
        )

    dd = ker1.copy()
    dd[:, 2] = ker1[:, 2] - ker2[:, 2]

    return dd


def main():
    error_dd = read_array(ERROR_DD_FILE, 4)

    if MAX_TASKS is not None:
        error_dd = error_dd[:MAX_TASKS]

    if error_dd.shape[0] == 0:
        raise ValueError("Error_dd_cs.txt contains no double-difference records")

    SENSITIVITY_IN_FILE.parent.mkdir(parents=True, exist_ok=True)
    SENSITIVITY_SUM_FILE.parent.mkdir(parents=True, exist_ok=True)

    sensitivity_sum = None
    num_grid = None

    with SENSITIVITY_IN_FILE.open("w", encoding="utf-8") as fout:
        for i, row in enumerate(error_dd):
            ns = int(round(row[0]))
            nr1 = int(round(row[1]))
            nr2 = int(round(row[2]))

            dd = build_common_source_dd_kernel(ns, nr1, nr2)

            if num_grid is None:
                num_grid = dd.shape[0]
                sensitivity_sum = np.zeros_like(dd)
                sensitivity_sum[:, :2] = dd[:, :2]
            elif dd.shape[0] != num_grid:
                raise ValueError(
                    f"Ns={ns}, Nr1={nr1}, Nr2={nr2}: "
                    f"grid point count {dd.shape[0]} != {num_grid}"
                )

            # 最终反演矩阵按 Error_dd_cs.txt 的原始顺序排列。
            np.savetxt(fout, dd[:, 2], fmt="%.15e")

            sensitivity_sum[:, 2] += dd[:, 2]

            if SAVE_DD_KERNELS:
                file_dd = Path(
                    KERNEL_DD_TEMPLATE.format(
                        Ns=ns,
                        Nr1=nr1,
                        Nr2=nr2,
                    )
                )
                file_dd.parent.mkdir(parents=True, exist_ok=True)
                np.savetxt(file_dd, dd, fmt=OUTPUT_FMT)

    np.savetxt(
        SENSITIVITY_SUM_FILE,
        sensitivity_sum,
        fmt=OUTPUT_FMT,
    )

    kernel_type = "scaled single kernels" if USE_SCALED_SINGLE_KERNELS else "original integral kernels"

    print(f"Common-source DD records : {error_dd.shape[0]}")
    print(f"Grid points per kernel   : {num_grid}")
    # print(f"Input kernel type        : {kernel_type}")
    # print(f"Sensitivity matrix       : {SENSITIVITY_IN_FILE}")
    # print(f"Summed DD kernel         : {SENSITIVITY_SUM_FILE}")


if __name__ == "__main__":
    main()
