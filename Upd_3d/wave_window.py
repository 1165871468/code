#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
对波形文件施加指定时间范围的 Tukey 窗。

输入文件格式：
    第1列：时间
    第2列及以后：波形分量

运行方式：
    python3 wave_window.py wave_s1_r1.txt 4.0 12.0

参数：
    第1个参数：输入波形文件名
    第2个参数：时间窗起始时间
    第3个参数：时间窗结束时间
"""

from pathlib import Path
import sys
import time
import numpy as np
from scipy.signal.windows import tukey


# ======================== 脚本参数 ========================

# Tukey窗参数：
# alpha=0.0：矩形窗
# alpha=1.0：Hann窗
# alpha=0.2：总窗长的20%用于两端渐变，即每端约10%
TUKEY_ALPHA = 1.0

# 输出文件名后缀。
OUTPUT_SUFFIX = "_win"

# 输出格式。
OUTPUT_FMT = "%.10e"

# 是否在最后增加一列窗函数值。
OUTPUT_WINDOW_COLUMN = False

# =========================================================


def read_waveform(file_name):
    file_path = Path(file_name).expanduser()

    if not file_path.is_file():
        raise FileNotFoundError(f"找不到波形文件：{file_path}")

    data = np.loadtxt(file_path, dtype=np.float64)

    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] < 2:
        raise ValueError("波形文件至少需要两列：时间列和一个波形分量。")

    if data.shape[0] < 2:
        raise ValueError("波形文件至少需要两个采样点。")

    if not np.all(np.isfinite(data)):
        raise ValueError("波形文件中存在 NaN 或 Inf。")

    if np.any(np.diff(data[:, 0]) <= 0.0):
        raise ValueError("时间列必须严格递增。")

    return file_path, data


def create_tukey_time_window(data_time, t_start, t_end, alpha):
    if t_end <= t_start:
        raise ValueError("结束时间必须大于起始时间。")

    if not 0.0 <= alpha <= 1.0:
        raise ValueError("TUKEY_ALPHA 必须位于 [0, 1]。")

    if t_start < data_time[0] or t_end > data_time[-1]:
        raise ValueError(
            f"时间窗 [{t_start}, {t_end}] 超出数据时间范围 "
            f"[{data_time[0]}, {data_time[-1]}]。"
        )

    mask = (data_time >= t_start) & (data_time <= t_end)
    num_window = int(np.count_nonzero(mask))

    if num_window < 2:
        raise ValueError("时间窗内的数据点少于2个。")

    window = np.zeros_like(data_time, dtype=np.float64)
    window[mask] = tukey(num_window, alpha=alpha, sym=True)

    return window, mask


def apply_window(data, window):
    output = data.copy()
    output[:, 1:] *= window[:, None]
    return output


def make_output_path(input_path):
    return input_path.with_name(
        f"{input_path.stem}{OUTPUT_SUFFIX}{input_path.suffix}"
    )


def parse_arguments():
    if len(sys.argv) != 4:
        program = Path(sys.argv[0]).name
        print(
            "用法：\n"
            f"  python3 {program} 输入波形文件 起始时间 结束时间\n\n"
            "示例：\n"
            f"  python3 {program} wave_s1_r1.txt 4.0 12.0",
            file=sys.stderr,
        )
        raise SystemExit(2)

    file_name = sys.argv[1]

    try:
        t_start = float(sys.argv[2])
        t_end = float(sys.argv[3])
    except ValueError as error:
        raise ValueError("起始时间和结束时间必须是数值。") from error

    return file_name, t_start, t_end


def main():
    file_name, t_start, t_end = parse_arguments()
    input_path, data = read_waveform(file_name)
    data_time = data[:, 0]

    if data_time[0] < 0.0:
        data_time[0] = 0.0


    # print(f"开始时间: {data_time[0]}")
    # print(f"结束时间: {data_time[-1]}")
    # print(f"时间窗: [{t_start}, {t_end}]")

    # # print(f"开始时间窗创建: {time.time()}")

    # time.sleep(1000)

    window, mask = create_tukey_time_window(
        data_time,
        t_start,
        t_end,
        TUKEY_ALPHA,
    )

    output_data = apply_window(data, window)

    if OUTPUT_WINDOW_COLUMN:
        output_data = np.column_stack((output_data, window))

    output_path = make_output_path(input_path)
    np.savetxt(output_path, output_data, fmt=OUTPUT_FMT)

    dt = float(np.median(np.diff(data_time)))

    # print("窗函数处理完成")
    # print(f"输入文件       : {input_path}")
    # print(f"输出文件       : {output_path}")
    # print(f"数据点数       : {data.shape[0]}")
    # print(f"波形分量数     : {data.shape[1] - 1}")
    # print(f"采样时间间隔   : {dt:.10e}")
    # print(f"时间窗         : [{t_start}, {t_end}]")
    # print(f"窗内采样点数   : {np.count_nonzero(mask)}")
    # print(f"Tukey alpha    : {TUKEY_ALPHA}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"Error: {error}", file=sys.stderr)
        raise SystemExit(1)
