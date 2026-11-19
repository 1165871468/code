#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import numpy as np
from scipy.signal import correlate, correlation_lags

# ===================== 参数 =====================
CONFIG_FILE = Path("../../Input/parameters_inversion.config")

ARRIVAL_FILE = Path("../../Input/arrival_p.txt")
OBS_ROOT = Path("../../Obs/Receivers_obs")
SYN_ROOT = Path("../../Output/Receivers_syn")

ERROR_FILE = Path("../../Output/Error/Error.txt")
SHIFT_FILE = Path("../../Output/Error/Error_1.txt")

# 共震源双差
RECEIVER_FILE = Path("../../Input/Receiver.txt")
ERROR_DD_FILE = Path("../../Output/Error/Error_dd_cs.txt")

# 与 C 程序一致：仅保留台站间距 <= 此阈值的台站对。
# 单位必须与 Receiver.txt 的坐标单位一致。
MIN_RECEIVER_DISTANCE = 20000.0
# ===============================================

def read_inversion_parameters(filename):
    """
    从 parameters_inversion.config 读取：
        第12个有效参数行 -> CROSS_CORRELATION_WINDOW_LEN
        第13个有效参数行 -> CC_THRESHOLD

    对应当前配置文件：
        4.0  -> 时间窗长度
        0.80 -> 最小波形互相关系数
    """
    filename = Path(filename)

    if not filename.is_file():
        raise FileNotFoundError(f"cannot open config file: {filename}")

    rows = []

    with filename.open("r", encoding="utf-8") as file:
        for line_number, line in enumerate(file, start=1):
            value_text = line.split("c:", 1)[0].strip()

            if not value_text:
                continue

            try:
                values = [float(x) for x in value_text.split()]
            except ValueError as exc:
                raise ValueError(
                    f"invalid numeric value in config line {line_number}: "
                    f"{line.rstrip()}"
                ) from exc

            rows.append((line_number, values))

    expected_counts = [2, 2, 2] + [1] * 12

    if len(rows) != len(expected_counts):
        raise ValueError(
            f"expected {len(expected_counts)} parameter lines, "
            f"but found {len(rows)}"
        )

    for (line_number, values), expected in zip(rows, expected_counts):
        if len(values) != expected:
            raise ValueError(
                f"config line {line_number} requires {expected} value(s), "
                f"but found {len(values)}"
            )

    cross_correlation_window_len = rows[11][1][0]
    cc_threshold = rows[12][1][0]

    if cross_correlation_window_len <= 0.0:
        raise ValueError(
            "cross-correlation window length must be greater than 0"
        )

    if not 0.0 <= cc_threshold <= 1.0:
        raise ValueError(
            "CC threshold must be within [0, 1]"
        )

    return cross_correlation_window_len, cc_threshold

def read_waveform(filename, component):
    """读取三列波形文件：第1列时间,第2列X,第3列Z。"""
    filename = Path(filename)

    if not filename.is_file():
        raise FileNotFoundError(f"cannot open waveform file: {filename}")

    data = np.loadtxt(filename, dtype=float, ndmin=2)

    if data.shape[0] < 3:
        raise ValueError(f"insufficient data in {filename}: {data.shape[0]} rows")

    if data.shape[1] < 3:
        raise ValueError(f"{filename} must contain at least 3 columns")

    if component not in (1, 2):
        raise ValueError("component must be 1 (X) or 2 (Z)")

    time_axis = data[:, 0]
    amplitude = data[:, component]

    if np.any(np.diff(time_axis) <= 0.0):
        raise ValueError(f"time column is not strictly increasing: {filename}")

    if not np.all(np.isfinite(time_axis)) or not np.all(np.isfinite(amplitude)):
        raise ValueError(f"{filename} contains NaN or Inf")

    return time_axis, amplitude


def cross_correlation(file_obs, file_syn, t_start, t_end, component):
    """
    计算观测波形和合成波形的互相关时间差。

    前提：
        1. Obs 与 Syn 使用完全相同的时间轴和采样率；
        2. 输入波形已经完成所需的 taper / window 处理；
        3. 此处不再重复 taper、不做振幅归一化、不做插值重采样。
    """
    if component not in (1, 2):
        raise ValueError("component must be 1 (X) or 2 (Z)")

    if t_end <= t_start:
        raise ValueError("invalid cross-correlation window")

    t_obs, y_obs_all = read_waveform(file_obs, component)
    t_syn, y_syn_all = read_waveform(file_syn, component)

    if t_obs.shape != t_syn.shape or not np.allclose(
        t_obs, t_syn, rtol=1.0e-10, atol=1.0e-12
    ):
        raise ValueError("Obs and Syn time axes are inconsistent")

    if t_start < t_obs[0] or t_end > t_obs[-1]:
        raise ValueError(
            f"window [{t_start:.6f}, {t_end:.6f}] is outside waveform range "
            f"[{t_obs[0]:.6f}, {t_obs[-1]:.6f}]"
        )

    dt_array = np.diff(t_obs)
    delta = float(np.median(dt_array))

    if delta <= 0.0:
        raise ValueError("invalid sampling interval")

    if not np.allclose(dt_array, delta, rtol=1.0e-10, atol=1.0e-12):
        raise ValueError("waveform sampling is not uniform")

    mask = (t_obs >= t_start) & (t_obs <= t_end)
    y_obs = y_obs_all[mask]
    y_syn = y_syn_all[mask]

    if y_obs.size < 2:
        raise ValueError("too few samples inside cross-correlation window")

    corr = correlate(y_syn, y_obs, mode="full", method="auto")
    lags = correlation_lags(y_syn.size, y_obs.size, mode="full")

    best_index = int(np.argmax(corr))
    best_shift = int(lags[best_index])
    max_corr = float(corr[best_index])

    norm = float(np.linalg.norm(y_obs) * np.linalg.norm(y_syn))
    coeff = 0.0 if norm == 0.0 else max_corr / norm

    time_shift = -best_shift * delta

    return float(time_shift), float(coeff)


def load_arrival_file(filename):
    filename = Path(filename)

    if not filename.is_file():
        raise FileNotFoundError(f"cannot open arrival file: {filename}")

    data = np.loadtxt(filename)

    if data.ndim == 1:
        data = data.reshape(1, -1)

    if data.shape[1] < 3:
        raise ValueError(f"{filename} must contain at least 3 columns")

    return data[:, :3]



def load_receiver_file(filename):
    """读取 Receiver.txt,前两列为台站 x、z/y 坐标；receiver_id 从1开始。"""
    filename = Path(filename)

    if not filename.is_file():
        raise FileNotFoundError(f"cannot open receiver file: {filename}")

    data = np.loadtxt(filename, dtype=float, ndmin=2)

    if data.shape[1] < 2:
        raise ValueError(f"{filename} must contain at least 2 columns")

    if not np.all(np.isfinite(data[:, :2])):
        raise ValueError(f"{filename} contains NaN or Inf")

    return data[:, :2]


def receiver_distance(receiver_data, id1, id2):
    """计算两个台站的二维距离。receiver_id 从1开始。"""
    p1 = receiver_data[id1 - 1]
    p2 = receiver_data[id2 - 1]
    return float(np.hypot(p1[0] - p2[0], p1[1] - p2[1]))


def calculate_common_source_double_difference(
    error_file,
    receiver_file,
    output_file,
    max_receiver_distance,
):
    """
    共震源双差：
        固定 source_id, 对该震源下所有有效 receiver 两两组合。

    单差记录:
        source_id receiver_id error cc

    双差:
        dd = error(source, receiver1) - error(source, receiver2)

    仅保存:
        distance(receiver1, receiver2) <= max_receiver_distance

    输出:
        source_id receiver_id_1 receiver_id_2 error_diff
    """
    error_file = Path(error_file)
    output_file = Path(output_file)

    if not error_file.is_file():
        raise FileNotFoundError(f"cannot open error file: {error_file}")

    error_data = np.loadtxt(error_file, dtype=float, ndmin=2)

    if error_data.size == 0:
        output_file.parent.mkdir(parents=True, exist_ok=True)
        output_file.write_text("", encoding="utf-8")
        print("Saved double differences: 0")
        return 0

    if error_data.shape[1] < 4:
        raise ValueError(f"{error_file} must contain at least 4 columns")

    receiver_data = load_receiver_file(receiver_file)

    max_receiver_id = int(np.max(error_data[:, 1]))
    if receiver_data.shape[0] < max_receiver_id:
        raise ValueError(
            f"Receiver file has {receiver_data.shape[0]} rows, "
            f"but Error.txt uses receiver_id up to {max_receiver_id}"
        )

    # {(source_id, receiver_id): error}
    error_map = {}
    source_receivers = {}

    for row in error_data:
        ns = int(round(row[0]))
        nr = int(round(row[1]))
        error = float(row[2])

        error_map[(ns, nr)] = error
        source_receivers.setdefault(ns, []).append(nr)

    output_file.parent.mkdir(parents=True, exist_ok=True)

    num_saved = 0
    num_skipped_distance = 0

    with output_file.open("w", encoding="utf-8") as fout:
        for ns in sorted(source_receivers):
            # 去重并排序,保证 id1 < id2,与 C 程序一致。
            receivers = sorted(set(source_receivers[ns]))

            for i in range(len(receivers) - 1):
                id1 = receivers[i]

                for j in range(i + 1, len(receivers)):
                    id2 = receivers[j]

                    distance = receiver_distance(receiver_data, id1, id2)

                    # C 程序实际逻辑是距离大于阈值则跳过。
                    if distance > max_receiver_distance:
                        num_skipped_distance += 1
                        continue

                    error1 = error_map[(ns, id1)]
                    error2 = error_map[(ns, id2)]
                    error_diff = error1 - error2

                    fout.write(
                        f"{ns:d} {id1:d} {id2:d} {error_diff:.16e}\n"
                    )
                    num_saved += 1

    print(f"Maximum receiver distance : {max_receiver_distance:.6f}")
    print(f"Saved double differences : {num_saved}")
    print(f"Skipped by distance      : {num_skipped_distance}")

    return num_saved

def main():
    cross_correlation_window_len, cc_threshold = read_inversion_parameters(
        CONFIG_FILE
    )

    arrival_data = load_arrival_file(ARRIVAL_FILE)

    ERROR_FILE.parent.mkdir(parents=True, exist_ok=True)
    SHIFT_FILE.parent.mkdir(parents=True, exist_ok=True)

    accepted_count = 0
    processed_count = 0

    print(f"cross-correlation window = {cross_correlation_window_len:.6f} s")
    print(f"CC threshold             = {cc_threshold:.6f}")

    # time.sleep(1000)

    with ERROR_FILE.open("w", encoding="utf-8") as fout_error, \
         SHIFT_FILE.open("w", encoding="utf-8") as fout_shift:

        for row in arrival_data:
            ns = int(row[0])
            nr = int(row[1])
            arrival_time = float(row[2])

            file_obs = OBS_ROOT / f"wave_s{ns}_r{nr}.txt"
            file_syn = SYN_ROOT / f"wave_s{ns}_r{nr}.txt"

            # file_obs = OBS_ROOT / f"wave_s{ns}_r{nr}_win.txt"
            # file_syn = SYN_ROOT / f"wave_s{ns}_r{nr}_win.txt"

            win_start = max(0.0, arrival_time - cross_correlation_window_len)
            win_end = arrival_time

            try:
                time_x, ratio_x = cross_correlation(file_obs, file_syn, win_start, win_end, component=1)
                time_z, ratio_z = cross_correlation(file_obs, file_syn, win_start, win_end, component=2)

            except (FileNotFoundError, ValueError) as exc:
                # 批处理中单个文件有问题时跳过,而不是终止全部计算
                print(f"Skip Ns={ns} Nr={nr}: {exc}")
                continue

            error = 0.5 * time_z**2
            error = time_z

            if ratio_z >= cc_threshold:
                fout_error.write(
                    f"{ns:d} {nr:5d} {error:20.12e} {ratio_z:15.4e}\n"
                )
                accepted_count += 1

            fout_shift.write(
                f"{ns:d} {nr:5d} {time_z:20.12e} {ratio_z:15.4e}\n"
            )

            processed_count += 1

    print(
        f"total={arrival_data.shape[0]} "
        f"processed={processed_count} "
        f"accepted={accepted_count}"
    )

    calculate_common_source_double_difference(
        ERROR_FILE,
        RECEIVER_FILE,
        ERROR_DD_FILE,
        MIN_RECEIVER_DISTANCE,
    )


if __name__ == "__main__":
    main()
