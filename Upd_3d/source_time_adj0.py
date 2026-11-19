#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import correlate, correlation_lags
from scipy.signal.windows import tukey

EPS = 1.0e-100
CROSS_CORRELATION_WINDOW_LEN = 4.0

def read_waveform(filename):
    """读取三列波形文件：time, X, Z。"""
    filename = Path(filename)

    if not filename.is_file():
        raise FileNotFoundError(f"cannot open waveform file: {filename}")

    data = np.loadtxt(filename, ndmin=2)

    if data.shape[0] < 5:
        raise ValueError(f"{filename} must contain at least 5 rows")

    if data.shape[1] < 3:
        raise ValueError(f"{filename} must contain at least 3 columns")

    data = data[:, :3]

    if np.any(np.diff(data[:, 0]) <= 0.0):
        raise ValueError(f"time column must be strictly increasing: {filename}")

    return data


def normalize_amplitude(data):
    """按最大绝对振幅归一化。"""
    data = np.asarray(data, dtype=float)
    max_abs = np.max(np.abs(data))

    if max_abs < EPS:
        return data.copy()

    return data / max_abs


def resample_to_common_time(t, amplitude, common_time):
    """使用 scipy.interpolate.interp1d 进行线性插值。"""
    interpolator = interp1d(
        t,
        amplitude,
        kind="linear",
        bounds_error=True,
        assume_sorted=True,
    )
    return interpolator(common_time)

def cross_correlation(file_obs, file_syn, t_start, t_end, component):
    obs = read_waveform(file_obs)
    syn = read_waveform(file_syn)

    t_obs = obs[:, 0]
    t_syn = syn[:, 0]

    a_obs = obs[:, component]
    a_syn = syn[:, component]

    if not np.allclose(t_obs, t_syn, rtol=1e-10, atol=1e-12):
        raise ValueError("Obs and Syn time axes are inconsistent")

    mask = (t_obs >= t_start) & (t_obs <= t_end)

    y_obs = a_obs[mask]
    y_syn = a_syn[mask]

    y_obs = normalize_amplitude(y_obs)
    y_syn = normalize_amplitude(y_syn)

    if y_obs.size < 2:
        raise ValueError("too few samples in cross-correlation window")

    delta = t_obs[1] - t_obs[0]

    corr = correlate(y_syn, y_obs, mode="full", method="auto")
    lags = correlation_lags(y_syn.size, y_obs.size, mode="full")

    best_index = np.argmax(corr)
    best_shift = lags[best_index]

    norm = np.linalg.norm(y_obs) * np.linalg.norm(y_syn)

    if norm < EPS:
        coeff = 0.0
    else:
        coeff = corr[best_index] / norm

    time_shift = -best_shift * delta

    return float(time_shift), float(coeff)

def first_order_derivative(data):
    """
    使用 numpy.gradient 计算一阶导数。

    edge_order=2：
        边界使用二阶精度；
        内部自动采用中心差分。
    """
    data = np.asarray(data, dtype=float)

    if data.shape[0] < 3:
        raise ValueError("at least 3 samples are required")

    time = data[:, 0]

    result = np.empty_like(data)
    result[:, 0] = time
    result[:, 1] = np.gradient(data[:, 1], time, edge_order=2)
    result[:, 2] = np.gradient(data[:, 2], time, edge_order=2)

    return result


def second_order_derivative(data):
    """使用两次 numpy.gradient 计算二阶导数。"""
    first = first_order_derivative(data)
    time = data[:, 0]

    result = np.empty_like(data)
    result[:, 0] = time
    result[:, 1] = np.gradient(first[:, 1], time, edge_order=2)
    result[:, 2] = np.gradient(first[:, 2], time, edge_order=2)

    return result


def source_adj(data_ori, time_res_x, time_res_z):
    """计算 X、Z 两个分量的伴随源时间函数。"""
    data_ori = np.asarray(data_ori, dtype=float)

    der1 = first_order_derivative(data_ori)
    der2 = second_order_derivative(data_ori)

    time = data_ori[:, 0]

    # 用 numpy.trapezoid 代替手工 sum*dt。
    # 若 NumPy 版本较老没有 trapezoid，则退回 trapz。
    integrate = getattr(np, "trapezoid", np.trapz)

    nrx_den = integrate(data_ori[:, 1] * der2[:, 1], x=time)
    nrz_den = integrate(data_ori[:, 2] * der2[:, 2], x=time)

    if abs(nrx_den) < EPS:
        raise ZeroDivisionError("X adjoint-source normalization denominator is zero")

    if abs(nrz_den) < EPS:
        raise ZeroDivisionError("Z adjoint-source normalization denominator is zero")

    nrx = time_res_x / nrx_den
    nrz = time_res_z / nrz_den

    result = np.column_stack((time, -der1[::-1, 1] * nrx, -der1[::-1, 2] * nrz))

    return result


def get_source_time_adj(ns, nr, dt, arrival_time):
    if dt <= EPS:
        raise ValueError("dt must be > 0")

    # file_obs = Path(f"../../Obs/Receivers_obs/wave_s{ns}_r{nr}.txt")
    # file_syn = Path(f"../../Output/Receivers_syn/wave_s{ns}_r{nr}.txt")
    file_obs = Path(f"../../Obs/Receivers_obs/wave_s{ns}_r{nr}_win.txt")
    file_syn = Path(f"../../Output/Receivers_syn/wave_s{ns}_r{nr}_win.txt")

    file_source = Path(f"../../Output/Source_adj/Source_times{ns}_r{nr}.txt")

    win_start = max(0.0, arrival_time - CROSS_CORRELATION_WINDOW_LEN)
    win_end = arrival_time

    # print(f"Ns             = {ns}")
    # print(f"Nr             = {nr}")
    # print(f"dt             = {dt}")
    # print(f"T              = {arrival_time}")
    # print(f"window         = {win_start:.6f} ~ {win_end:.6f}")
    # print(f"observed       = {file_obs}")
    # print(f"synthetic      = {file_syn}")

    data_ori = read_waveform(file_syn)

    file_dt = np.median(np.diff(data_ori[:, 0]))
    if not np.isclose(file_dt, dt, rtol=1.0e-6, atol=EPS):
        print(
            f"Warning: input dt={dt:.12g}, "
            f"waveform dt={file_dt:.12g}"
        )

    time_x, ratio_x = cross_correlation(
        file_obs, file_syn, win_start, win_end, component=1
    )
    time_z, ratio_z = cross_correlation(
        file_obs, file_syn, win_start, win_end, component=2
    )

    # print(f"timeX          = {time_x:.6e}")
    # print(f"ratioX         = {ratio_x:.6e}")
    # print(f"timeZ          = {time_z:.6e}")
    # print(f"ratioZ         = {ratio_z:.6e}")

    adj = source_adj(data_ori, time_x, time_z)

    file_source.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(file_source, adj, fmt="%.4f %.4e %.4e")

    # print(f"output         = {file_source}")

    return adj


def main():
    if len(sys.argv) != 5:
        print("Usage:")
        print("  python3 SourceTime.py Ns Nr dt T")
        print("")
        print("Example:")
        print("  python3 SourceTime.py 1 5 0.01 12.35")
        sys.exit(1)

    try:
        ns = int(sys.argv[1])
        nr = int(sys.argv[2])
        dt = float(sys.argv[3])
        arrival_time = float(sys.argv[4])

        get_source_time_adj(ns, nr, dt, arrival_time)

    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
