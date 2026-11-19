#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from pathlib import Path

import numpy as np
from scipy.signal import correlate, correlation_lags
from scipy.signal.windows import tukey


EPS = 1.0e-50
CROSS_CORRELATION_WINDOW_LEN = 4.0
ADJOINT_WINDOW_ALPHA = 1.0


def read_waveform(filename):
    filename = Path(filename)

    if not filename.is_file():
        raise FileNotFoundError(f"cannot open waveform file: {filename}")

    data = np.loadtxt(filename, dtype=float, ndmin=2)

    if data.shape[0] < 5:
        raise ValueError(f"{filename} must contain at least 5 rows")

    if data.shape[1] < 3:
        raise ValueError(f"{filename} must contain at least 3 columns")

    data = data[:, :3]

    if not np.all(np.isfinite(data)):
        raise ValueError(f"{filename} contains NaN or Inf")

    if np.any(np.diff(data[:, 0]) <= 0.0):
        raise ValueError(f"time column must be strictly increasing: {filename}")

    return data


def cross_correlation(file_obs, file_syn, t_start, t_end, component):
    if component not in (1, 2):
        raise ValueError("component must be 1 (X) or 2 (Z)")

    if t_end <= t_start:
        raise ValueError("invalid cross-correlation window")

    obs = read_waveform(file_obs)
    syn = read_waveform(file_syn)

    t_obs = obs[:, 0]
    t_syn = syn[:, 0]

    if t_obs.shape != t_syn.shape or not np.allclose(
        t_obs, t_syn, rtol=1.0e-10, atol=1.0e-12
    ):
        print(file_obs)
        print(file_syn)
        raise ValueError("Obs and Syn time axes are inconsistent")

    if t_start < t_obs[0] or t_end > t_obs[-1]:
        raise ValueError(
            f"window [{t_start:.6f}, {t_end:.6f}] is outside waveform range "
            f"[{t_obs[0]:.6f}, {t_obs[-1]:.6f}]"
        )

    dt_array = np.diff(t_obs)
    delta = float(np.median(dt_array))

    if not np.allclose(dt_array, delta, rtol=1.0e-10, atol=1.0e-12):
        raise ValueError("waveform sampling is not uniform")

    mask = (t_obs >= t_start) & (t_obs <= t_end)

    y_obs = obs[mask, component]
    y_syn = syn[mask, component]

    if y_obs.size < 2:
        raise ValueError("too few samples in cross-correlation window")

    corr = correlate(y_syn, y_obs, mode="full", method="auto")
    lags = correlation_lags(y_syn.size, y_obs.size, mode="full")

    best_index = int(np.argmax(corr))
    best_shift = int(lags[best_index])
    max_corr = float(corr[best_index])

    norm = float(np.linalg.norm(y_obs) * np.linalg.norm(y_syn))
    coeff = 0.0 if norm < EPS else max_corr / norm

    time_shift = -best_shift * delta

    return float(time_shift), float(coeff)


def first_order_derivative(data):
    data = np.asarray(data, dtype=float)

    if data.ndim != 2 or data.shape[1] < 3:
        raise ValueError("data must contain at least 3 columns")

    if data.shape[0] < 3:
        raise ValueError("at least 3 samples are required")

    time = data[:, 0]

    result = np.empty_like(data)
    result[:, 0] = time
    result[:, 1] = np.gradient(data[:, 1], time, edge_order=2)
    result[:, 2] = np.gradient(data[:, 2], time, edge_order=2)

    return result


def create_phase_window(time, t_start, t_end, alpha=1.0):
    time = np.asarray(time, dtype=float)

    if t_end <= t_start:
        raise ValueError("window end must be greater than window start")

    if not 0.0 <= alpha <= 1.0:
        raise ValueError("window alpha must be within [0, 1]")

    if t_start < time[0] or t_end > time[-1]:
        raise ValueError(
            f"adjoint window [{t_start:.6f}, {t_end:.6f}] is outside "
            f"waveform range [{time[0]:.6f}, {time[-1]:.6f}]"
        )

    mask = (time >= t_start) & (time <= t_end)
    n = int(np.count_nonzero(mask))

    if n < 2:
        raise ValueError("too few samples inside adjoint window")

    window = np.zeros_like(time)
    window[mask] = tukey(n, alpha=alpha, sym=True)

    return window


def source_adj(data_syn, window, time_res_x, time_res_z):
    data_syn = np.asarray(data_syn, dtype=float)
    time = data_syn[:, 0]

    # 先对未加窗的原始合成波形求导
    der1 = first_order_derivative(data_syn)

    du_x = der1[:, 1]
    du_z = der1[:, 2]

    integrate = getattr(np, "trapezoid", np.trapz)

    # N = integral[ w(t) * (du/dt)^2 dt ]
    nx = float(integrate(window * du_x**2, x=time))
    nz = float(integrate(window * du_z**2, x=time))

    if nx < EPS:
        raise ZeroDivisionError("X adjoint-source normalization denominator is zero")

    if nz < EPS:
        raise ZeroDivisionError("Z adjoint-source normalization denominator is zero")

    # f†(t) = (DeltaT/N) * w(t) * du/dt
    # adj_x = (time_res_x / nx) * window * du_x
    # adj_z = (time_res_z / nz) * window * du_z
    adj_x = (1.0 / nx) * window * du_x
    adj_z = (1.0 / nz) * window * du_z

    # 最后做时间反转
    result = np.column_stack((time, adj_x[::-1], adj_z[::-1]))

    return result


def get_source_time_adj(ns, nr, dt, arrival_time):
    if dt <= 0.0:
        raise ValueError("dt must be > 0")

    # 互相关与求导均使用未预加窗的原始波形
    file_obs = Path(f"../../Obs/Receivers_obs/wave_s{ns}_r{nr}.txt")
    file_syn = Path(f"../../Output/Receivers_syn/wave_s{ns}_r{nr}.txt")

    file_source = Path(
        f"../../Output/Source_adj/Source_times{ns}_r{nr}.txt"
    )

    win_start = max(0.0, arrival_time - CROSS_CORRELATION_WINDOW_LEN)
    win_end = arrival_time

    data_syn = read_waveform(file_syn)

    file_dt = float(np.median(np.diff(data_syn[:, 0])))

    if not np.isclose(file_dt, dt, rtol=1.0e-6, atol=1.0e-12):
        print(
            f"Warning: input dt={dt:.12g}, "
            f"waveform dt={file_dt:.12g}"
        )

    # 使用未加窗波形测量走时差
    time_x, ratio_x = cross_correlation(
        file_obs, file_syn, win_start, win_end, component=1
    )

    time_z, ratio_z = cross_correlation(
        file_obs, file_syn, win_start, win_end, component=2
    )

    # 单独构造伴随源震相窗
    window = create_phase_window(
        data_syn[:, 0],
        win_start,
        win_end,
        ADJOINT_WINDOW_ALPHA,
    )

    adj = source_adj(
        data_syn,
        window,
        time_x,
        time_z,
    )

    file_source.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(file_source, adj, fmt="%.4f %.4e %.4e")

    return adj


def main():
    if len(sys.argv) != 5:
        print("Usage:")
        print("  python3 source_time_adj.py Ns Nr dt T")
        print("")
        print("Example:")
        print("  python3 source_time_adj.py 1 5 0.01 12.35")
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
