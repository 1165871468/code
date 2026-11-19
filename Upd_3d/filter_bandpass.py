#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from pathlib import Path

import numpy as np
from scipy.signal import butter, sosfiltfilt


LOW_CUTOFF_FREQUENCY_HZ = 0.1
HIGH_CUTOFF_FREQUENCY_HZ = 4.0
FILTER_ORDER = 2

TIME_STEP_RTOL = 1.0e-8
OUTPUT_FMT = "%.15e"


def read_waveform(filename):
    filename = Path(filename)

    if not filename.is_file():
        raise FileNotFoundError(f"cannot open file: {filename}")

    data = np.loadtxt(filename, dtype=np.float64, ndmin=2)

    if data.shape[1] < 3:
        raise ValueError(f"{filename} must contain at least 3 columns")

    if data.shape[0] < 10:
        raise ValueError(f"{filename} contains too few samples")

    if not np.all(np.isfinite(data[:, :3])):
        raise ValueError(f"{filename} contains NaN or Inf")

    return data[:, 0], data[:, 1], data[:, 2]


def get_sampling_rate(time):
    dt_array = np.diff(time)

    if np.any(dt_array <= 0.0):
        raise ValueError("time column must be strictly increasing")

    dt = float(np.median(dt_array))

    if not np.allclose(
        dt_array,
        dt,
        rtol=TIME_STEP_RTOL,
        atol=TIME_STEP_RTOL * abs(dt),
    ):
        raise ValueError("time sampling is not uniform")

    return dt, 1.0 / dt


def bandpass_zero_phase(signal, fs, fmin, fmax, order=2):
    signal = np.asarray(signal, dtype=np.float64)

    if np.max(np.abs(signal)) == 0.0:
        return np.zeros_like(signal)

    nyquist = 0.5 * fs

    if not (0.0 < fmin < fmax < nyquist):
        raise ValueError(
            f"invalid bandpass range: {fmin} - {fmax} Hz; "
            f"Nyquist = {nyquist} Hz"
        )

    sos = butter(
        N=order,
        Wn=[fmin, fmax],
        btype="bandpass",
        fs=fs,
        output="sos",
    )

    return sosfiltfilt(sos, signal)


def parse_command_line():
    if len(sys.argv) not in (2, 3):
        name = Path(sys.argv[0]).name
        print(f"Usage: python3 {name} input_file [output_file]")
        sys.exit(1)

    input_file = Path(sys.argv[1])

    if len(sys.argv) == 3:
        output_file = Path(sys.argv[2])
    else:
        suffix = input_file.suffix if input_file.suffix else ".txt"
        output_file = input_file.with_name(
            f"{input_file.stem}_bandpass{suffix}"
        )

    return input_file, output_file


def main():
    input_file, output_file = parse_command_line()

    time, x, z = read_waveform(input_file)
    dt, fs = get_sampling_rate(time)

    x_bp = bandpass_zero_phase(
        x, fs,
        LOW_CUTOFF_FREQUENCY_HZ,
        HIGH_CUTOFF_FREQUENCY_HZ,
        FILTER_ORDER,
    )

    z_bp = bandpass_zero_phase(
        z, fs,
        LOW_CUTOFF_FREQUENCY_HZ,
        HIGH_CUTOFF_FREQUENCY_HZ,
        FILTER_ORDER,
    )

    output_file.parent.mkdir(parents=True, exist_ok=True)

    np.savetxt(
        output_file,
        np.column_stack((time, x_bp, z_bp)),
        fmt=OUTPUT_FMT,
    )

    # print(f"dt       = {dt:.12g} s")
    # print(f"fs       = {fs:.12g} Hz")
    # print(
    #     f"bandpass = {LOW_CUTOFF_FREQUENCY_HZ} - "
    #     f"{HIGH_CUTOFF_FREQUENCY_HZ} Hz"
    # )
    # print(f"order    = {FILTER_ORDER}")
    # print(f"output   = {output_file}")


if __name__ == "__main__":
    main()
