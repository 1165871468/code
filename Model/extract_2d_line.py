#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import time
from pathlib import Path
from read_file import read_inversion_config3d

######################################################
############## read configuration file ###############
######################################################
config_file = Path("../../Input/parameters_inversion3d.config")
config = read_inversion_config3d(config_file)

######################################################
############## read source and receiver file #########
######################################################
SOURCE_FILE = Path("../../Input/source_3d.txt")
RECEIVER_FILE = Path("../../Input/receiver_3d.txt")
OUTPUT_FILE = Path("../../Output/Model/profile2d/line.txt")

X_MIN, X_MAX = config.x1, config.x2
Y_MIN, Y_MAX = config.y1, config.y2

X_COL, Y_COL, Z_COL = 0, 1, 2
EPS = 1.0e-8
######################################################

def read_xyz(filename):
    data = np.loadtxt(filename, dtype=np.float64, ndmin=2)
    return data[:, [X_COL, Y_COL, Z_COL]]


def project_to_z0(point):
    return np.array([point[0], point[1], 0.0], dtype=np.float64)


def unique_points(points):
    out = []
    for p in points:
        if not any(np.allclose(p, q, rtol=0.0, atol=EPS) for q in out):
            out.append(p)
    return out


def line_boundary_intersections(p1, p2):
    x1, y1 = p1[:2]
    x2, y2 = p2[:2]
    dx, dy = x2 - x1, y2 - y1

    points = []

    if abs(dx) > EPS:
        for xb in (X_MIN, X_MAX):
            t = (xb - x1) / dx
            y = y1 + t * dy
            if Y_MIN - EPS <= y <= Y_MAX + EPS:
                points.append(np.array([xb, np.clip(y, Y_MIN, Y_MAX), 0.0]))

    if abs(dy) > EPS:
        for yb in (Y_MIN, Y_MAX):
            t = (yb - y1) / dy
            x = x1 + t * dx
            if X_MIN - EPS <= x <= X_MAX + EPS:
                points.append(np.array([np.clip(x, X_MIN, X_MAX), yb, 0.0]))

    points = unique_points(points)

    if len(points) != 2:
        raise ValueError(f"boundary intersections = {len(points)}, expected 2")

    return points[0], points[1]


def order_endpoints(p1, p2):
    # 优先从左到右；若 x 相同，则从上到下（y 大 -> y 小）
    if p1[0] < p2[0] - EPS:
        return p1, p2
    if p2[0] < p1[0] - EPS:
        return p2, p1
    return (p1, p2) if p1[1] >= p2[1] else (p2, p1)


def build_line(source, receiver):
    src0 = project_to_z0(source)
    rec0 = project_to_z0(receiver)

    # 投影位置相同：采用从左到右的水平测线
    if np.allclose(src0[:2], rec0[:2], rtol=0.0, atol=EPS):
        y = src0[1]
        return (
            np.array([X_MIN, y, 0.0]),
            np.array([X_MAX, y, 0.0]),
        )

    p1, p2 = line_boundary_intersections(src0, rec0)
    return order_endpoints(p1, p2)

def distance_3d(p1, p2):
    p1 = np.asarray(p1, dtype=np.float64)
    p2 = np.asarray(p2, dtype=np.float64)

    return np.linalg.norm(p2 - p1)

def main():
    sources = read_xyz(SOURCE_FILE)
    receivers = read_xyz(RECEIVER_FILE)

    num_sources = sources.shape[0]
    num_receivers = receivers.shape[0]

    if num_sources != config.num_sources:
        raise ValueError(f"number of sources = {num_sources}, expected {config.num_sources}")
    if num_receivers != config.num_receivers:
        raise ValueError(f"number of receivers = {num_receivers}, expected {config.num_receivers}")

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

    with OUTPUT_FILE.open("w", encoding="utf-8") as f:
        for ns, source in enumerate(sources, start=1):
            for nr, receiver in enumerate(receivers, start=1):
                start, end = build_line(source, receiver)
                distance = distance_3d(start, end)

                f.write(
                    f"{ns:d} {nr:d} "
                    f"{start[0]:.6f} {start[1]:.6f} {start[2]:.6f} "
                    f"{end[0]:.6f} {end[1]:.6f} {end[2]:.6f} {distance:.6f}\n"
                )
    print('\n--------------------------------')
    print(f"line for source and receiver: {OUTPUT_FILE}")
    print('--------------------------------\n')

if __name__ == "__main__":
    main()
