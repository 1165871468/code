#!/usr/bin/env python3
"""三维模型初始化, 调用同目录下用户提供的 read_file.py;

依赖 numpy、scipy; 运行 python Initial.py;需要 Python >= 3.10 (read_file.py) ;
所有相对路径均相对于运行时工作目录;
输入: parameters_inversion3d.config; 正演 vel_initial.dat, 前六列 x y z vp vs rho;
输出: Forver/vel.dat (原文件字节复制) 、Forver/per_vel.dat、
      Inver/vel_initial.dat、Inver/vel.dat、Inver/per_vel.dat;
模型输出六列 x y z vp vs rho; 扰动输出六列 x y z delta_vp delta_vs delta_rho;
注意: read_file.read_perturbation 只读取前四列, 不会返回三个扰动分量;
反演节点按 x 外层、y 中层、z 内层排列;单位不变, 最近邻按三维欧氏距离;
等距时选择输入文件中最先出现的节点;查询不限制距离, 域外节点取最近点;
默认保留原程序4位小数; 精细网格可自行修改 OUTPUT_FORMAT;
"""

import shutil
import numpy as np
import time
from pathlib import Path
from scipy.spatial import cKDTree
from read_file import read_inversion_config3d, read_initial_model

# 输入和输出路径
CONFIG_FILE = Path('../../Input/parameters_inversion3d.config')
FORWARD_DIR = Path('../../Output/Model/Forver')
INVERSION_DIR = Path('../../Output/Model/Inver')
BATCH_SIZE = 10000
OUTPUT_FORMAT = '%.4f'  # 列之间显式留空格, 防止大数值粘连

# File_ini3d = "vel_obs3d.txt"
File_ini3d = "vel_ini3d.txt"
File_vel3d = "vel3d.txt"
File_per_vel3d = "per_vel3d.txt"

def grid_batches(config, batch_size):
    """分批生成 xyz, 内存不随完整反演网格节点数增长;"""
    nx, ny, nz = config.nx_inv, config.ny_inv, config.nz_inv
    lower = np.array([config.x1, config.y1, config.z1])
    step = (np.array([config.x2, config.y2, config.z2]) - lower) / (np.array([nx, ny, nz]) - 1)
    for start in range(0, nx * ny * nz, batch_size):
        ids = np.arange(start, min(start + batch_size, nx * ny * nz), dtype=np.int64)
        indices = np.column_stack((ids // (ny * nz), (ids // nz) % ny, ids % nz))
        yield lower + indices * step


def nearest_indices(tree, points, workers):
    """一次查询同时服务 vp、vs、rho; 并列最近点使用确定的输入顺序;"""
    if tree.n == 1:
        return np.zeros(len(points), dtype=np.intp)
    distances, indices = tree.query(points, k=2, workers=workers)
    chosen = indices[:, 0].copy()
    for row in np.flatnonzero(distances[:, 0] == distances[:, 1]):
        candidates = np.asarray(tree.query_ball_point(
            points[row], np.nextafter(distances[row, 0], np.inf)), dtype=np.intp)
        difference = tree.data[candidates] - points[row]
        squared = np.einsum('ij,ij->i', difference, difference)
        chosen[row] = candidates[squared == squared.min()].min()
    return chosen


def initialize(config_file, forward_dir, inversion_dir, batch_size=BATCH_SIZE):
    """初始化文件并返回正演、反演节点数; 会覆盖指定目录下同名输出;"""
    if isinstance(batch_size, bool) or not isinstance(batch_size, (int, np.integer)) or batch_size < 1:
        raise ValueError('batch_size must be a positive integer')

    config = read_inversion_config3d(config_file)
    bounds = np.array([config.x1, config.x2, config.y1, config.y2, config.z1, config.z2])

    if not np.isfinite(bounds).all():
        raise ValueError('model bounds must be finite')

    # read_file.py 已校验范围大小、节点数和 num_cores;
    forward_dir, inversion_dir = Path(forward_dir), Path(inversion_dir)

    if forward_dir.resolve() == inversion_dir.resolve():
        raise ValueError('forward and inversion directories must be different')

    source = forward_dir / File_ini3d
    model = read_initial_model(source)

    if len(model) == 0:
        raise ValueError(f'{source}: model must not be empty')

    coordinates = np.ascontiguousarray(model[:, :3])
    tree = cKDTree(coordinates)
    inversion_dir.mkdir(parents=True, exist_ok=True)

    # 零扰动使用单个批次缓冲区, 避免创建整个模型大小的零数组;
    buffer = np.zeros((batch_size, 6))
    with (forward_dir / File_per_vel3d).open('w', encoding='utf-8') as fp:
        for start in range(0, len(model), batch_size):
            xyz = coordinates[start:start + batch_size]
            block = buffer[:len(xyz)]
            block[:, :3] = xyz
            np.savetxt(fp, block, fmt=OUTPUT_FORMAT)

    with (inversion_dir / File_ini3d).open('w', encoding='utf-8') as model_fp, \
         (inversion_dir / File_per_vel3d).open('w', encoding='utf-8') as perturbation_fp:
        for xyz in grid_batches(config, batch_size):
            block = buffer[:len(xyz)]
            block[:, :3] = xyz
            block[:, 3:] = 0.0
            np.savetxt(perturbation_fp, block, fmt=OUTPUT_FORMAT)
            nearest = nearest_indices(tree, xyz, config.num_cores)
            block[:, 3:] = model[nearest, 3:6]
            np.savetxt(model_fp, block, fmt=OUTPUT_FORMAT)

    shutil.copyfile(source, forward_dir / File_vel3d)
    shutil.copyfile(inversion_dir / File_ini3d, inversion_dir / File_vel3d)
    count = config.nx_inv * config.ny_inv * config.nz_inv
    print('\n--------------------------------')
    print(f'Done: forward model nodes = {len(model)}, inversion nodes = {count}')
    print(f'Inversion grid: {config.nx_inv} x {config.ny_inv} x {config.nz_inv}')
    print(f'Forward output: {forward_dir.resolve()}')
    print(f'Inversion output: {inversion_dir.resolve()}')
    print('--------------------------------\n')
    return len(model), count


if __name__ == '__main__':
    initialize(CONFIG_FILE, FORWARD_DIR, INVERSION_DIR)
