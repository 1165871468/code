"""按 line.txt 批量生成二维四节点网格, 使用配套 read_file.py 
运行 python build_topo_mesh.py; 依赖 NumPy; 路径相对于运行时工作目录 
line.txt 九列: source_id receiver_id xs ys zs xe ye ze length 
长度以第九列为准, 单位为米, 模型 x 从 0 开始, z 向下为正 
物理区长度向上取整为 dis 的整数倍; PML 另加, 不包含在该长度内 
深度保留配置值, nz=ceil(depth/dis), 基准 dz=depth/nz 
正弦地形使每列物理区层厚变化, 底部边界保持水平 
模型: Mesh_s1_r1.txt 等; 汇总: mesh_summary.txt 
"""

import math
import time
import numpy as np
from pathlib import Path
from read_file import read_forward_config

CONFIG_FILE = Path('../../Input/parameters_forward.config')
LINE_FILE = Path('../../Output/Model/profile2d/line.txt')
OUTPUT_DIR = Path('../../Output/Mesh')
AMPLITUDE = 0.0
HALF_WAVES = 6
WRITE_INTERMEDIATE = False

def build_mesh(nx, nz, pmlx, pmlz, xmin, xmax, ztop, zbottom, amplitude, half_waves):
    """返回节点坐标、1 起始的连接表、材料号及 x 优先的中间节点表 """
    for name, value, minimum in [('nx', nx, 1), ('nz', nz, 1), ('pmlx', pmlx, 0), ('pmlz', pmlz, 0)]:
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < minimum:
            raise ValueError(f'{name} 必须是 >= {minimum} 的整数')
    values = [xmin, xmax, ztop, zbottom, amplitude, half_waves]
    if not np.isfinite(values).all() or xmax <= xmin or zbottom <= ztop:
        raise ValueError('坐标及地形参数必须有限, 且 xmax > xmin、zbottom > ztop')
    if amplitude < 0 or half_waves < 0 or amplitude >= zbottom - ztop:
        raise ValueError('振幅和半波数须非负, 振幅须小于物理区厚度')

    dx = (xmax - xmin) / nx
    dz = (zbottom - ztop) / nz
    x = xmin + np.arange(-pmlx, nx + pmlx + 1) * dx
    fraction = np.clip((x - xmin) / (xmax - xmin), 0.0, 1.0)
    surface = ztop - amplitude * np.sin(half_waves * np.pi * fraction)

    # 每列物理区均分为 nz 层, 顶面贴合地形、底面固定为 zbottom 
    t = np.linspace(0.0, 1.0, nz + 1)[:, None]
    z = (1.0 - t) * surface[None, :] + t * zbottom
    if pmlz:
        bottom = zbottom + np.arange(1, pmlz + 1)[:, None] * dz
        z = np.vstack((z, np.broadcast_to(bottom, (pmlz, x.size))))
    if not np.all(np.diff(x) > 0) or not np.all(np.diff(z, axis=0) > 0):
        raise ValueError('出现退化或倒置网格, 请检查坐标范围和单元数')
    grid = np.stack((np.broadcast_to(x, z.shape), z), axis=-1)
    nodes = grid.reshape(-1, 2)

    # 顺序与旧程序完全一致: [下一层右、下一层左、本层左、本层右] 
    ids = np.arange(1, nodes.shape[0] + 1).reshape(z.shape)
    elements = np.stack((ids[1:, 1:], ids[1:, :-1], ids[:-1, :-1], ids[:-1, 1:]), axis=-1).reshape(-1, 4)
    material = np.full((nz + pmlz, nx + 2 * pmlx), 2, dtype=int)
    material[:nz, pmlx:pmlx + nx] = 1
    return nodes, elements, material.ravel(), grid.transpose(1, 0, 2).reshape(-1, 2)

def write_mesh(path, nodes, elements, material):
    """物理区先写、PML 后写, 保留原始单元 ID; 坐标保留充分精度 """
    with Path(path).open('w', encoding='utf-8') as fp:
        fp.write(f'number of elements & nodes & ndime:\n{len(elements)} {len(nodes)} 2\ncoordinates\n')
        for index, (x, z) in enumerate(nodes, 1):
            fp.write(f'{index} {x:.15g} {z:.15g}\n')
        fp.write('End coordinates\nelements\n')
        for region in (1, 2):
            for index in np.flatnonzero(material == region):
                a, b, c, d = elements[index]
                fp.write(f'{index + 1} {a} {b} {c} {d} {region}\n')
        fp.write('End elements\n')

def read_mesh_parameters(path):
    params = read_forward_config(path)
    missing = [key for key in ('dis', 'Layx', 'Layz', 'model_size') if key not in params]

    if missing:
        raise ValueError(f'{path}: 缺少参数 {missing}; 模型长度、深度行可标注 // model_size')

    # print(params)
    # print(missing)
    # time.sleep(1000)

    spacing = float(params['dis'])
    size = np.asarray(params['model_size'], dtype=float)
    if size.shape != (2,) or not np.isfinite(size).all() or (size <= 0).any():
        raise ValueError('model_size 必须为两个正数: 参考长度、模型深度')
    depth = float(size[1])
    if not np.isfinite(spacing) or spacing <= 0:
        raise ValueError('dis 必须为有限正数')
    layers = []
    for key in ('Layx', 'Layz'):
        value = float(params[key])
        if not np.isfinite(value) or value < 0 or value != int(value):
            raise ValueError(f'{key} 必须为非负整数')
        layers.append(int(value))
    return spacing, depth, *layers

def covering_cells(length, spacing):
    """向上取整; 修正浮点运算使整倍数意外多出一层的情况 """
    count = max(1, math.ceil(length / spacing))
    if count > 1 and (count - 1) * spacing >= length:
        count -= 1
    if count * spacing < length:
        count += 1
    return count

def generate_meshes(config_file, line_file, output_dir):
    spacing, depth, pmlx, pmlz = read_mesh_parameters(config_file)
    lines = np.loadtxt(line_file, ndmin=2)

    # print(spacing, depth, pmlx, pmlz)
    # print(lines)
    # time.sleep(1000)

    ### Check the validity of the line.txt file ###
    if lines.shape[0] == 0 or lines.shape[1] != 9 or not np.isfinite(lines).all():
        raise ValueError('line.txt 必须为非空、有限数值的九列数据')
    ids = lines[:, :2]
    if (ids < 1).any() or (ids != np.floor(ids)).any():
        raise ValueError('震源和台站编号必须是正整数')
    pairs = [(int(s), int(r)) for s, r in ids]
    if len(set(pairs)) != len(pairs):
        raise ValueError('line.txt 存在重复震源—台站编号, 会导致输出文件重名')
    if (lines[:, 8] <= 0).any():
        raise ValueError('测线长度必须大于零')
    if not np.isfinite([AMPLITUDE, HALF_WAVES]).all() or not 0 <= AMPLITUDE < depth or HALF_WAVES < 0:
        raise ValueError('地形参数无效: 要求 0 <= AMPLITUDE < depth, HALF_WAVES >= 0')
    ### End of checking ###
    nz = covering_cells(depth, spacing)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = []
    for row, (source_id, receiver_id) in zip(lines, pairs):
        ### Get the parameters for the mesh ###
        length = float(row[8])
        nx = covering_cells(length, spacing)
        model_length = nx*spacing
        nodes, elements, material, profile = build_mesh(nx, nz, pmlx, pmlz, 0.0, model_length, 0.0, depth, AMPLITUDE, HALF_WAVES)
        ### End of getting parameters ###
        ### Write the mesh file ###
        suffix = f's{source_id}_r{receiver_id}'
        mesh_file = output_dir / f'Mesh_{suffix}.txt'
        write_mesh(mesh_file, nodes, elements, material)
        ### End of writing ###
        if WRITE_INTERMEDIATE:
            x = np.linspace(0.0, model_length, 1001)
            height = AMPLITUDE * (np.sin(HALF_WAVES * np.pi * x / model_length) + 1.0)
            np.savetxt(output_dir / f'topo3_{suffix}.txt', np.column_stack((x, height)), fmt='%.15g')
            np.savetxt(output_dir / f'topo_resple1_{suffix}.txt', profile, fmt='%.15g')
            np.savetxt(output_dir / f'topo_resple_{suffix}.txt',
                       np.column_stack((profile[:, 0], depth - profile[:, 1])), fmt='%.15g')

        summary.append((source_id, receiver_id, length, model_length, depth, nx, nz,
                        spacing, depth / nz, pmlx, pmlz, len(nodes), len(elements)))
        print(f'{mesh_file.name}: L={length:.6f} -> {model_length:.6f}, cells={nx} x {nz}')
    header = ('source_id receiver_id line_length model_length depth nx nz dx dz '
              'pml_x pml_z num_nodes num_elements')
    np.savetxt(output_dir / 'mesh_summary.txt', np.asarray(summary), header=header,
               fmt=['%d', '%d', '%.15g', '%.15g', '%.15g', '%d', '%d',
                    '%.15g', '%.15g', '%d', '%d', '%d', '%d'])
    print(f'Done: {len(summary)} meshes in {output_dir.resolve()}')
    return summary

if __name__ == '__main__':
    print("\n--------------------------------")
    print("Building topography mesh: \n")

    generate_meshes(CONFIG_FILE, LINE_FILE, OUTPUT_DIR)

    print("--------------------------------\n")