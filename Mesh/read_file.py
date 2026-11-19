#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from dataclasses import dataclass
from pathlib import Path
import numpy as np
from dataclasses import dataclass
from pathlib import Path


CONFIG_FILE = Path("../../Input/parameters_inversion3d.config")


@dataclass
class ModelDomain3D:
    """三维模型空间范围。"""
    x1: float
    x2: float
    y1: float
    y2: float
    z1: float
    z2: float


@dataclass
class GridSize3D:
    """三维网格节点数。"""
    nx: int
    ny: int
    nz: int


@dataclass
class DataInfo:
    """震源和台站信息。"""
    num_sources: int
    num_receivers: int


@dataclass
class KernelConfig:
    """敏感核计算相关参数。"""
    integration_interval: int
    wavelength: float


@dataclass
class MeasurementConfig:
    """波形互相关和数据筛选参数。"""
    cross_correlation_window: float
    min_waveform_cc: float


@dataclass
class IterationConfig:
    """反演迭代控制参数。"""
    num_iterations: int
    max_perturbation: float
    damping_factor: float
    smoothing_factor: float


@dataclass
class RuntimeConfig:
    """程序运行相关参数。"""
    inversion_method: int
    python_command: str
    num_cores: int


@dataclass
class ForwardConfig3D:
    """三维正演总配置。"""
    Num: int
    Dt: float
    Time: float
    Nstep: int
    Nstep: int
    Domain: int
    Vp: float
    Vs: float
    Rho: float
    Nys: int



@dataclass
class InversionConfig3D:
    x1: float
    x2: float
    y1: float
    y2: float
    z1: float
    z2: float

    nx_inv: int
    ny_inv: int
    nz_inv: int

    nx_for: int
    ny_for: int
    nz_for: int

    num_sources: int
    num_receivers: int

    integration_interval: int
    wavelength: float

    cross_correlation_window: float
    min_waveform_cc: float

    num_iterations: int
    max_perturbation: float

    damping_factor: float
    smoothing_factor: float

    inversion_method: int
    python_command: str
    num_cores: int

######################################################
############# read configuration file ################
######################################################
def read_inversion_config3d(filename: str | Path) -> InversionConfig3D:
    filename = Path(filename)

    if not filename.is_file():
        raise FileNotFoundError(f"cannot open config file: {filename}")

    rows = []

    with filename.open("r", encoding="utf-8") as file:
        for line_number, line in enumerate(file, start=1):
            value_text = line.split("c:", 1)[0].strip()

            if not value_text:
                continue

            values = value_text.split()
            rows.append((line_number, values))

    expected_counts = [
        2,  # x1 x2
        2,  # y1 y2
        2,  # z1 z2
        3,  # nx_inv ny_inv nz_inv
        3,  # nx_for ny_for nz_for
        1,  # num_sources
        1,  # num_receivers
        1,  # integration_interval
        1,  # wavelength
        1,  # cross_correlation_window
        1,  # min_waveform_cc
        1,  # num_iterations
        1,  # max_perturbation
        1,  # damping_factor
        1,  # smoothing_factor
        1,  # inversion_method
        1,  # python_command
        1,  # num_cores
    ]

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

    try:
        config = InversionConfig3D(
            x1=float(rows[0][1][0]),
            x2=float(rows[0][1][1]),

            y1=float(rows[1][1][0]),
            y2=float(rows[1][1][1]),

            z1=float(rows[2][1][0]),
            z2=float(rows[2][1][1]),

            nx_inv=int(rows[3][1][0]),
            ny_inv=int(rows[3][1][1]),
            nz_inv=int(rows[3][1][2]),

            nx_for=int(rows[4][1][0]),
            ny_for=int(rows[4][1][1]),
            nz_for=int(rows[4][1][2]),

            num_sources=int(rows[5][1][0]),
            num_receivers=int(rows[6][1][0]),

            integration_interval=int(rows[7][1][0]),
            wavelength=float(rows[8][1][0]),

            cross_correlation_window=float(rows[9][1][0]),
            min_waveform_cc=float(rows[10][1][0]),

            num_iterations=int(rows[11][1][0]),
            max_perturbation=float(rows[12][1][0]),

            damping_factor=float(rows[13][1][0]),
            smoothing_factor=float(rows[14][1][0]),

            inversion_method=int(rows[15][1][0]),
            python_command=rows[16][1][0],
            num_cores=int(rows[17][1][0]),
        )

    except ValueError as exc:
        raise ValueError(
            f"invalid parameter type in config file: {filename}"
        ) from exc

    validate_inversion_config3d(config)

    return config


def validate_inversion_config3d(config: InversionConfig3D) -> None:
    if config.x2 <= config.x1:
        raise ValueError("x2 must be greater than x1")

    if config.y2 <= config.y1:
        raise ValueError("y2 must be greater than y1")

    if config.z2 <= config.z1:
        raise ValueError("z2 must be greater than z1")

    if config.nx_inv < 2 or config.ny_inv < 2 or config.nz_inv < 2:
        raise ValueError("nx_inv, ny_inv and nz_inv must all be >= 2")

    if config.nx_for < 2 or config.ny_for < 2 or config.nz_for < 2:
        raise ValueError("nx_for, ny_for and nz_for must all be >= 2")

    if config.num_sources < 1:
        raise ValueError("num_sources must be >= 1")

    if config.num_receivers < 1:
        raise ValueError("num_receivers must be >= 1")

    if config.integration_interval < 1:
        raise ValueError("integration_interval must be >= 1")

    if config.wavelength <= 0.0:
        raise ValueError("wavelength must be greater than 0")

    if config.cross_correlation_window <= 0.0:
        raise ValueError("cross_correlation_window must be greater than 0")

    if not 0.0 <= config.min_waveform_cc <= 1.0:
        raise ValueError("min_waveform_cc must be within [0, 1]")

    if config.num_iterations < 1:
        raise ValueError("num_iterations must be >= 1")

    if config.max_perturbation < 0.0:
        raise ValueError("max_perturbation must be >= 0")

    if config.damping_factor < 0.0:
        raise ValueError("damping_factor must be >= 0")

    if config.smoothing_factor < 0.0:
        raise ValueError("smoothing_factor must be >= 0")

    if config.inversion_method not in (1, 2, 3, 4):
        raise ValueError("inversion_method must be 1, 2, 3 or 4")

    if not config.python_command:
        raise ValueError("python_command cannot be empty")

    if config.num_cores < 1:
        raise ValueError("num_cores must be >= 1")


def print_inversion_config3d(config: InversionConfig3D) -> None:
    dx_inv = (config.x2 - config.x1) / (config.nx_inv - 1)
    dy_inv = (config.y2 - config.y1) / (config.ny_inv - 1)
    dz_inv = (config.z2 - config.z1) / (config.nz_inv - 1)

    dx_for = (config.x2 - config.x1) / (config.nx_for - 1)
    dy_for = (config.y2 - config.y1) / (config.ny_for - 1)
    dz_for = (config.z2 - config.z1) / (config.nz_for - 1)

    num_model_parameters = config.nx_inv * config.ny_inv * config.nz_inv

    print(f"X range                  : {config.x1:.6f} - {config.x2:.6f}")
    print(f"Y range                  : {config.y1:.6f} - {config.y2:.6f}")
    print(f"Z range                  : {config.z1:.6f} - {config.z2:.6f}")

    print(f"Grid nodes               : {config.nx_inv} x {config.ny_inv} x {config.nz_inv}")
    print(f"Grid nodes               : {config.nx_for} x {config.ny_for} x {config.nz_for}")
    print(f"Grid spacing             : dx_inv={dx_inv:.6f}, dy_inv={dy_inv:.6f}, dz_inv={dz_inv:.6f}")
    print(f"Grid spacing             : dx_for={dx_for:.6f}, dy_for={dy_for:.6f}, dz_for={dz_for:.6f}")
    print(f"Model parameters         : {num_model_parameters}")

    print(f"Number of sources        : {config.num_sources}")
    print(f"Number of receivers      : {config.num_receivers}")

    print(f"Integration interval     : {config.integration_interval}")
    print(f"Dominant wavelength      : {config.wavelength:.6f}")

    print(f"Cross-correlation window : {config.cross_correlation_window:.6f}")
    print(f"Minimum waveform CC      : {config.min_waveform_cc:.6f}")

    print(f"Maximum iterations       : {config.num_iterations}")
    print(f"Maximum perturbation     : {config.max_perturbation:.6e}")

    print(f"Damping factor           : {config.damping_factor:.6e}")
    print(f"Smoothing factor         : {config.smoothing_factor:.6e}")

    print(f"Inversion method         : {config.inversion_method}")
    print(f"Python command           : {config.python_command}")
    print(f"CPU cores                : {config.num_cores}")

######################################################
############# read initial model file ################
######################################################
def read_initial_model(filename):
    """
    读取输入速度模型。

    文件格式：
        x  y  z  vp  vs  rho
    """
    filename = Path(filename)

    if not filename.is_file():
        raise FileNotFoundError(f"cannot open initial model file: {filename}")

    data = np.loadtxt(filename, dtype=np.float64, ndmin=2)

    if data.shape[1] < 6:
        raise ValueError(
            f"{filename} must contain at least 6 columns: "
            "x y z vp vs rho"
        )

    data = data[:, :6]

    if not np.all(np.isfinite(data)):
        raise ValueError(f"{filename} contains NaN or Inf")

    return data

######################################################
############# read perturbation file #################
######################################################
def read_perturbation(filename):
    """
    读取反演模型扰动。

    文件格式：
        x  y  z  perturbation
    """
    filename = Path(filename)

    if not filename.is_file():
        raise FileNotFoundError(f"cannot open perturbation file: {filename}")

    data = np.loadtxt(filename, dtype=np.float64, ndmin=2)

    if data.shape[1] < 4:
        raise ValueError(
            f"{filename} must contain at least 4 columns: "
            "x y z perturbation"
        )

    data = data[:, :4]

    if not np.all(np.isfinite(data)):
        raise ValueError(f"{filename} contains NaN or Inf")

    return data


######################################################
############# read forward config file ###############
######################################################
def parse_value(text):
    """把字符串参数自动转换为 int、float 或 str。"""
    text = text.strip()

    # 尝试转成整数
    try:
        value = int(text)
        return value
    except ValueError:
        pass

    # 尝试转成浮点数
    try:
        value = float(text)
        return value
    except ValueError:
        pass

    # 不能转换就保留字符串
    return text
######################################################

def read_forward_config(config_file):
    """
    读取 parameters_forward.config 配置文件。

    文件格式示例：
        4        c: Interportation order for the spectral-element method //N
        0.01     c: Time interval //dt
        100000.0 40000.0 c: The horizontal and vertical scales ...

    返回：
        params: dict
    """
    config_file = Path(config_file)
    params = {}

    with config_file.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()

            # 跳过空行
            if not line:
                continue

            # 去掉注释前后的内容
            if "c:" in line:
                value_part, comment_part = line.split("c:", 1)
            else:
                value_part = line
                comment_part = ""

            value_part = value_part.strip()
            comment_part = comment_part.strip()

            # 提取变量名：优先使用 // 后面的内容
            if "//" in comment_part:
                _, name_part = comment_part.split("//", 1)
                name = name_part.strip()
            elif "horizontal and vertical scales" in comment_part.lower():
                # 兼容旧正演配置：末行模型长度、深度没有 // 变量名。
                name = "model_size"
            else:
                # 其他没有变量名的行保持原来的跳过行为。
                continue

            # 去掉变量名中的空格，例如 "// F" -> "F"
            name = name.replace(" ", "")

            # 读取参数值
            values = value_part.split()

            if len(values) == 1:
                params[name] = parse_value(values[0])
            else:
                params[name] = [parse_value(v) for v in values]

    return params