import numpy as np

def is_point_in_triangle(p, a, b, c):
    """
    判断点 p 是否在三角形 (a, b, c) 内（向量叉积法）
    """
    v0 = np.array(c) - np.array(a)
    v1 = np.array(b) - np.array(a)
    v2 = np.array(p) - np.array(a)
    
    dot00 = np.dot(v0, v0)
    dot01 = np.dot(v0, v1)
    dot02 = np.dot(v0, v2)
    dot11 = np.dot(v1, v1)
    dot12 = np.dot(v1, v2)
    
    inv_denom = 1.0 / (dot00 * dot11 - dot01 * dot01)
    u = (dot11 * dot02 - dot01 * dot12) * inv_denom
    v = (dot00 * dot12 - dot01 * dot02) * inv_denom
    
    return (u >= 0) & (v >= 0) & (u + v <= 1)

def barycentric_interpolation(p, a, b, c, val_a, val_b, val_c):
    """
    在三角形 (a, b, c) 内，对点 p 进行重心坐标插值
    """
    # 计算重心坐标（这里简化，实际可用面积比例法或上述叉积法）
    # 这里直接调用 is_point_in_triangle 的计算结果（实际需要单独计算权重）
    # 更合理的实现：
    v0 = np.array(c) - np.array(a)
    v1 = np.array(b) - np.array(a)
    v2 = np.array(p) - np.array(a)
    
    dot00 = np.dot(v0, v0)
    dot01 = np.dot(v0, v1)
    dot02 = np.dot(v0, v2)
    dot11 = np.dot(v1, v1)
    dot12 = np.dot(v1, v2)
    
    inv_denom = 1.0 / (dot00 * dot11 - dot01 * dot01)
    u = (dot11 * dot02 - dot01 * dot12) * inv_denom
    v = (dot00 * dot12 - dot01 * dot02) * inv_denom
    w = 1 - u - v
    
    return w * val_a + u * val_b + v * val_c

def interpolate_quad(p, v1, v2, v3, v4, val1, val2, val3, val4):
    """
    在任意四边形内插值（分割为两个三角形）
    """
    # 分割为三角形1 (v1, v2, v3) 和 三角形2 (v1, v3, v4)
    if is_point_in_triangle(p, v1, v2, v3):
        return barycentric_interpolation(p, v1, v2, v3, val1, val2, val3)
    elif is_point_in_triangle(p, v1, v3, v4):
        return barycentric_interpolation(p, v1, v3, v4, val1, val3, val4)
    else:
        # 点不在四边形内（理论上对凸四边形不应发生）
        return None

if __name__ == "__main__":
    # 示例使用
    v1 = (0, 0)
    v2 = (2, 2)
    v3 = (2, 2)
    v4 = (0, 2)

    val1 = 10
    val2 = 20
    val3 = 30
    val4 = 40

    p = (1, 1)
    result = interpolate_quad(p, v1, v2, v3, v4, val1, val2, val3, val4)
    print(f"插值结果: {result}")

