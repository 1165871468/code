#include <stdio.h>
#include <math.h>
#include <stdbool.h>  // 添加必要的头文件

// 四边形顶点结构
typedef struct {
    double x;
    double y;
} Vec2;
double Cross(Vec2 a, Vec2 b);
Vec2 Invbilinear(Vec2 p, Vec2 a, Vec2 b, Vec2 c, Vec2 d) ;
double Interplot_2din(Vec2 p, Vec2 a, Vec2 b, Vec2 c, Vec2 d, 
    double fa, double fb, double fc, double fd);


int main()
{
    int i, j, k;
    double ff, fa, fb, fc, fd;
    Vec2 p, a, b, c, d;

    // a.x = 0.0;  a.y = 2.0;  fa = 1.0;
    // b.x = 1.0;  b.y = 1.0;  fb = 1.0;
    // c.x = 1.0;  c.y = 0.0;  fc = 0.0;
    // d.x = 0.0;  d.y = 0.0;  fd = 0.0;

    a.x = 1.0;  a.y = 2.0;  fa = 1.0;
    b.x = 1.0;  b.y = 0.0;  fb = 0.0;
    c.x = 0.0;  c.y = 0.0;  fc = 0.0;
    d.x = 0.0;  d.y = 1.0;  fd = 1.0;

    p.x = 0.5;  p.y = 0.5;


    ff = Interplot_2din(p, a, b, c, d, fa, fb, fc, fd);


    printf("%.4f \n", ff);
}



// 计算二维向量叉积 (a.x*b.y - a.y*b.x)
double Cross(Vec2 a, Vec2 b)
{
    return a.x * b.y - a.y * b.x;
}

// 计算逆双线性变换，返回(u, v)参数
Vec2 Invbilinear(Vec2 p, Vec2 a, Vec2 b, Vec2 c, Vec2 d)
{
    Vec2 res = {-1.0, -1.0};
    Vec2 e = {b.x - a.x, b.y - a.y};       // 边向量a->b
    Vec2 f = {d.x - a.x, d.y - a.y};       // 边向量a->d
    Vec2 g = {c.x - b.x - d.x + a.x, c.y - b.y - d.y + a.y};  // 扭曲向量
    Vec2 h = {p.x - a.x, p.y - a.y};       // 点p相对于a的偏移

    double k2 = Cross(g, f);
    double k1 = Cross(e, f) + Cross(h, g);
    double k0 = Cross(h, e);

    // 处理接近线性的情况 (k2很小)
    const double EPS = 1e-9;
    if (fabs(k2) < EPS) {
        // 线性方程: k1*v + k0 = 0
        if (fabs(k1) < EPS) {
            return res;  // 无解
        }
        
        double v = -k0 / k1;
        // 计算u: h = e*u + f*v + g*u*v，当k2≈0时g≈0，简化为h ≈ e*u + f*v
        double denominator = e.x;
        if (fabs(denominator) < EPS) {
            denominator = e.y;  // 避免除以0，使用y分量
            if (fabs(denominator) < EPS) {
                return res;  // 退化边
            }
            double u = (h.y - f.y * v) / denominator;
            res.x = u;
            res.y = v;
        } else {
            double u = (h.x - f.x * v) / denominator;
            res.x = u;
            res.y = v;
        }
    } else {
        // 二次方程求解: k2*v² + k1*v + k0 = 0
        double discriminant = k1 * k1 - 4.0 * k0 * k2;
        if (discriminant < 0.0) {
            return res;  // 无实根
        }
        discriminant = sqrt(discriminant);
        double ik2 = 0.5 / k2;
        
        // 尝试第一个解
        double v = (-k1 - discriminant) * ik2;
        double denominator = e.x + g.x * v;
        double u;
        
        if (fabs(denominator) < EPS) {
            denominator = e.y + g.y * v;  // 使用y分量计算u
            if (fabs(denominator) < EPS) {
                return res;  // 退化情况
            }
            u = (h.y - f.y * v) / denominator;
        } else {
            u = (h.x - f.x * v) / denominator;
        }
        
        // 检查是否在[0,1]范围内，否则尝试第二个解
        if (u < -EPS || u > 1.0 + EPS || v < -EPS || v > 1.0 + EPS) {
            v = (-k1 + discriminant) * ik2;
            denominator = e.x + g.x * v;
            
            if (fabs(denominator) < EPS) {
                denominator = e.y + g.y * v;
                if (fabs(denominator) < EPS) {
                    return res;
                }
                u = (h.y - f.y * v) / denominator;
            } else {
                u = (h.x - f.x * v) / denominator;
            }
        }
        
        res.x = u;
        res.y = v;
    }
    
    return res;
}

// 二维四边形插值主函数
double Interplot_2din(Vec2 p, Vec2 a, Vec2 b, Vec2 c, Vec2 d, 
    double fa, double fb, double fc, double fd)
{
    Vec2 uv = Invbilinear(p, a, b, c, d);
    const double EPS = 1e-9;
    
    // 检查参数是否在有效范围内 [0, 1]
    if (uv.x < -EPS || uv.x > 1.0 + EPS || uv.y < -EPS || uv.y > 1.0 + EPS) {
        printf("Can not find this point! -- %.4f %.4f\n", p.x, p.y);  return -12345.0;
    }
    
    // 确保u和v在[0,1]范围内（处理浮点误差）
    double u = fmax(0.0, fmin(1.0, uv.x));
    double v = fmax(0.0, fmin(1.0, uv.y));
    
    // 双线性插值公式
    return (1.0 - u) * (1.0 - v) * fa +  // a点权重
           u * (1.0 - v) * fb +          // b点权重
           u * v * fc +                  // c点权重
           (1.0 - u) * v * fd;           // d点权重
}

    
