#include <stdio.h>
#include <math.h>
#include <stdbool.h>  // 添加必要的头文件

// 四边形顶点结构
typedef struct {
    double x;
    double y;
} Vec2;

int Nearest_point2d(double *target, double **points, int count);
int Nearest_point2d_lu(double *target, double **points, int count);
int Nearest_point2d_ru(double *target, double **points, int count);
int Nearest_point2d_ld(double *target, double **points, int count);
int Nearest_point2d_rd(double *target, double **points, int count);

double Cross(Vec2 a, Vec2 b);
Vec2 Invbilinear(Vec2 p, Vec2 a, Vec2 b, Vec2 c, Vec2 d) ;
double Interplot_2din(Vec2 p, Vec2 a, Vec2 b, Vec2 c, Vec2 d, 
    double fa, double fb, double fc, double fd);


