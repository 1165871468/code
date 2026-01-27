#ifndef _Interpolation_h
#define _Interpolation_h

typedef struct
{
    double x;
    double y;
} Vec2;

int Get_num_xz(int count, double x1, double x2, double xx);
int Nearest_point1d(double target, double *points, int count);
int Nearest_point2d(double *target, double **points, int count);
int Nearest_point3d(double *target, double **points, int count);
int Nearest_point2d_lu(double *target, double **points, int count);
int Nearest_point2d_ru(double *target, double **points, int count);
int Nearest_point2d_ld(double *target, double **points, int count);
int Nearest_point2d_rd(double *target, double **points, int count);

double Cross(Vec2 a, Vec2 b);
double Interplot_1d(double xx[3], double yy[3]);
double Interplot_2d(double xx[3], double yy[3], double zz[2][2]);
double Interplot_2din(Vec2 p, Vec2 a, Vec2 b, Vec2 c, Vec2 d, double fa, double fb, double fc, double fd);
double Interplot_3d(double xx[3], double yy[3], double zz[3], double vv[2][2][2]);

void Find_node(int numData, double **arrayData, double *tarPoints, double arrayPoints[][3]);
void Fine_node_cube(int numY, int numZ, int idx1, int idy1, int idz1, int idx2, int idy2, int idz2, 
    double **arrayV, double xx[3], double yy[3], double zz[3], double vv[2][2][2]);
int Fine_node_square(int numData, int id[4], double **arrayData, double *locNode);

Vec2 Invbilinear(Vec2 p, Vec2 a, Vec2 b, Vec2 c, Vec2 d);

#endif
