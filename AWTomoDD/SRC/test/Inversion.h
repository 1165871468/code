#ifndef _CG_h
#define _CG_h

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

typedef struct
{
	int row;
	int col;
	float val;
}Sarray;


int countfile(char *FL);
int *Arrayint1(long Nx);
long *Arraylong1(long Nx);
float *Arrayfloat1(long Nx);

double *Arraydouble1(long Nx);
double **Arraydouble2(long Nx,long Nz);
double ***Arraydouble3(long N1,long N2,long N3);
double lagrange_value(int n, int i, double x, double *data_x);
double interpolate(int n, double x, double y, double *arrfx, double *arrfz, double **arrfy);
double Max_abs_vec(long len, double *vector0, double *vector);
double Vector_product(long len, double *V1, double *V2);

void CG(long nz, long row, long col, double eta, Sarray *Fd, Sarray *Ft, double *Unc, double *data, double *deltaX);
void Delta_update(long col, double vt, double *Pcg, double *deltaX, double *deltaX1);
void Rdfile(char *FL, int count, int n, double **ary);
void Initial_sarray(long nz, Sarray *Fd);
void Freedouble2(long Nx, double **Array);
void Freedouble3(long N1, long N2, double ***Array);
void Sensitivity_matrix(long row, long col, long nz, double **array_a, Sarray *Fd);
void Transpose(long row, long col, long nz, Sarray *Fd, Sarray *Ft);
void Gradient(long row, long col, long nz, double eta, double *Ta, Sarray *Fd, Sarray *Ft, double *data, double *Unc, double *deltaX, double *Gr);
void Total_error(long nz, long row, long col, double eta, Sarray *Fd, double *data, double *deltaX, double *Unc, double *error);
void Error0(long row, double *data, double *error);
void Model_update(int row, int col, int nz, double *deltaX, double *data, Sarray *Fd);


#endif