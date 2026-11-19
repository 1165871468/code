#ifndef _common_h
#define _common_h
///////////////////////////////////////////////////////////////////////
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include <unistd.h>
///////////////////////////////////////////////////////////////////////
int Count_file(char *FL);
int *Arrayint1(long Nx);
long *Arraylong1(long Nx);
float *Arrayfloat1(long Nx);

double *Arraydouble1(long Nx);
double **Arraydouble2(long Nx,long Nz);
double ***Arraydouble3(long N1,long N2,long N3);
double Lagrange_value(int n, int i, double x, double *data_x);
double Interpolate(int n, double x, double y, double *arrfx, double *arrfz, double **arrfy);
double Vector_product(long len, double *V1, double *V2);

void Rdfile(char *FL, int count, int n, double **ary);
void Matrix_multi(double **A, int m1, int n, double **B, int m2, double **C);
void Freedouble2(long Nx, double **Array);
void Freedouble3(long N1, long N2, double ***Array);


///////////////////////////////////////////////////////////////////////
#endif