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
	double val;
}Sarray;


int Count_file(const char *FL);
int *Arrayint1(long Nx);
long *Arraylong1(long Nx);
float *Arrayfloat1(long Nx);

double *Arraydouble1(long Nx);
double **Arraydouble2(long Nx,long Nz);
double ***Arraydouble3(long N1,long N2,long N3);
double lagrange_value(int n, int i, double x, double *data_x);
double interpolate(int n, double x, double y, double *arrfx, double *arrfz, double **arrfy);
double interpolate2D(int ix, int jy, double x, double y, double ***arrayData);
double Max_abs_vec(long len, double *vector);
double Vector_product(long len, double *V1, double *V2);

FILE *Open_file(const char *flnm, const char *mode);

static int Read_config_rows(char rows[18][500]);
static void Model_grid_infor(int *nx, int *ny, int *nz, double lo[3], double hi[3]);
static void Read_python_command(char command[128]);

void CG(long nz, long row, long col, double eta, Sarray *Fd, Sarray *Ft, double *Unc, double *data, double *deltaX);
void Collect_matrix(int step, int Ns, int Nr, int row, int col, double *data, double **A);
void Delta_update(long col, double vt, double *Pcg, double *deltaX, double *deltaX1);
void Error0(long row, double *data, double *error);
void Freedouble2(long Nx, double **Array);
void Freedouble3(long N1, long N2, double ***Array);
void Gradient(long row, long col, long nz, double eta, double *Ta, Sarray *Fd, Sarray *Ft, double *data, double *Unc, double *deltaX, double *Gr);
void Get_gradient(int num);
void Initial_sarray(long nz, Sarray *Fd);
void Initial_velocity0(int col, double *velc);
void Initial_velocity1(int col, double *velc);
void Matrix_infor(int *num_cor, int *Ns, int *Nr, int *row, int *col, int *nz, int *step, double *eta);
void Model_update(int step, int col, double *deltaX, double *velc);
void Output_delta(int step, int nz, int row, int col, int Nz, double *deltaX);
void Rdfile(char *flnm, int row, int col, double **array_data);
void Sensitivity_matrix(int row, int col, int nz, double **array_a, Sarray *Fd);
void Transpose(int row, int col, int nz, Sarray *Fd, Sarray *Ft);
void Total_error(long nz, long row, long col, double eta, Sarray *Fd, double *data, double *deltaX, double *Unc, double *error);
void modInforIn(int *numx, int *numz, double *x0, double *x1, double *z0, double *z1);
void Update_matrix(int *row, int *col, int *nz);


#endif
