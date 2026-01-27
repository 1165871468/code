#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#define pi acos(-1.0)
#define earth_r 6371393.0

int *Arrayint1(long nx);
int **Arrayint2(long nx, long ny);
int Count_file(const char *flnm);
int File_exist(const char* flnm);
double *Arraydouble1(long nx);
double **Arraydouble2(long nx, long ny);
double ***Arraydouble3(long N1,long N2,long N3);
FILE *Open_file(const char *flnm, const char *mode);
void Freeint2(long nx, int **Array);
void Average_rank(int rank, int size, int Ns, int *t1, int *t2);
void Freedouble2(long nx, double **Array);
void Freedouble3(long N1, long N2, double ***Array);
void Get_model_infor(FILE *fp1, int num_s, int num_r, double loc_s[], double loc_r[], double loc_new[],
    double rot_theta[], double y_o, double model_dis, double model_dep);
void Model_parameter(int *Ns, int *Nr, double *elong_factor);
void Read_file(char *flnm, int row, int col, double **array_data);
