#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include <unistd.h>

int *Arrayint1(long Nx);
int Count_file(char *FL);
long *Arraylong1(long Nx);
float *Arrayfloat1(long Nx);
double *Arraydouble1(long Nx);
double **Arraydouble2(long Nx,long Nz);
double ***Arraydouble3(long N1,long N2,long N3);
double Lagrange_value(int n, int i, double x, double *data_x);
double Interpolate(int n, double x, double y, double *arrfx, double *arrfz, double **arrfy);
void Rdfile(char *FL, int count, int n, double **ary);
void Freedouble2(long Nx, double **Array);
void Freedouble3(long N1, long N2, double ***Array);
void Matrix_infor(int *numx, int *numz, double *x0, double *x1, double *z0, double *z1);

int main(int argc, char *argv[])
{
    int i, j, k, l, n, count1, count2, nmx, nmz, flgi1, fs;
	int num_cor, Ns, Nr, row, col, nz, step, step_max;

    char flnm1[100], flnm2[100], flnm3[100], flnm4[100], flnm5[100], FF[250];

    double xx1, xx2, x0, x1, zz1, zz2, z0, z1, dx, dz, dis, rr, eta, ff[3];
    double **arrdm1;

    FILE *fp1, *fp2, *fp3;

	sprintf(flnm1, "../Output/Model/Forver/vel_initial.dat");
	sprintf(flnm2, "../Output/Model/Forver/vel.dat");
	sprintf(flnm3, "../Output/Model/Forver/per_vel.dat");
	sprintf(flnm4, "../Output/Model/Inver/per_vel.dat");
	sprintf(flnm5, "../Output/Model/Inver/vel_initial.dat");
	sprintf(FF, " cp %s %s", flnm1, flnm2);

	fs = system(FF);

	fp1 = fopen(flnm3, "w"); fp2 = fopen(flnm4, "w"); fp3 = fopen(flnm5, "w");
	if (fp1 == NULL || fp2 == NULL || fp3 == NULL)
	{
		printf("Open file filed! -- Initialize...\n"), exit(1);
	}

	Matrix_infor(&nmx, &nmz, &x0, &x1, &z0, &z1);

	dx = (x1 - x0) / (nmx - 1);
	dz = (z1 - z0) / (nmz - 1);
	dis = dx;

	count2 = Count_file(flnm1);
	arrdm1 = Arraydouble2(count2, 5);

	Rdfile(flnm1, count2, 5, arrdm1);

	for (i = 0; i < count2; i++)
	{
		xx1 = arrdm1[i][0];
		zz1 = arrdm1[i][1];
		fprintf(fp1, "%.1f%10.1f%10.1f%10.1f%10.1f\n", xx1, zz1, 0.0, 0.0, 0.0), fflush(fp1);
	}

	for (i = 0; i < nmx; i++)
	{
		for (j = 0; j < nmz; j++)
		{
			xx1 = x0 + i*dx;
			zz1 = z0 + j*dz;
			fprintf(fp2, "%10.1f%10.1f%10.1f%10.1f%10.1f\n", xx1, zz1, 0.0, 0.0, 0.0), fflush(fp2);

			for (k = 0; k < count2; k++)
			{
				xx2 =  arrdm1[k][0];
				zz2 =  arrdm1[k][1];

				rr = sqrt((xx1 - xx2)*(xx1 - xx2) + (zz1 - zz2)*(zz1 - zz2));

				if (rr <= dis)
				{
					ff[0] = arrdm1[k][2];
					ff[1] = arrdm1[k][3];
					ff[2] = arrdm1[k][4];
					flgi1 = 1;
				}
			}

			if (flgi1 == 0)
			{
				printf("This node can't be find %10.1f%10.1f\n", xx1, zz1);
			}
			else
			{
				fprintf(fp3, "%10.1f%10.1f%10.1f%10.1f%10.1f\n", xx1, zz1, ff[0], ff[1], ff[2]), fflush(fp3);
			}
		}
	}

//    printf("node number: %10d%10d%10d\n", nmx, nmz, count);
////////////////////////////////////////////////////////////////
	Freedouble2(count2, arrdm1);
	fclose(fp1); fclose(fp2); fclose(fp3);
	return 0;
}

double Interpolate(int n, double x, double y, double *arrfx, double *arrfz, double **arrfy)
{
    int i, j;
    double shapx, shapz, result;
    result = 0.0;

	for (i = 0; i < n - 1; i++)
	{
		if (arrfx[i] == arrfx[i + 1])
			arrfx[i + 1] = arrfx[i] + 0.000001;

		if (arrfz[i] == arrfz[i + 1])
			arrfz[i + 1] = arrfz[i] + 0.000001;
	}
	
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            shapx = Lagrange_value(n, i, x, arrfx);
            shapz = Lagrange_value(n, j, y, arrfz);

            result += shapx * shapz * arrfy[i][j];
        }
    }

    return result;
}

double Lagrange_value(int n, int i, double x, double *data_x)
{
    int j;
    double term;

    term = 1.0;

    for (j = 0; j < n; j++)
    {
        if (i != j)
        {
            term *= (x - data_x[j]) / (data_x[i] - data_x[j]);
        }
    }

    return term;
}

int Count_file(char *FL)
{
	int count;
	char FF[100],*ff;
	FILE *fp;

	fp = fopen(FL,"r");

	if (fp == NULL)
		printf("can't open file --- %s\n", FL), exit(1);

	count = 0;

	while((ff = fgets(FF, 100, fp))!=NULL)
	{
		count++;
	}

	rewind(fp), fclose(fp);

	return count;
}

void Rdfile(char *FL, int count, int n, double **ary)
{
	int i, fs;
	char FF[100],*ff;
	float tmp1,tmp2,tmp3,tmp4,tmp5;
	FILE *fp;

	fp = fopen(FL,"r");
	
	if (fp == NULL)
		printf("can't open file --- %s\n", FL), exit(1);
	
	
	if (n == 1)
	{
		for ( i = 0; i < count; i++)
		{
			fs = fscanf(fp,"%f",&tmp1);
			ary[i][0]=tmp1;
		}
	}

	if (n == 2)
	{
		for ( i = 0; i < count; i++)
		{
			fs = fscanf(fp,"%f%f",&tmp1,&tmp2);
			ary[i][0]=tmp1;
			ary[i][1]=tmp2;
		}
	}
	
	if (n == 3)
	{
		for ( i = 0; i < count; i++)
		{
			fs = fscanf(fp,"%f%f%f",&tmp1,&tmp2,&tmp3);
			ary[i][0]=tmp1;
			ary[i][1]=tmp2;
			ary[i][2]=tmp3;
		}
	}

	if (n == 4)
	{
		for ( i = 0; i < count; i++)
		{
			fs = fscanf(fp,"%f%f%f%f",&tmp1,&tmp2,&tmp3,&tmp4);
			ary[i][0]=tmp1;
			ary[i][1]=tmp2;
			ary[i][2]=tmp3;
			ary[i][3]=tmp4;
		}
	}

	if (n == 5)
	{
		for ( i = 0; i < count; i++)
		{
			fs = fscanf(fp,"%f%f%f%f%f",&tmp1,&tmp2,&tmp3,&tmp4,&tmp5);
			ary[i][0]=tmp1;
			ary[i][1]=tmp2;
			ary[i][2]=tmp3;
			ary[i][3]=tmp4;
			ary[i][4]=tmp5;
		}
	}

	rewind(fp), fclose(fp);
}

int *Arrayint1(long Nx)
{
	long i;
	int *Array;
	Array=(int *)malloc(Nx*sizeof(int));
	for(i=0;i<Nx;i++)
		Array[i]=0;
	return Array;
}

double *Arraydouble1(long Nx)
{
	long i;
	double *Array;
	Array=(double *)malloc(Nx*sizeof(double));
	for(i=0;i<Nx;i++)
		Array[i]=0.0;
	return Array;
}

float *Arrayfloat1(long Nx)
{
	long i;
	float *Array;
	Array=(float *)malloc(Nx*sizeof(float));
	for(i=0;i<Nx;i++)
		Array[i]=0.0;
	return Array;
}

long *Arraylong1(long Nx)
{
	long i;
	long *Array;
	Array=(long *)malloc(Nx*sizeof(long));
	for(i=0;i<Nx;i++)
		Array[i]=0;
	return Array;
}

double **Arraydouble2(long Nx,long Nz)
{
	long i,j;
	double **Array;
	Array=(double **)malloc(Nx*sizeof(double *));
	for(i=0;i<Nx;i++)
		*(Array+i)=(double *)malloc(Nz*sizeof(double));

	for(i=0;i<Nx;i++)
		for(j=0;j<Nz;j++)
			Array[i][j]=0.0;

	return Array;
}

double ***Arraydouble3(long N1,long N2,long N3)
{
	long i,j,k;
	double ***Array;

	Array=(double ***)malloc(N1*sizeof(double **));

	for(i=0;i<N1;i++)
	{
		Array[i]=(double **)malloc(N2*sizeof(double *));
		for(j=0;j<N2;j++)
		{
			Array[i][j]=(double *)malloc(N3*sizeof(double *));
			for(k=0;k<N3;k++)
			{
				Array[i][j][k]=0.0;
			}
		}
	}

	return Array;
}

void Freedouble2(long Nx, double **Array)
{
	long i;
	for(i=0;i<Nx;i++)
		free(Array[i]);
	free(Array);
}

void Freedouble3(long N1, long N2, double ***Array)
{
	long i, j;
	for(i = 0; i < N1; i++)
	{
		for (j = 0; j < N2; j++)
		{
			free(Array[i][j]);
		}
		
		free(Array[i]);
	}

	free(Array);
}

void Matrix_infor(int *numx, int *numz, double *x0, double *x1, double *z0, double *z1)
{
	int i, j, ii1, ii2, sf1;
	double dd1, dd2;
	char cc1[500], flnm1[100], flnm2[100], *sf2;
	FILE *fp1;

	sprintf(flnm1, "../Input/parameters_inversion");

	if ((fp1 = fopen(flnm1, "r")) == NULL)
	{
		printf("Can not read this file -- %s", flnm1); exit(1);
	}
	////////////////////////////////////////////////////////////////////////////
	sf1 = fscanf(fp1, "%lf%lf", &dd1, &dd2);  sf2 = fgets(cc1, 500, fp1); // 1
	*x0 = dd1; *x1 = dd2;
	sf1 = fscanf(fp1, "%lf%lf", &dd1, &dd2);  sf2 = fgets(cc1, 500, fp1); // 2
	*z0 = dd1; *z1 = dd2;
	sf1 = fscanf(fp1, "%d%d", &ii1, &ii2);  sf2 = fgets(cc1, 500, fp1);   // 3
	*numx = ii1; *numz = ii2;
}
