#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include <unistd.h>

int countfile(char *FL);
int *Arrayint1(long Nx);
long *Arraylong1(long Nx);
float *Arrayfloat1(long Nx);
double *Arraydouble1(long Nx);
double **Arraydouble2(long Nx,long Nz);
double ***Arraydouble3(long N1,long N2,long N3);
double lagrange_value(int n, int i, double x, double *data_x);
double interpolate(int n, double x, double y, double *arrfx, double *arrfz, double **arrfy);
void rdfile(char *FL, int count, int n, double **ary);
void freedouble2(long Nx, double **Array);
void freedouble3(long N1, long N2, double ***Array);

int main(int argc, char *argv[])
{
    int i, j, k, l, n, count;

    double a, b, dd, rr1, rr2, xx1, xx2, zz1, zz2, ff1, ff2, pi = acos(-1.0);

    double **arrdm1, **arrdm2, **arrdm3, **arrdm4, **arrdm5;

    char flnm1[100], flnm2[100], flnm3[100], flnm4[100];

    FILE *fp1;

	dd = 2000.0;

	sprintf(flnm1, "../S1/Kernel/Kernel_lam.dat");
	sprintf(flnm2, "../S2/Kernel/Kernel_lam.dat");
	sprintf(flnm3, "../S3/Kernel/Kernel_lam.dat");
	sprintf(flnm4, "Kernel_lam_smooth.dat");

	fp1 = fopen(flnm4, "w");
	if (fp1 == NULL) {
		printf("Can't open file-- %s !\n", flnm4), exit(1);
	}

	count = countfile(flnm1);

	arrdm1 = Arraydouble2(count, 3);
	arrdm2 = Arraydouble2(count, 3);
	arrdm3 = Arraydouble2(count, 3);
	arrdm4 = Arraydouble2(count, 3);
	arrdm5 = Arraydouble2(count, 3);

	rdfile(flnm1, count, 3, arrdm1);
	rdfile(flnm2, count, 3, arrdm2);
	rdfile(flnm3, count, 3, arrdm3);


	for (i = 0; i < count; i++) {
		arrdm4[i][0] = arrdm1[i][0];
		arrdm4[i][1] = arrdm1[i][1];
		arrdm4[i][2] = arrdm1[i][2] + arrdm2[i][2] + arrdm3[i][2];

		arrdm5[i][0] = arrdm4[i][0];
		arrdm5[i][1] = arrdm4[i][1];
		arrdm5[i][2] = arrdm4[i][2];
	}

	a = 4.0/(pi*dd*dd);

	for (i = 0; i < count; i++) {
		xx1 = arrdm5[i][0];
		zz1 = arrdm5[i][1];

		ff1 = 0.0;
		ff2 = 0.0;

		for (j = 0; j < count; j++) {
			xx2 = arrdm4[j][0];
			zz2 = arrdm4[j][1];

			rr1 = fabs(xx2 - xx1) / dd;
			rr2 = fabs(zz2 - zz1) / dd;

			if (rr1 >= 3.0 && rr2 >= 3.0) {
				continue;
			}
			
			b = -4.0*((xx2 - xx1)*(xx2 - xx1) + (zz2 - zz1)*(zz2 - zz1)) / (dd*dd);

			ff1 += arrdm4[j][2]*a*exp(b);
			ff2 += a*exp(b);
		}

		ff1 = ff1 / ff2;

		fprintf(fp1, "%.1f%10.1f%20.8e\n", arrdm5[i][0], arrdm5[i][1], ff1), fflush(fp1);
	}

////////////////////////////////////////////////////////////////
	freedouble2(count, arrdm1);
	freedouble2(count, arrdm2);
	freedouble2(count, arrdm3);
	freedouble2(count, arrdm4);
	freedouble2(count, arrdm5);
	fclose(fp1);
	return 0;
}

double interpolate(int n, double x, double y, double *arrfx, double *arrfz, double **arrfy)
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
            shapx = lagrange_value(n, i, x, arrfx);
            shapz = lagrange_value(n, j, y, arrfz);

            result += shapx * shapz * arrfy[i][j];
        }
    }

    return result;
}

double lagrange_value(int n, int i, double x, double *data_x)
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

void Find_invernode(int nmx, int nmz, long elem_num,double len,double **FDM_array_x, double **FDM_array_z,
	long *Elem1_count,long **Elem_array,double **node_array,double **Fii_fdm,double **Fax_fdm,double **Faz_fdm, double **Fii_sem)
{
/*
	int i, j, k, l, m, N, it, i1, i2, i3, i4, id, flgi, flgi1, flgi2, Cin[4];
	
	int *arr_int;

	float a1, a2, a3, a4, fx, fz, x1, x2, x3, x4, xx, zz, z1, z2, z3, z4, ds, tmp;

	double *P, c1, c2;

	N = 200 - 1;

	ds = sqrt(2) * len / (N + 1);

	arr_int = Arrayint1(2);
	P = Arraydouble1(N);

	Cin[0] = 0;
	Cin[1] = num * (num - 1);
	Cin[2] = num * num - 1;
	Cin[3] = num - 1;

	for(i = 0; i < N; i++)
	{
		P[i] = -1.0 + 2.0 / (N-1) * i;
	}
	
	for(i = L; i <= Nx + L; i++)
	{		
		for(j = 0; j <= zf; j++)
		{
			fx = FDM_array_x[i][j];
			fz = FDM_array_z[i][j];

			flgi = 0;

			for(k = 0; k < elem_num; k++)
			{
				it = Elem1_count[k] - 1;
				i1 = Elem_array[it][Cin[0]] - 1;
				x1 = node_array[i1][0];
				z1 = node_array[i1][1];

				i2 = Elem_array[it][Cin[1]]-1;
				x2 = node_array[i2][0];
				z2 = node_array[i2][1];

				i3 = Elem_array[it][Cin[2]]-1;
				x3 = node_array[i3][0];
				z3 = node_array[i3][1];

				i4 = Elem_array[it][Cin[3]]-1;
				x4 = node_array[i4][0];
				z4 = node_array[i4][1];

				xx = (x1 + x2 + x3 + x4) / 4.0 - fx;
				zz = (z1 + z2 + z3 + z4) / 4.0 - fz;

				if (fabs(xx) <= len && fabs(zz) <= len)
				{
					for(l = 0; l < N; l++)
					{
						for(m = 0; m < N ; m++)
						{
							a1 = (1.0 - P[l]) * (1.0 - P[m]) * 0.25;
							a2 = (1.0 + P[l]) * (1.0 - P[m]) * 0.25;
							a3 = (1.0 + P[l]) * (1.0 + P[m]) * 0.25;
							a4 = (1.0 - P[l]) * (1.0 + P[m]) * 0.25;
												
							xx = a1 * x1 + a2 * x2 + a3 * x3 + a4 * x4 - fx;
							zz = a1 * z1 + a2 * z2 + a3 * z3 + a4 * z4 - fz;
							// c1 = a1 * x1 + a2 * x2 + a3 * x3 + a4 * x4;
							// c2 = a1 * z1 + a2 * z2 + a3 * z3 + a4 * z4;

							tmp = sqrt(xx * xx + zz * zz);
						
							if(tmp <= ds)
							{
								Fii_fdm[i][j] = k;
								Fax_fdm[i][j] = P[l];
								Faz_fdm[i][j] = P[m];

								flgi = 1;
								break;
							}
						}

						if (flgi == 1)
							break;
					}

					if (flgi == 1)
						break;
				}
			}
		}
	}

	free(P), free(arr_int);
*/
}

int countfile(char *FL)
{
	int count;
	char FF[100],*ff;
	FILE *fp;

	fp = fopen(FL,"r");

	if (fp == NULL)
		printf("can't open file --- %s\n", FL), exit(1);

	count = 0;

	while((ff = fgets(FF, 200, fp))!=NULL)
	{
		count++;
	}

	rewind(fp), fclose(fp);

	return count;
}

void rdfile(char *FL, int count, int n, double **ary)
{
	int i;
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
			fscanf(fp,"%f",&tmp1);
			ary[i][0]=tmp1;
		}
	}

	if (n == 2)
	{
		for ( i = 0; i < count; i++)
		{
			fscanf(fp,"%f%f",&tmp1,&tmp2);
			ary[i][0]=tmp1;
			ary[i][1]=tmp2;
		}
	}
	
	if (n == 3)
	{
		for ( i = 0; i < count; i++)
		{
			fscanf(fp,"%f%f%f",&tmp1,&tmp2,&tmp3);
			ary[i][0]=tmp1;
			ary[i][1]=tmp2;
			ary[i][2]=tmp3;
		}
	}

	if (n == 4)
	{
		for ( i = 0; i < count; i++)
		{
			fscanf(fp,"%f%f%f%f",&tmp1,&tmp2,&tmp3,&tmp4);
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
			fscanf(fp,"%f%f%f%f%f",&tmp1,&tmp2,&tmp3,&tmp4,&tmp5);
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

void freedouble2(long Nx, double **Array)
{
	long i;
	for(i=0;i<Nx;i++)
		free(Array[i]);
	free(Array);
}

void freedouble3(long N1, long N2, double ***Array)
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
