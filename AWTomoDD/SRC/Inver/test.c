#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
void Average_rank(int rank, int size, int Ns, int *t1, int *t2);

int main(int argc, char *argv[])
{
	clock_t start, finish;

	int i, j, k, l, n, id1, s1, s2, rank, size, Nr, Ns, thd, nmx, nmz, count1, count2, count3, count4, fs;
	double a, b, dd, dx, dz, rr1, rr2, xx1, xx2, zz1, zz2, ff1, ff2, tt, coff, pi = acos(-1.0);
	double **arrdm1, **arrdm2, **arrdm3, **arrdm4, **arrdm5;
	char flnm1[100], flnm2[100], flnm3[100], flnm4[100];
	FILE *fp1;

	Ns = 10;
	Nr = 10;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	Average_rank(rank, size, Nr, &s1, &s2);

	for (i = 0; i < Ns; i++)
	{
		if (rank == 0)
		{
			printf("%10.1f%%\n", 100.0*((float)(i+1)/Ns));
		}


		for (j = s1; j <= s2; j++)
		{    k = 1;
			// sleep(2);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == 0)
	{
		printf("yes\n");
	}

	MPI_Finalize();
	return 0;
}

int Count_file(char *FL)
{
	int count1;
	char FF[200],*ff;
	FILE *fp;

	fp = fopen(FL,"r");

	if (fp == NULL)
	{
		printf("can't open file --- %s\n", FL), exit(1);
	}


	count1 = 0;

	while((ff = fgets(FF, sizeof FF, fp))!=NULL)
	{
		count1++;
	}

	rewind(fp), fclose(fp);

	return count1;
}

void Rdfile(char *FL, int count1, int n, double **ary)
{
	int i, fs;
	char FF[100],*ff;
	float tmp1,tmp2,tmp3,tmp4,tmp5;
	FILE *fp;

	fp = fopen(FL,"r");
	
	if (fp == NULL)
	{
		printf("can't open file --- %s\n", FL), exit(1);
	}

	if (n == 1)
	{
		for ( i = 0; i < count1; i++)
		{
			fs = fscanf(fp,"%f",&tmp1);
			ary[i][0]=tmp1;
		}
	}

	if (n == 2)
	{
		for ( i = 0; i < count1; i++)
		{
			fs = fscanf(fp,"%f%f",&tmp1,&tmp2);
			ary[i][0]=tmp1;
			ary[i][1]=tmp2;
		}
	}
	
	if (n == 3)
	{
		for ( i = 0; i < count1; i++)
		{
			fs = fscanf(fp,"%f%f%f",&tmp1,&tmp2,&tmp3);
			ary[i][0]=tmp1;
			ary[i][1]=tmp2;
			ary[i][2]=tmp3;
		}
	}

	if (n == 4)
	{
		for ( i = 0; i < count1; i++)
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
		for ( i = 0; i < count1; i++)
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
	{
		Array[i]=0;
	}

	return Array;
}

double *Arraydouble1(long Nx)
{
	long i;
	double *Array;

	Array=(double *)malloc(Nx*sizeof(double));

	for(i=0;i<Nx;i++)
	{
		Array[i]=0.0;
	}

	return Array;
}

float *Arrayfloat1(long Nx)
{
	long i;
	float *Array;

	Array=(float *)malloc(Nx*sizeof(float));

	for(i=0;i<Nx;i++)
	{
		Array[i]=0.0;
	}

	return Array;
}

long *Arraylong1(long Nx)
{
	long i;
	long *Array;

	Array=(long *)malloc(Nx*sizeof(long));

	for(i=0;i<Nx;i++)
	{
		Array[i]=0;
	}

	return Array;
}

double **Arraydouble2(long Nx,long Nz)
{
	long i,j;
	double **Array;

	Array=(double **)malloc(Nx*sizeof(double *));

	for(i=0;i<Nx;i++)
	{
		*(Array+i)=(double *)malloc(Nz*sizeof(double));
	}

	for(i=0;i<Nx;i++)
	{
		for(j=0;j<Nz;j++)
		{
			Array[i][j]=0.0;
		}
	}

	return Array;
}

void Freedouble2(long Nx, double **Array)
{
	long i;

	for(i=0;i<Nx;i++)
	{
		free(Array[i]);
	}

	free(Array);
}

void Average_rank(int rank, int size, int Ns, int *t1, int *t2)
{
	int i, es, rs, ss, s1, s2;

	es = Ns / size;
	rs = Ns - es*size;
	ss = 0; s1 = 0;

	for(i = 0; i <= rank; i++)
	{
		ss += es;
		if (i < rs)
		{
			ss++;
		}
		if (i == rank - 1)
		{
			s1 = ss;
		}
	}

	s2 = ss;
	*t1 = s1;
	*t2 = s2;
}

void Matrix_infor(int *Ns, int *Nr, int *numx, int *numz, double *x0, double *x1, double *z0, double *z1, double *dx, double *dz)
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
	sf1 = fscanf(fp1, "%d", &ii1);  sf2 = fgets(cc1, 500, fp1);           // 4
	*Ns = ii1;
	sf1 = fscanf(fp1, "%d", &ii1);  sf2 = fgets(cc1, 500, fp1);           // 5
	*Nr = ii1;

	*dx = (*x1 - *x0) / (*numx - 1);
	*dz = (*z1 - *z0) / (*numz - 1);
}
