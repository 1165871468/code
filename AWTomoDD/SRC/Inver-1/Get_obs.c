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
    int i, j, k, count, num_s, num_r;
    double a, b, dd, rr1, rr2, xx1, xx2, zz1, zz2, ff1, ff2, pi = acos(-1.0);
    double **arrdm1, **arrdm2, **arrdm3, **arrdm4, **arrdm5;
    char flnm1[100], flnm2[100], flnm3[100], flnm4[100];
    FILE *fp1;

	dd = 2000.0;

	sprintf(flnm1, "../Input/arrival_p.txt");

	count = countfile(flnm1);
	arrdm1 = Arraydouble2(count, 3);
	rdfile(flnm1, count, 3, arrdm1);

	for (i = 0; i < count; i++)
	{
		num_s = arrdm1[i][0];
		num_r = arrdm1[i][1];
		sprintf(flnm2, "../Obs/Output/Receivers/source%d/D%d", num_s, num_r);
		sprintf(flnm3, "../Obs/Receivers/source%d/D%d", num_s, num_r);

		fp1 = fopen(flnm3, "w");
		if (fp1 == NULL)
		{
			printf("Can not open this file --%s Get_obs\n", flnm3); exit(1);
		}
		



		printf("%s %.1f %.1f\n", flnm2, arrdm1[i][0], arrdm1[i][1]);

	}


		// fprintf(fp1, "%.1f%10.1f%20.8e\n", arrdm1[i][0], arrdm1[i][1], ff1), fflush(fp1);


////////////////////////////////////////////////////////////////
	freedouble2(count, arrdm1);
	// fclose(fp1);
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
		{
			arrfx[i + 1] = arrfx[i] + 0.000001;
		}

		if (arrfz[i] == arrfz[i + 1])
		{
			arrfz[i + 1] = arrfz[i] + 0.000001;
		}
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
