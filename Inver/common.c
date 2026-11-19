#include "common.h"

void Matrix_multi(double **A, int m1, int m2, double **B, int n2, double **C)
{
    int i, j, k;
    double flgd1;

    for (i = 0; i < n2; i++)
    {
        for (j = 0; j < m1; j++)
        {
            flgd1 = 0;

            for (k = 0; k < m2; k++)
            {
                flgd1 += A[j][k] * B[k][i];
            }
            
            k--;
            C[j][i] = flgd1;
        }
    }
}

double Vector_product(long len, double *V1, double *V2)
{
	long i;
	double tmp;

	tmp = 0.0;
	for(i = 0; i < len; i++)
	{
		tmp += V1[i]*V2[i];
	}

	return tmp;
}

double Interpolate(int n, double x, double y, double *arrfx, double *arrfz, double **arrfy)
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
	char FF[200],*ff;
	FILE *fp;

	fp = fopen(FL,"r");

	if (fp == NULL)
	{
		printf("can't open file --- %s\n", FL), exit(1);
	}


	count = 0;

	while((ff = fgets(FF, sizeof FF, fp))!=NULL)
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
	{
		printf("can't open file --- %s\n", FL), exit(1);
	}

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
	{
		free(Array[i]);
	}

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