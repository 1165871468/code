#include "Inversion.h"

int main(int argc, char *argv[])
{
    int i, j, k, l, m, n, nz, nmx, nmz, row, col, count1, count2, step, flgi1, Ns, Nr;

    char flnm1[100], flnm2[100], flnm3[100], flnm4[100];

    double xx1, xx2, x0, x1, zz1, zz2, z0, z1, dx, dz, dis, rr;
	double error, error1, ratio, beta, vt, vt1, g1, gama, tmp1, tmp2, tmp3, maxx, eta, eta1, eta2;
	double *Ta, *Ta1, *Ta2, *Pcg, *Gr, *Gr1, *deltaX, *deltaX1;
    double **arrd_a, *arrd_b, *arrd_x, *data, *Unc, ***arrin;
	Sarray *Fd, *Ft;
    FILE *fp1, *fp2;

	Ns = 3;   Nr = 21;
	m = 3;    n = 3;
	row = m;  col = n;
	nz = m*n; eta = 0.1;

	arrd_a = Arraydouble2(m, n); arrd_b = Arraydouble1(n);  arrd_x = Arraydouble1(n);
	deltaX = Arraydouble1(m);    deltaX1 = Arraydouble1(m);
	data = Arraydouble1(n);      Unc = Arraydouble1(n);

	arrd_a[0][0] = 4.0; arrd_a[0][1] = 3.0;  arrd_a[0][2] = 0.0; /// 111
	arrd_a[1][0] = 3.0; arrd_a[1][1] = 4.0;  arrd_a[1][2] = -1.0; /// 111
	arrd_a[2][0] = 0.0; arrd_a[2][1] = -1.0; arrd_a[2][2] = 4.0; /// 111

	arrd_b[0] = 0.3;    arrd_b[1] = 0.5;    arrd_b[2] = -0.5;
	arrd_x[0] = 0.0;    arrd_x[1] = 0.0;     arrd_x[2] = 0.0;

	Fd=(Sarray *)malloc(nz*sizeof(Sarray));
	Ft=(Sarray *)malloc(nz*sizeof(Sarray));
	Initial_sarray(nz, Fd);
	Initial_sarray(nz, Ft);

	for (i = 0; i < row; i++)
	{
		data[i] = arrd_b[i];
		Unc[i] = 1.0;
	}
	
	for (i = 0; i < col; i++)
	{
		deltaX[i] = 0.0;
	}
	
	for(step = 0; step < 20; step++)
	{
		fprintf(stdout,"Inversion step = %d\n",step + 1); fflush(stdout);

		Sensitivity_matrix(row, col, nz, arrd_a, Fd);
		Transpose(row, col, nz, Fd, Ft);

		CG(nz, row, col, eta, Fd, Ft, Unc, data, deltaX);
	}

	for (i = 0; i < col; i++)
	{
		printf("%10.4f\n", deltaX[i]);
	}

	free(Fd);
}

void CG(long nz, long row, long col, double eta, Sarray *Fd, Sarray *Ft, double *Unc, double *data, double *deltaX)
{
	int i,j, k, step;
	double error, error1, fixe, ratio, beta, vt, vt1, g1, gama, maxx, eta1, eta2, tmp1, tmp2, tmp3;
	double *Ta,*Pcg,*Gr,*Gr1,*deltaX0,*deltaX1;

	Ta = Arraydouble1(row);      Pcg = Arraydouble1(col);
	Gr = Arraydouble1(col);      Gr1 = Arraydouble1(col);
	deltaX0 = Arraydouble1(col);
	deltaX1 = Arraydouble1(col);
	for(i = 0; i < col; i++)
	{
		deltaX0[i] = deltaX[i];
	}

	Error0(row, data, &error);

	fprintf(stdout, "RMS=%e(s)\n", sqrt(error / row)); fflush(stdout);

	eta1 = 10.0; eta2 = 0.0; maxx = 1.0;

	step = 0;
//	while((maxx > eta) || (maxx < 0.2*eta && eta1 >= 1.0))
	while(step < 20)
	{
		ratio = 1.0;
		while(step < 20 && ratio > 0.5)
		{
			Total_error(nz, row, col, eta1, Fd, data, deltaX, Unc, &error);
			Gradient(row, col, nz, eta1, Ta, Fd, Ft, data, Unc, deltaX, Gr);
			eta1 = 0.0;

			if(step == 0)
			{
				fixe = error;

				for(i = 0; i < col; i++)
				{
					Pcg[i] = -1.0*Gr[i];
				}
			}
			else
			{
				tmp1 = 0.0; tmp2 = 0.0; tmp3 = 0.0;
				for(i = 0; i < col; i++)
				{
					tmp1 += Gr[i]*Gr[i];
					tmp2 += Gr[i]*Gr1[i];
					tmp3 += Gr1[i]*Gr1[i];
				}
				beta = (tmp1 - tmp2) / (tmp3 + 1.0e-20);

				if(beta < 0.0)
				{
					beta = 0.0;
				}

				for(i = 0; i < col; i++)
				{
					Pcg[i] = -1.0*Gr[i] + beta*Pcg[i];
				}
			}

			g1 = Vector_product(col, Gr, Pcg);

			if(fabs(g1) < 1.0e-32)
			{
				vt = 0.0;
			}
			else
			{
				vt = -2.0*error / g1;
			}

			Delta_update(col, vt, Pcg, deltaX, deltaX1);
			Total_error(nz, row, col, eta1, Fd, data, deltaX1, Unc, &error1);

			while (error1 >= (1.0 + 1.0e-6)*error)
			{
				if(error1 + error < 1.0e-10)
				{
					vt1 = 0.0;
				}
				else
				{
					vt1 = error / (error1 + error)*vt;
				}

				Delta_update(col, vt1, Pcg, deltaX, deltaX1);
				Total_error(nz, row, col, eta1, Fd, data, deltaX1, Unc, &error1);
				vt = vt1;
			}
			// printf("******************************************\n");

			gama = ((error1 - error) - g1*vt) / (vt*vt + 1.0e-25);

			if (vt*vt < 1.0e-30)
			{
				vt = 0.0;
			}
			else if (gama >= 1.0e-10)
			{
				vt = -1.0*g1 / (2.0*gama);
			}

			Delta_update(col, vt, Pcg, deltaX, deltaX);
			Total_error(nz, row, col, eta1, Fd, data, deltaX, Unc, &error);
			ratio = error / (fixe + 1.0e-8);

			for(i = 0; i < col; i++)
			{
				Gr1[i] = Gr[i];
			}

			step++;
		}

		maxx = Max_abs_vec(col, deltaX0, deltaX);

		if(maxx > eta)
		{
			eta2 = eta1;
			eta1 = 2.0*eta1;
			printf("111\n");
		}
		else if(maxx < 0.1*eta)
		{	
			eta1 = (eta1 + eta2) / 2.0;
			printf("222\n");
		}
		else
		{
			printf("333\n");
			break;
		}

		eta1 = 0.0;
	}

	fprintf(stdout,"damping = %e; maximal perturbation = %e; RMS = %e(s)\n", eta1, maxx, sqrt(error/row)); fflush(stdout);
	
	free(Ta);      free(Pcg);
	free(Gr);      free(Gr1);
	free(deltaX1);
}

void Model_update(int row, int col, int nz, double *deltaX, double *data, Sarray *Fd)
{
	int i;
	double *arrd1;

	arrd1 = Arraydouble1(row);

// 	for(i = 0; i < col; i++)
// 	{
// //		arrd_x[i] += deltaX[i];
// 		deltaX[i] += deltaX[i];
// 	}

	for(i = 0; i < nz; i++)
	{
		arrd1[Fd[i].row] += Fd[i].val*deltaX[Fd[i].col];
	}

	for (i = 0; i < row; i++)
	{
		data[i] = data[i] - arrd1[i];
	}
	
	free(arrd1);
}

void Delta_update(long col, double vt, double *Pcg, double *deltaX, double *deltaX1)
{
	long i;
	for(i = 0; i < col; i++)
	{
		deltaX1[i] = deltaX[i] + vt*Pcg[i];
	}	
}

double Vector_product(long len, double *V1, double *V2)
{
	long i;
	double tmp;

	tmp=0.0;
	for(i=0;i<len;i++)
	{
		tmp+=V1[i]*V2[i];
	}

	return tmp;
}

double Max_abs_vec(long len, double *vector0, double *vector)
{
	long i;
	double tmp;

	tmp = 0.0;
	for (i = 0; i < len; i++)
	{
		if (tmp < fabs(vector[i] - vector0[i]))
		{
			tmp = fabs(vector[i] - vector0[i]);
		}
	}
	return tmp;
}

void Error0(long row, double *data, double *error)
{
	long i;
	double err;

	err = 0.0;
	for(i = 0; i < row; i++)
	{
		err += data[i]*data[i];
	}
	*error = err / 2.0;
}

void Gradient(long row, long col, long nz, double eta, double *Ta, Sarray *Fd, Sarray *Ft, double *data, double *Unc, double *deltaX, double *Gr)
{
	long i;

	for(i = 0; i < row; i++)
	{
		Ta[i] = 0.0;
	}
	for(i = 0; i < col; i++)
	{
		Gr[i] = 0.0;
	}

	for(i = 0; i < nz; i++)
	{
		Ta[Fd[i].row] += Fd[i].val*deltaX[Fd[i].col];
		// printf("Ax: %10.4f\n", Ta[Fd[i].row]);
	}

	for(i = 0; i < row; i++)
	{
		Ta[i] -= data[i];
		Ta[i] = Ta[i] / Unc[i];
		// printf("Ta: %10.4f\n", Ta[i]);
	}

	for(i = 0; i < nz; i++)
	{
		Gr[Ft[i].row] += Ft[i].val*Ta[Ft[i].col];
		// printf("%10.4f\n", Gr[Ft[i].row]);
	}

	for(i = 0; i < col; i++)
	{
//		Gr[i] = 2.0*Gr[i] + eta*deltaX[i];
		Gr[i] = Gr[i] + eta*deltaX[i];
		// printf("Gr: %10.4f\n", Gr[i]);
	}
}

void Total_error(long nz, long row, long col, double eta, Sarray *Fd, double *data, double *deltaX, double *Unc, double *error)
{
	long i;
	double err, err1, *tmp;

	tmp = Arraydouble1(row);

	for(i = 0; i < nz; i++)
	{
		tmp[Fd[i].row] += Fd[i].val*deltaX[Fd[i].col];
	}

	err = 0.0; err1 = 0.0;

	for(i = 0; i < row; i++)
	{
		err += (tmp[i] - data[i])*(tmp[i] - data[i]) / Unc[i];
	}

	for(i = 0; i < col; i++)
	{
//		err1 += 0.5*eta*deltaX[i]*deltaX[i];
		err1 += eta*deltaX[i]*deltaX[i];
//		printf("%10.4f %10.4f\n", eta, deltaX[i]);
	}
	*error = (err + err1) / 2.0;
	// *error = err / 2.0;

//	printf("error+++++++++++++++++: %20.4e%20.4e%20.4e\n", err, err1, *error);

	free(tmp);
}

void Sensitivity_matrix(long row, long col, long nz, double **array_a, Sarray *Fd)
{
	long i, j, i1, i2, count1;

	count1 = 0;
	for(i = 0; i < row; i++)
	{	i1 = i;
		for (j = 0; j < col; j++)
		{
			i2 = j;
			Fd[count1].row = i1;
			Fd[count1].col = i2;
			Fd[count1].val = array_a[i][j];
			count1++;
		}
	}
}

void Transpose(long row, long col, long nz, Sarray *Fd, Sarray *Ft)
{
	long i, trans, *tmp1, *tmp2;

	tmp1=Arraylong1(col);
	tmp2=Arraylong1(col);
	
	for(i = 0; i < nz; i++)
	{
		tmp1[Fd[i].col]++;
	}

	for(i = 1; i < col; i++)
	{
		tmp2[i] = tmp2[i-1] + tmp1[i-1];
	}

	for(i = 0; i < nz; i++)
	{
		trans = tmp2[Fd[i].col]++;
		Ft[trans].col = Fd[i].row;
		Ft[trans].row = Fd[i].col;
		Ft[trans].val = Fd[i].val;
	}

	free(tmp1); free(tmp2);
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

int countfile(char *FL)
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
	int i, ff1;
	char FF[100],*ff;
	float tmp1,tmp2,tmp3,tmp4,tmp5;
	FILE *fp;

	if ((fp = fopen(FL,"r")) == NULL)
	{
		printf("can't open file --- %s\n", FL), exit(1);
	}

	if (n == 1)
	{
		for ( i = 0; i < count; i++)
		{
			ff1 = fscanf(fp,"%f",&tmp1);
			ary[i][0]=tmp1;
		}
	}

	if (n == 2)
	{
		for ( i = 0; i < count; i++)
		{
			ff1 = fscanf(fp,"%f%f",&tmp1,&tmp2);
			ary[i][0]=tmp1;
			ary[i][1]=tmp2;
		}
	}
	
	if (n == 3)
	{
		for ( i = 0; i < count; i++)
		{
			ff1 = fscanf(fp,"%f%f%f",&tmp1,&tmp2,&tmp3);
			ary[i][0]=tmp1;
			ary[i][1]=tmp2;
			ary[i][2]=tmp3;
		}
	}

	if (n == 4)
	{
		for ( i = 0; i < count; i++)
		{
			ff1 = fscanf(fp,"%f%f%f%f",&tmp1,&tmp2,&tmp3,&tmp4);
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
			ff1 = fscanf(fp,"%f%f%f%f%f",&tmp1,&tmp2,&tmp3,&tmp4,&tmp5);
			ary[i][0]=tmp1;
			ary[i][1]=tmp2;
			ary[i][2]=tmp3;
			ary[i][3]=tmp4;
			ary[i][4]=tmp5;
		}
	}

	rewind(fp), fclose(fp);
}

void Initial_sarray(long nz, Sarray *Fd)
{
	long i;
	for(i = 0; i < nz; i++)
	{
		Fd[i].row = 0;
		Fd[i].col = 0;
		Fd[i].val = 0.0;
	}
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
