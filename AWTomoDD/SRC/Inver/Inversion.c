#include "Inversion.h"

int main(int argc, char *argv[])
{
    int i, j, k, l, m, n, nz, nmx, nmz, num_cor, row, col, count1, count2, step, step_max, flgi1, Ns, Nr;
    double xx1, xx2, x0, x1, zz1, zz2, z0, z1, dx, dz, dis, rr;
	double error, error1, ratio, beta, vt, vt1, g1, gama, tmp1, tmp2, tmp3, maxx, eta, eta1, eta2;
	double *velc, *Ta, *Ta1, *Ta2, *Pcg, *Gr, *Gr1, *deltaX, *deltaX1;
    double **arrd_A, *arrd_b, *arrd_x, *data, *Unc, ***arrin;
    char flnm1[100], flnm2[100], flnm3[100], flnm4[100];
	Sarray *Fd, *Ft;
    FILE *fp1, *fp2;

	Matrix_infor(&num_cor, &Ns, &Nr, &row, &col, &nz, &step_max, &eta);

	row = Ns*Nr*(Nr-1)/2;
	nz = row*col;
	arrd_A = Arraydouble2(row, col);     arrd_b = Arraydouble1(row);      arrd_x = Arraydouble1(col);
	deltaX = Arraydouble1(col);          deltaX1 = Arraydouble1(col);
	data = Arraydouble1(row);            Unc = Arraydouble1(row);         velc = Arraydouble1(nz);
	Fd=(Sarray *)malloc(nz*sizeof(Sarray));
	Ft=(Sarray *)malloc(nz*sizeof(Sarray));

	Initial_sarray(nz, Fd);
	Initial_sarray(nz, Ft);
	Initial_velocity(nz, velc);
	Initial_velocity1(nz, velc);

	for (i = 0; i < row; i++) {
		Unc[i] = 1.0;
	}
	for (i = 0; i < col; i++) {
		deltaX[i] = 0.0;
	}

	for(step = 0; step < step_max; step++) {
		printf("Inversion step = %d\n", step+1);

		if (step > 0) {
			Get_gradient(num_cor, Ns, Nr);

			Initial_velocity1(nz, velc);
		}

		Update_matrix(&row, &col, &nz);

		Collect_matrix(step, Ns, Nr, row, col, data, arrd_A);

		Sensitivity_matrix(row, col, nz, arrd_A, Fd);

		Transpose(row, col, nz, Fd, Ft);

		CG(nz, row, col, eta, Fd, Ft, Unc, data, deltaX);

		Model_update(step+1, col, deltaX, velc);

		Output_delta(step+1, nz, row, col, nz, deltaX);
	}

	free(Fd);
	free(Ft);
}

void CG(long nz, long row, long col, double eta, Sarray *Fd, Sarray *Ft, double *Unc, double *data, double *deltaX)
{
	int i, step;
	double error, error1, fixe, ratio, beta, vt, vt1, g1, gama, tmp1, tmp2, tmp3, maxx, eta1, eta2;
	double *Ta, *Pcg, *Gr, *Gr1, *deltaX1;

	Ta = Arraydouble1(row);      Pcg = Arraydouble1(col);
	Gr = Arraydouble1(col);      Gr1 = Arraydouble1(col);
	deltaX1 = Arraydouble1(col);

	for (i = 0; i < col; i++) {
		deltaX[i] = 0.0;
	}

	Error0(row, data, &error);
	fprintf(stdout, "RMS=%e(s)\n", sqrt(error/row)); fflush(stdout);

	if (error < 1.0e-100) {
		printf("The error is too small!!! \n"); exit(-1);
	}

	eta1 = 1.0e-6; eta2 = 0.0; maxx = 10.0;

	while ((maxx > eta) || (maxx < 0.2*eta && eta1 >= 1.0)) {
		step=0; ratio=1.0;

		for (i = 0; i < col; i++) {
			deltaX[i] = 0.0;
		}

		while (step < 20 && ratio > 0.01) {
			Total_error(nz, row, col, eta1, Fd, data, deltaX, Unc, &error);
			Gradient(row, col, nz, eta1, Ta, Fd, Ft, data, Unc, deltaX, Gr);

			// eta1 = 0.0;

			if (step == 0) {
				fixe = error;
				for (i = 0; i < col; i++) {
					Pcg[i] = -1.0*Gr[i];
				}
			}
			else {
				tmp1 = 0.0; tmp2=0.0; tmp3=0.0;

				for (i = 0; i < col; i++) {
					tmp1 += Gr[i]*Gr[i];
					tmp2 += Gr[i]*Gr1[i];
					tmp3 += Gr1[i]*Gr1[i];
				}
			
				beta = (tmp1 - tmp2) / (tmp3 + 1.0e-40);

				if (beta < 0.0) {
					beta = 0.0;
				}

				for (i = 0; i < col; i++) {
					Pcg[i] = -1.0*Gr[i] + beta*Pcg[i];
				}
			}

			g1 = Vector_product(col, Gr, Pcg);

			if(fabs(g1) < 1.0e-30) {
				vt = 0.0;
			}
			else {
				vt = -2.0*error / g1;
			}

			Delta_update(col, vt, Pcg, deltaX, deltaX1);
			Total_error(nz, row, col, eta1, Fd, data, deltaX1, Unc, &error1);

			while (error1 >= error*(1.0 + 1.0e-6)) {
				if (error1 + error < 1.0e-20) {
					vt1 = 0.0;
				}
				else {
					vt1 = error / (error1 + error)*vt;
				}

				Delta_update(col, vt1, Pcg, deltaX, deltaX1);
				Total_error(nz, row, col, eta1, Fd, data, deltaX1, Unc, &error1);
				vt = vt1;
			}

			gama = 2.0*((error1 - error) - g1*vt) / (vt*vt + 1.0e-30);

			if (vt*vt < 1.0e-30) {
				vt = 0.0;
			}
			else if (gama >= 1.0e-10) {
				vt = (-1.0*g1)/gama;
			}

			Delta_update(col, vt, Pcg, deltaX, deltaX);
			Total_error(nz, row, col, eta1, Fd, data, deltaX, Unc, &error);

			ratio = error / ((1.0 + 1.0e-3)*fixe);
			step++;

			for (i = 0; i < col; i++) {
				Gr1[i] = Gr[i];
			}
		}
		maxx = Max_abs_vec(col, deltaX);

		if (maxx > eta) {
			eta2 = eta1;
			eta1 = 2.0*eta1;
		}
		else if (maxx < 0.1*eta) {	
			eta1 = 0.5*(eta1 + eta2);
		}
		else {
			break;
		}
	}
	printf("before: %.4e\nafter: %.4e\n", fixe, error);
	fprintf(stdout, "damping=%e; maximal perturbation=%e\n\n", eta1, maxx); fflush(stdout);
	
	free(Ta);   free(Pcg);
	free(Gr);   free(Gr1);
	free(deltaX1);
}

void Model_update(int step, int col, double *deltaX, double *velc)
{
	int i, j, k, ix, jz, n, id1, count1, count2, flgi1, flgi2, nmx, nmz, fs;
	char flnm1[100], flnm2[100], flnm3[100], flnm4[100], flnm5[100], flnm6[100], flnm7[100];
	double x0, z0, x1, z1, xx1, zz1, xx2, zz2, xx3, zz3, dx, dz, ff, ff1;
	double *arrdx, *arrdz, **arrdy, **arrd_for1, **arrd_for2, **arrd_in1, **arrd_in2, ***arrd_in3;
	FILE *fp1, *fp2, *fp3, *fp4;

	modInforIn(&nmx, &nmz, &x0, &x1, &z0, &z1);

	dx = (x1 - x0)/(nmx - 1);
	dz = (z1 - z0)/(nmz - 1);

	sprintf(flnm1, "../Output/Result/vel%d.dat", step);
	sprintf(flnm2, "../Output/Model/Forver/vel_initial.dat");
	sprintf(flnm3, "../Output/Model/Forver/vel.dat");
	sprintf(flnm4, "../Output/Model/Inver/vel.dat");
	sprintf(flnm5, "../Output/Model/Inver/vel_initial.dat");
	sprintf(flnm6, "../Output/Result/velForwardPer%d.dat.dat", step);

	fp1 = fopen(flnm1, "w");
	fp2 = fopen(flnm3, "w");
	fp3 = fopen(flnm4, "w");
	fp4 = fopen(flnm6, "w");

	if (fp1 == NULL || fp2 == NULL || fp3 == NULL) {
		printf("Can't open this file-- %s !\n", flnm1), exit(1);
	}

	count1 = col;
	count2 = Count_file(flnm2);

	arrd_for1 = Arraydouble2(count2, 5);
	arrd_for2 = Arraydouble2(count2, 5);
	arrd_in1 = Arraydouble2(count1, 1);
	arrd_in2 = Arraydouble2(count1, 5);
	arrd_in3 = Arraydouble3(nmx, nmz, 3);
	///
	//////
	for(i = 0; i < count1; i++) {
		velc[i] = (1.0 - deltaX[i])*velc[i];
		arrd_in1[i][0] = velc[i];

		fprintf(fp1, "%.16e\n", velc[i]);
	}
	///
	Rdfile(flnm5, count1, 5, arrd_in2);

	for(i = 0; i < count1; i++) {
		fprintf(fp3, "%.4f %.4f %.4f %.4f %.4f\n", arrd_in2[i][0], arrd_in2[i][1], velc[i], arrd_in2[i][3], arrd_in2[i][4]);
	}
	///
	//////
	id1 = 0;
	for (i = 0; i < nmx; i++) {
		for (j = 0; j < nmz; j++) {
			arrd_in3[i][j][0] = x0 + i*dx;
			arrd_in3[i][j][1] = z0 + j*dz;
			arrd_in3[i][j][2] = arrd_in1[id1][0];

			id1++;
		}
	}

	Rdfile(flnm2, count2, 5, arrd_for1);
	///
	//////
	for (i = 0; i < count2; i++) {
		xx1 = arrd_for1[i][0];
		zz1 = arrd_for1[i][1];

		flgi1 = 0; 
		for (j = 0; j < nmx; j++) {
			for (k = 0; k < nmz; k++) {
				xx2 = arrd_in3[j][k][0];
				zz2 = arrd_in3[j][k][1];

				xx3 = xx2 + dx;
				zz3 = zz2 + dz;

				if (xx1 >= xx2 && xx1 < xx3 && zz1 >= zz2 && zz1 < zz3 && j < nmx - 1 && k < nmz - 1) {
					ix = j;  jz = k;
					flgi1 = 1;
					break;
				}
			}

			if (flgi1 == 1) {
				break;
			}
		}

		if (flgi1 == 1) {
			ff = interpolate2D(ix, jz, xx1, zz1, arrd_in3);
		}
		else {
			ff = arrd_for1[i][2];
		}

		ff1 = (ff - arrd_for1[i][2])/arrd_for1[i][2];
		//////

		fprintf(fp2, "%.4f %.4f %.4f %.4f %.4f\n", arrd_for1[i][0], arrd_for1[i][1], ff, arrd_for1[i][3], arrd_for1[i][4]);
		fprintf(fp4, "%.4f %.4f %.4f %.4f %.4f\n", arrd_for1[i][0], arrd_for1[i][1], ff1, 0.0, 0.0);
	}

	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);

	Freedouble2(count1, arrd_in1);
	Freedouble2(count1, arrd_in2);
	Freedouble2(count2, arrd_for1);
	Freedouble2(count2, arrd_for2);
	Freedouble3(nmx, nmz, arrd_in3);
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

double Max_abs_vec(long len, double *vector)
{
	long i;
	double tmp;

	tmp = 0.0;
	for (i = 0; i < len; i++)
	{
		if (tmp < fabs(vector[i]))
		{
			tmp = fabs(vector[i]);
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

double interpolate2D(int ix, int jy, double x, double y, double ***arrayData)
{
    int i, j;
    double cof1, cof2, ff;

    ff = 0.0;
	for (i = 0; i < 2; i++) {
		for (j = 0; j < 2; j++) {
			cof1 = (x - arrayData[ix+i][jy+j][0]) / (arrayData[ix+1][jy+1][0] - arrayData[ix][jy][0]);
			cof2 = (y - arrayData[ix+i][jy+j][1]) / (arrayData[ix+1][jy+1][1] - arrayData[ix][jy][1]);
			cof1 = 1.0 - fabs(cof1);
			cof2 = 1.0 - fabs(cof2);

			ff += arrayData[ix+i][jy+j][2]*cof1*cof2;
		}
	}

    return ff;
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

FILE *Open_file(const char *flnm, const char *mode)
{
    FILE *fp1;

    if ((fp1 = fopen(flnm, mode)) == NULL)
    {
        perror(flnm); exit(1);
    }
    return fp1;
}

void Rdfile(char *flnm, int row, int col, double **array_data)
{
    int i, j, sf;
    double tt[100];
    char ch, FF[100], *ff;
    FILE *fp;

    fp = Open_file(flnm, "r");
    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            sf = fscanf(fp, "%lf", &tt[j]);
            array_data[i][j] = tt[j];
        }
    }
    fclose(fp);
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

void Output_delta(int step, int nz, int row, int col, int Nz, double *deltaX)
{
	int i, j, fs;
	char flnm1[100];
	double xx, yy, zz, ff;
	FILE *fp1;

	sprintf(flnm1,"../Output/Result/vel_per_in%d.dat", step);

	if ((fp1 = fopen(flnm1, "w")) == NULL) {
		printf("Cannot open V.dat %s\n", flnm1); exit(1);
	}
	for (i = 0; i < col; i++) {
		fprintf(fp1, "%.16e\n", deltaX[i]);
	}
	fclose(fp1);
}

void Collect_matrix(int step, int Ns, int Nr, long row, long col, double *data, double **A)
{
	int i, j, ii1, ii2, ii3, id1, id2, fs, count1, count2, count3;
	double tt1, ff1, **arrd1;
	char flnm1[100], flnm2[100], cf[250];
	FILE *fp1;

	sprintf(flnm1, "../Output/Model/Inver/Sensitivity_in_dd.dat");
	count1 = Count_file(flnm1);
	arrd1 = Arraydouble2(count1, 1);
	Rdfile(flnm1, count1, 1, arrd1);

	sprintf(flnm2, "../Output/Error/Error_dd.txt");
	count2 = Count_file(flnm2);

	if (count2 != row) {
		printf("The data is not marching! -- Collect_matrix\n"); exit(1);
	}
	else {
		sprintf(flnm2, "../Output/Error/Error_dd.txt");

		fp1 = fopen(flnm2, "r");
		if (fp1 == NULL) {
			printf("Can not open tisn file -- %s\n", flnm2); exit(1);
		}

		for (i = 0; i < row; i++) {
			fs = fscanf(fp1, "%d %d %d %lf", &ii1, &ii2, &ii3, &tt1);
			data[i] = -1.0*tt1;
		}
		fclose(fp1);
	}

	id2 = 0;
	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			A[i][j] = arrd1[id2][0];
			id2++;
		}
	}

	Freedouble2(count1, arrd1);
	sprintf(flnm1, "../Output/Error/Error.txt");
	sprintf(flnm2, "../Output/Error/Error%d.txt", step);
	sprintf(cf, "cp %s %s", flnm1, flnm2);
	fs = system(cf);
}


void Matrix_infor(int *num_cor, int *Ns, int *Nr, int *row, int *col, int *nz, int *step, double *eta)
{
	int i, j, ii1, ii2, sf1;
	double dd1, dd2;
	char cc1[500], flnm1[100], flnm2[100], *sf2;
	FILE *fp1;

	sprintf(flnm1, "../Input/parameters_inversion.config");

	if ((fp1 = fopen(flnm1, "r")) == NULL)
	{
		printf("Can not read this file -- %s", flnm1); exit(1);
	}
	////////////////////////////////////////////////////////////////////////////
	sf1 = fscanf(fp1, "%lf%lf", &dd1, &dd2);  sf2 = fgets(cc1, 500, fp1); // 1

	sf1 = fscanf(fp1, "%lf%lf", &dd1, &dd2);  sf2 = fgets(cc1, 500, fp1); // 2

	sf1 = fscanf(fp1, "%d%d", &ii1, &ii2);  sf2 = fgets(cc1, 500, fp1);   // 3
	*col = ii1*ii2;
	sf1 = fscanf(fp1, "%d", &ii1);  sf2 = fgets(cc1, 500, fp1);           // 4
	*Ns = ii1;
	sf1 = fscanf(fp1, "%d", &ii1);  sf2 = fgets(cc1, 500, fp1);           // 5
	*Nr = ii1;
	sf1 = fscanf(fp1, "%d", &ii1);  sf2 = fgets(cc1, 500, fp1);           // 6
	*step = ii1;
	sf1 = fscanf(fp1, "%lf", &dd1);  sf2 = fgets(cc1, 500, fp1);          // 7
	*eta = dd1;
	sf1 = fscanf(fp1, "%lf", &dd1);  sf2 = fgets(cc1, 500, fp1);          // 8

	sf1 = fscanf(fp1, "%lf", &dd1);  sf2 = fgets(cc1, 500, fp1);          // 9

	sf1 = fscanf(fp1, "%d", &ii1);  sf2 = fgets(cc1, 500, fp1);          // 10
	*num_cor = ii1;
	////////////////////////////////////////////////////////////////////////////
	*row = Count_file("../Output/Error/Error_dd.txt");

	dd1 = *row; dd2 = *col;
	*nz = dd1*dd2;
}

void Update_matrix(int *row, int *col, int *nz)
{
	double d1, d2;
	d1 = *row;
	d2 = *col;

	d1 = Count_file("../Output/Error/Error_dd.txt");
	d2 = *col;

	*row = d1;
	*col = d2;
	*nz = d1*d2;

}

void Get_gradient(int num, int Ns, int Nr)
{
	int i, fs;
	char ff[100];
 
	printf("Calculating forward wavefiled ...\n");
	sprintf(ff, "mpirun -np %d ./Update", num);
	fs = system(ff);

	printf("Sorting error ...\n");
	sprintf(ff, "./timeDiff");
	fs = system(ff);
	sprintf(ff, "./Error");
	fs = system(ff);

	printf("Calculating sensitivity matrix ...\n");
	sprintf(ff, "mpirun -np %d ./Sensitivity", num);
	fs = system(ff);
}

void Initial_velocity(int nz, double *velc)
{
	int i, fs;
	char flnm1[100];
	double tt1, tt2, tt3, tt4, tt5;
	FILE *fp1;

	sprintf(flnm1, "../Output/Model/Inver/vel_initial.dat");
	
	if ((fp1 = fopen(flnm1, "r")) == NULL) {
		printf("Can not open this file -- %s -- Initial_velocity\n", flnm1); exit(1);
	}
	
	for (i = 0; i < nz; i++) {
		fs = fscanf(fp1, "%lf%lf%lf%lf%lf", &tt1, &tt2, &tt3, &tt4, &tt5);
		velc[i] = tt3;
	}
	
	fclose(fp1);
}

void Initial_velocity1(int nz, double *velc)
{
	int i, fs;
	char flnm1[100];
	double tt1, tt2, tt3, tt4, tt5;
	FILE *fp1;

	sprintf(flnm1, "../Output/Model/Inver/vel.dat");
	
	if ((fp1 = fopen(flnm1, "r")) == NULL) {
		printf("Can not open this file -- %s -- Initial_velocity1\n", flnm1); exit(1);
	}
	
	for (i = 0; i < nz; i++) {
		fs = fscanf(fp1, "%lf%lf%lf%lf%lf", &tt1, &tt2, &tt3, &tt4, &tt5);
		velc[i] = tt3;
	}
	
	fclose(fp1);
}

void modInforIn(int *numx, int *numz, double *x0, double *x1, double *z0, double *z1)
{
	int i, j, ii1, ii2, sf1;
	double dd1, dd2;
	char cc1[500], flnm1[100], flnm2[100], *sf2;
	FILE *fp1;

	sprintf(flnm1, "../Input/parameters_inversion.config");

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
