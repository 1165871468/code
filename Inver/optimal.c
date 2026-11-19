#include "common.h"
#include "optimal.h"

void Output_model(int step, int count1, double *deltaX2)
{
	int i;
	char flnm1[100], flnm2[100];
	double xx, zz, per, **arrayd1;
	FILE *fp1;

	sprintf(flnm1, "../Output/Model/Inver/Kernel_lam_in.dat");
	sprintf(flnm2, "../Output/Model/Result/Kernel_lam_in%d.dat", step);

	arrayd1 = Arraydouble2(count1, 3);
	Rdfile(flnm1, count1, 3, arrayd1);

	fp1 = fopen(flnm2, "w");
	if (fp1 == NULL) {
		printf("can't open this file! -- Output_model\n"), exit(1);
	}

	for (i = 0; i < count1; i++)
	{
		xx = arrayd1[i][0];
		zz = arrayd1[i][1];
		per = deltaX2[i];
		fprintf(fp1, "%.1f%10.1f%20.8e\n", xx, zz, per), fflush(fp1);
	}

	Freedouble2(count1, arrayd1);
	fclose(fp1);
}

void Gradient(int Ns, int Nr, int count1, double *Gr)
{
	int i, ss;
	double tmp1, **arrayd1;
	char ff[100], flnm1[100];

	printf("Calcualting adjoint source ...\n");

	sprintf(ff, "./Gradient");
	ss = system(ff);
	sprintf(ff, "./Invergrid_lam");
	ss = system(ff);
	sprintf(ff, "./Smooth");
	ss = system(ff);
}

void Get_gradient(int count1, double *Gr)
{
	int i, count2, count3;
	char flnm1[100], flnm2[100];
	double dd, k, tmp1;
	double **arrayd1, **arrayd2;

	// k = 1.40919992e+10;
	dd = 1.0e6;
	// dd = 1.0;
	sprintf(flnm1, "../Output/Model/Inver/Kernel_lam_in_smooth.dat");
	sprintf(flnm2, "../Output/Model/Inver/lam_in.dat");

	count2 = Count_file(flnm1);
	count3 = Count_file(flnm2);

	if ((count1 != count2) || (count1 != count3))
	{
		printf("The data is not marching! -- Get_gradient\n"); exit(1);
	}

	arrayd1 = Arraydouble2(count1, 3);
	arrayd2 = Arraydouble2(count1, 3);

	Rdfile(flnm1, count1, 3, arrayd1);
	Rdfile(flnm2, count1, 3, arrayd2);

	// for (i = 0; i < count1; i++)
	// {
	// 	arrayd1[i][2] = 1.0;
	// 	arrayd2[i][2] = 1.0;
	// }

	for (i = 0; i < count1; i++)
	{
		Gr[i] = dd*arrayd1[i][2]*arrayd2[i][2];
	}

	Freedouble2(count1, arrayd1);
	Freedouble2(count1, arrayd2);
}

void Model_update(int Ns, int Nr, double *error)
{
    int ss;
	char ff[100];
	double tmp1;

	printf("Forword modeling ...\n");

	sprintf(ff, "./Update");
	ss = system(ff);
	tmp1 = Total_error(Ns, Nr);
	*error = tmp1;
}

void Delta_update(int col,double vt, double *Pcg, double *deltaX1, double *deltaX2)
{
    int i, j, k, l, n, fs, count1, count2, nmx, nmz, id1, ix, jz, flgi1;
    char flnm1[100], flnm2[100], flnm3[100], flnm4[100], flnm5[100], flnm6[100];
    double ff, xx, xx1, xx2, x0, x1, zz, zz1, zz2, z0, z1, dx, dz;
    double *arrdx, *arrdz, **arrdy, **arrin, **arrin1, ***arrin2, **arrin3;
	FILE *fp1, *fp2, *fp3;

	sprintf(flnm1, "../Output/Model/Forver/lam_initial.txt");
	sprintf(flnm2, "../Output/Model/Forver/per_lam_for.dat");
	sprintf(flnm3, "../Output/Model/Forver/lam_syn.txt");
	sprintf(flnm4, "../Output/Model/Inver/lam_in_intial.dat");
	sprintf(flnm5, "../Output/Model/Inver/lam_in.dat");
	sprintf(flnm6, "./test.dat");

	count1 = col;

	fp3 = fopen(flnm6, "w");

	for(i = 0; i < count1; i++)
	{
		deltaX2[i] = deltaX1[i] + vt*Pcg[i];

//		printf("111 %20.8e\n", deltaX2[i]);

		fprintf(fp3, "%20.8e\n", deltaX2[i]), fflush(fp3);
	}

	fclose(fp3);
// 	printf("ccc4 %10d\n", count1), sleep(1000);

	fp1 = fopen(flnm2, "w"); fp2 = fopen(flnm3, "w"); fp3 = fopen(flnm5, "w");
	if (fp1 == NULL || fp2 == NULL || fp3 == NULL)
	{
		printf("can't open this file! -- Delta_update\n"), exit(1);
	}

	x0 = 0.0, x1 = 100000.0, dx = 1000.0;
	z0 = 0.0, z1 = 40000.0,  dz = 1000.0;

	nmx = (x1 - x0) / dx + 1;
	nmz = (z1 - z0) / dz + 1;

	arrin2 = Arraydouble3(nmx, nmz, 3);
	arrin3 = Arraydouble2(count1, 3);
	Rdfile(flnm4, count1, 3, arrin3);

	count1 = col;
	count2 = nmx*nmz;
	if (count1 != count2)
	{
		printf("Invergrid is error! -- Delta_update\n"), exit(1);
	}
	
	count2 = Count_file(flnm1);
	arrin = Arraydouble2(count2, 5);
	arrin1 = Arraydouble2(count2, 1);
	Rdfile(flnm1, count2, 5, arrin);

	for(i = 0; i < count1; i++)
	{
		deltaX2[i] = deltaX1[i] + vt*Pcg[i];
	}

	id1 = 0;
	for (i = 0; i < nmx; i++)
	{
		for (j = 0; j < nmz; j++)
		{
			arrin2[i][j][0] = x0 + i * dx;
			arrin2[i][j][1] = z0 + j * dz;
			arrin2[i][j][2] = deltaX2[id1];
			ff = arrin3[id1][2]*(1.0 + deltaX2[id1]);

			fprintf(fp3, "%10.1f%10.1f%20.8e\n", arrin2[i][j][0], arrin2[i][j][1], ff), fflush(fp3);
			id1++;
		}
	}

	n = 2;
	arrdx = Arraydouble1(n);
	arrdz = Arraydouble1(n);
	arrdy = Arraydouble2(n, n);

	for (i = 0; i < count2; i++)
	{
		xx = arrin[i][0];
		zz = arrin[i][1];

		flgi1 = 0; 
		for (j = 0; j < nmx; j++)
		{
			for (k = 0; k < nmz; k++)
			{
				xx1 = arrin2[j][k][0];
				zz1 = arrin2[j][k][1];

				xx2 = xx1 + dx;
				zz2 = zz1 + dz;

				if (xx >= xx1 && xx < xx2 && zz >= zz1 && zz < zz2 && j < nmx - 1 && k < nmz - 1)
				{
					ix = j;  jz = k;
					flgi1 = 1;
					break;
				}
				
			}
			if (flgi1 == 1)
			{
				break;
			}
		}

		if (flgi1 == 1)
		{
			for (j = 0; j < n; j++)
			{
				for (k = 0; k < n; k++)
				{
					arrdx[j] = arrin2[ix+j][jz+k][0];
					arrdz[k] = arrin2[ix+j][jz+k][1];
					arrdy[j][k] =  arrin2[ix+j][jz+k][2];
				}
			}

			ff = Interpolate(n, xx, zz, arrdx, arrdz, arrdy);
		}
		else
		{
			ff = 0.0;
		}

		arrin1[i][0] = ff;
		fprintf(fp1, "%.1f%10.1f%20.8e%20.8e%20.8e\n", xx, zz, ff, 0.0, 0.0), fflush(fp1);
	}

	for (i = 0; i < count2; i++)
	{
		ff = arrin[i][2]*(1.0 + arrin1[i][0]);
		// printf("%10d%10d%20.4e\n", count2, i, ff); sleep(1);
		fprintf(fp2, "%.1f%10.1f%20.8e%20.8e%20.8e\n", arrin[i][0], arrin[i][1], ff, arrin[i][3], arrin[i][4]), fflush(fp2);
	}

	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	free(arrdx);
	free(arrdz);
	Freedouble2(n, arrdy);
	Freedouble2(count2, arrin);
	Freedouble2(count2, arrin1);
	Freedouble3(nmx, nmz, arrin2);
	Freedouble2(count1, arrin3);
}

double Total_error(int Ns, int Nr)
{
	int i, j, sf, count1;
	char flnm1[100], flnm2[100], flnm3[100];
	double ff, error;
	FILE *fp1;

	error = 0.0;
	for (i = 0; i < Ns; i++)
	{
		sprintf(flnm1, "../Output/Error/source%d/error.txt", i+1);

		fp1 = fopen(flnm1, "r");
		if (fp1 == NULL)
		{
			printf("Cannot read this file-- %s!\n", flnm1), exit(1);
		}

		sf = fscanf(fp1,"%lf",&ff);
		fclose(fp1);

		error +=  ff;
	}

	return error;
}
