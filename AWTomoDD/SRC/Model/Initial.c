#include "Common.h"
#include "Rfile.h"
#include "Interpolation.h"

double Lagrange_value(int n, int i, double x, double *data_x);
double Interpolate(int n, double x, double y, double *arrfx, double *arrfz, double **arrfy);

void modInforIn(int *numx, int *numz, double *x0, double *x1, double *z0, double *z1);

int main(int argc, char *argv[])
{
    int i, j, k, l, n, id1, id2, id3, id4, count1, count2, nmx, nmz, flgi1, fs;
	int num_cor, Ns, Nr, row, col, nz, step, step_max;

    char flnm1[100], flnm2[100], flnm3[100], flnm4[100], flnm5[100], flnm6[100], FF[250];

    double xx1, xx2, x0, x1, zz1, zz2, z0, z1, dx, dz, dis, dis1, rr, eta, target[2], ff[3];
    double fa, fb, fc, fd, **arrdm1, **arrayVel1, **arrayVel2, **arrayVel3;

    FILE *fp1, *fp2, *fp3;

	Vec2 p, a, b, c, d;



	sprintf(flnm1, "../Output/Model/Forver/vel_initial.dat");
	sprintf(flnm2, "../Output/Model/Forver/vel.dat");
	sprintf(flnm3, "../Output/Model/Forver/per_vel.dat");
	sprintf(flnm4, "../Output/Model/Inver/per_vel.dat");
	sprintf(flnm5, "../Output/Model/Inver/vel_initial.dat");
	sprintf(flnm6, "../Output/Model/Inver/vel.dat");
	
	sprintf(FF, " cp %s %s", flnm1, flnm2);   fs = system(FF);

	fp1 = Open_file(flnm3, "w");
	fp2 = Open_file(flnm4, "w");
	fp3 = Open_file(flnm5, "w");

	modInforIn(&nmx, &nmz, &x0, &x1, &z0, &z1);

	dx = (x1 - x0) / (nmx - 1);
	dz = (z1 - z0) / (nmz - 1);
	dis = sqrt(dx*dx + dz*dz);

	count2 = Count_file(flnm1);
	arrdm1 = Arraydouble2(count2, 5);
	arrayVel1 = Arraydouble2(count2, 3);
	arrayVel2 = Arraydouble2(count2, 3);
	arrayVel3 = Arraydouble2(count2, 3);

	Read_file(flnm1, count2, 5, arrdm1);


	for (i = 0; i < count2; i++) {
		arrayVel1[i][0] = arrdm1[i][0];
		arrayVel1[i][1] = arrdm1[i][1];
		arrayVel1[i][2] = arrdm1[i][2];

		arrayVel2[i][0] = arrdm1[i][0];
		arrayVel2[i][1] = arrdm1[i][1];
		arrayVel2[i][2] = arrdm1[i][3];

		arrayVel3[i][0] = arrdm1[i][0];
		arrayVel3[i][1] = arrdm1[i][1];
		arrayVel3[i][2] = arrdm1[i][4];
	}

	for (i = 0; i < count2; i++) {
		xx1 = arrdm1[i][0];
		zz1 = arrdm1[i][1];
		fprintf(fp1, "%.1f%10.1f%10.1f%10.1f%10.1f\n", xx1, zz1, 0.0, 0.0, 0.0);
	}

	for (i = 0; i < nmx; i++) {
		for (j = 0; j < nmz; j++) {
			xx1 = x0 + i*dx;
			zz1 = z0 + j*dz;

			fprintf(fp2, "%10.1f%10.1f%10.1f%10.1f%10.1f\n", xx1, zz1, 0.0, 0.0, 0.0);

			target[0] = xx1;
			target[1] = zz1;
			p.x = xx1;
			p.y = zz1;
/*
			id1 = Nearest_point2d_lu(target, arrayVel1, count2);
			id2 = Nearest_point2d_ru(target, arrayVel1, count2);
			id3 = Nearest_point2d_rd(target, arrayVel1, count2);
			id4 = Nearest_point2d_ld(target, arrayVel1, count2);

			if (id1 != -1 && id2 != -1 && id3 != -1 && id4 != -1 ) {
				a.x = arrayVel1[id1][0];   a.y = arrayVel1[id1][1];  fa = arrayVel1[id1][2];
				b.x = arrayVel1[id2][0];   b.y = arrayVel1[id2][1];  fb = arrayVel1[id2][2];
				c.x = arrayVel1[id3][0];   c.y = arrayVel1[id3][1];  fc = arrayVel1[id3][2];
				d.x = arrayVel1[id4][0];   d.y = arrayVel1[id4][1];  fd = arrayVel1[id4][2];
				ff[0] = Interplot_2din(p, a, b, c, d, fa, fb, fc, fd);

				a.x = arrayVel2[id1][0];   a.y = arrayVel2[id1][1];  fa = arrayVel2[id1][2];
				b.x = arrayVel2[id2][0];   b.y = arrayVel2[id2][1];  fb = arrayVel2[id2][2];
				c.x = arrayVel2[id3][0];   c.y = arrayVel2[id3][1];  fc = arrayVel2[id3][2];
				d.x = arrayVel2[id4][0];   d.y = arrayVel2[id4][1];  fd = arrayVel2[id4][2];
				ff[1] = Interplot_2din(p, a, b, c, d, fa, fb, fc, fd);

				a.x = arrayVel3[id1][0];   a.y = arrayVel3[id1][1];  fa = arrayVel3[id1][2];
				b.x = arrayVel3[id2][0];   b.y = arrayVel3[id2][1];  fb = arrayVel3[id2][2];
				c.x = arrayVel3[id3][0];   c.y = arrayVel3[id3][1];  fc = arrayVel3[id3][2];
				d.x = arrayVel3[id4][0];   d.y = arrayVel3[id4][1];  fd = arrayVel3[id4][2];
				ff[2] = Interplot_2din(p, a, b, c, d, fa, fb, fc, fd);
			}
			else {
				printf("!!! %d %d\n", i, j);
				ff[0] = Nearest_point2d(target, arrayVel1, count2);
				ff[1] = Nearest_point2d(target, arrayVel2, count2);
				ff[2] = Nearest_point2d(target, arrayVel3, count2);
			}
*/
			id1 = Nearest_point2d(target, arrayVel1, count2);
			id2 = Nearest_point2d(target, arrayVel2, count2);
			id3 = Nearest_point2d(target, arrayVel3, count2);

			ff[0] = arrayVel1[id1][2];
			ff[1] = arrayVel2[id1][2];
			ff[2] = arrayVel3[id1][2];


			fprintf(fp3, "%10.1f%10.1f%10.1f%10.1f%10.1f\n", xx1, zz1, ff[0], ff[1], ff[2]);
		}
	}



//    printf("node number: %10d%10d%10d\n", nmx, nmz, count);
////////////////////////////////////////////////////////////////
	Freedouble2(count2, arrdm1);
	Freedouble2(count2, arrayVel1);
	Freedouble2(count2, arrayVel2);
	Freedouble2(count2, arrayVel3);
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	
	sprintf(FF, " cp %s %s", flnm5, flnm6);
	fs = system(FF);
	
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

void modInforIn(int *numx, int *numz, double *x0, double *x1, double *z0, double *z1)
{
	int i, j, ii1, ii2, sf1;
	double dd1, dd2;
	char cc1[500], flnm1[100], flnm2[100], *sf2;
	FILE *fp1;

	sprintf(flnm1, "../Input/parameters_inversion.config");

	if ((fp1 = fopen(flnm1, "r")) == NULL) {
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
