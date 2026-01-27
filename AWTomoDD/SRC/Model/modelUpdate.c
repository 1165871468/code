#include "Common.h"
#include "Rfile.h"

double interpolate(int n, double x, double y, double *arrfx, double *arrfz, double **arrfy);
double interpolate2D(int ix, int jy, double x, double y, double ***arrayData);
double lagrange_value(int n, int i, double x, double *data_x);
void modInforIn(int *numx, int *numz, double *x0, double *x1, double *z0, double *z1);

int main(int argc, char *argv[])
{
	int i, j, k, ix, jz, n, id1, count1, count2, flgi1, flgi2, nmx, nmz, step, fs;
	int numFor, numInv;
	char filePerResInv[100], fileVelResInv[100], fileVelIniFor[100], fileVelResFor[100];
	char fileVelResInpInv[100], fileVelIniInpInv[100], PS[250];
	double x0, z0, x1, z1, xx1, zz1, xx2, zz2, xx3, zz3, dx, dz, ff;
	double ***arrayVel1, **arrayVel2;
	double **arrayPerResInv, **arrayVelResInv, **arrayVelIniFor, **arrayVelResFor, **arrayVelResInpInv, **arrayVelIniInpInv;

	FILE *fp1, *fp2, *fp3, *fp4;

    step = 2;

	//////
	modInforIn(&nmx, &nmz, &x0, &x1, &z0, &z1);

	dx = (x1 - x0)/(nmx - 1);
	dz = (z1 - z0)/(nmz - 1);
	///
	//////
	sprintf(fileVelIniFor, "../Output/Model/Forver/vel_initial.dat");
	sprintf(fileVelResFor, "../Output/Model/Forver/vel.dat");
	sprintf(fileVelResInv, "../Output/BackUp/Result/vel%d.dat", step);
	sprintf(fileVelResInpInv, "../Output/Model/Inver/vel.dat");
	sprintf(fileVelIniInpInv, "../Output/Model/Inver/vel_initial.dat");
	///
	numFor = Count_file(fileVelIniFor);
	numInv = Count_file(fileVelResInv);

	arrayVelIniFor = Arraydouble2(numFor, 5);
	arrayVelResFor = Arraydouble2(numFor, 5);
	arrayPerResInv = Arraydouble2(numInv, 2);
	arrayVelResInv = Arraydouble2(numInv, 2);
	arrayVelResInpInv = Arraydouble2(numInv, 5);
	arrayVelIniInpInv = Arraydouble2(numInv, 5);
	arrayVel1 = Arraydouble3(nmx, nmz, 3);
	arrayVel2 = Arraydouble2(numInv, 2);


    Read_file(fileVelIniFor, numFor, 5, arrayVelIniFor);

    Read_file(fileVelResInv, numInv, 1, arrayVelResInv);

    Read_file(fileVelIniInpInv, numInv, 5, arrayVelIniInpInv);
	///
	//////
	// printf("\n%s \n", fileVelResInv);
	// id1 = 0;
	// for (i = 0; i < nmx; i++) {
	// 	for (j = 0; j < nmz; j++) {
	// 		arrayVel1[i][j][0] = x0 + i*dx;
	// 		arrayVel1[i][j][1] = z0 + j*dz;
	// 		arrayVel1[i][j][2] = arrayVelResInv[id1][0];

	// 		id1++;
	// 	}
	// }

	for (i = 1; i <= step; i++) {

		sprintf(filePerResInv, "../Output/BackUp/Result/vel_per_in%d.dat", i);

		printf("\n %s \n", filePerResInv);

		Read_file(filePerResInv, numInv, 1, arrayPerResInv);
		for (j = 0; j < numInv; j++) {

			arrayVel2[j][0] += arrayPerResInv[j][0];
		}
	}

	id1 = 0;
	for (i = 0; i < nmx; i++) {
		for (j = 0; j < nmz; j++) {
			arrayVel1[i][j][0] = x0 + i*dx;
			arrayVel1[i][j][1] = z0 + j*dz;
			arrayVel1[i][j][2] = (1.0 - arrayVel2[id1][0])*arrayVelIniInpInv[id1][2];

			id1++;
		}
	}
	///
	//////

	fp2 = fopen(fileVelResFor, "w");

	for (i = 0; i < numFor; i++) {
		xx1 = arrayVelIniFor[i][0];
		zz1 = arrayVelIniFor[i][1];

		flgi1 = 0; 
		for (j = 0; j < nmx; j++) {
			for (k = 0; k < nmz; k++) {
				xx2 = arrayVel1[j][k][0];
				zz2 = arrayVel1[j][k][1];

				xx3 = xx2 + dx;
				zz3 = zz2 + dz;

				if (xx1 >= xx2 && xx1 < xx3 && zz1 >= zz2 && zz1 < zz3 && j < nmx - 1 && k < nmz - 1) {
					ix = j;
					jz = k;
					flgi1 = 1;
					break;
				}
			}

			if (flgi1 == 1) {
				break;
			}
		}

		if (flgi1 == 1) {
			ff = interpolate2D(ix, jz, xx1, zz1, arrayVel1);
			// ff = arrayVelIniFor[i][2]*(1.0 - ff);
		}
		else {
			ff = arrayVelIniFor[i][2];
		}

		// printf("%.4f %.4f \n", ff, arrayVelIniFor[i][2]), sleep(1);

		fprintf(fp2, "%.4f %.4f %.4f %.4f %.4f\n", arrayVelIniFor[i][0], arrayVelIniFor[i][1], ff, arrayVelIniFor[i][3], arrayVelIniFor[i][4]);
	}
	fclose(fp2);
	///

	fp3 = fopen(fileVelResInpInv, "w");   // velIn

	id1 = 0;
	for (i = 0; i < nmx; i++) {
		for (j = 0; j < nmz; j++) {
			fprintf(fp3, "%.4f %.4f %.4f %.4f %.4f\n", arrayVelIniInpInv[id1][0], arrayVelIniInpInv[id1][1], arrayVel1[i][j][2], arrayVelIniInpInv[id1][3], arrayVelIniInpInv[id1][4]);
			
			id1++;
		}
	}
	fclose(fp3);

	///
	//////
	Freedouble2(numFor, arrayVelIniFor);
	Freedouble2(numFor, arrayVelResFor);
	Freedouble2(numInv, arrayPerResInv);
	Freedouble2(numInv, arrayVelResInv);
	Freedouble2(numInv, arrayVelResInpInv);
	Freedouble2(numInv, arrayVelIniInpInv);
	Freedouble2(numInv, arrayVel2);
	Freedouble3(nmx, nmz, arrayVel1);
	//////
    return 0;
}

void modInforIn(int *numx, int *numz, double *x0, double *x1, double *z0, double *z1)
{
	int i, j, ii1, ii2, sf1;
	double dd1, dd2;
	char cc1[500], flnm1[100], flnm2[100], *sf2;
	FILE *fp1;

	sprintf(flnm1, "../Input/parameters_inversion");

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
