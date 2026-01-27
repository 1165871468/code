#include "Rfile.h"
#include "Common.h"

double Read_par2d(int Nr, int Ns, double *X, double *Z, double *Rx, double *Rz, double *Sx, double *Sz)
{
	int i, id, count, Nr1, Ns1;
	double Rx1, Rz1, Sx1, Sz1, X1, Z1, **modelInfors;
	char flnm[100];

	sprintf(flnm, "../Output/Mesh/Model_infor.txt");
	count = Count_file(flnm);
	modelInfors = Arraydouble2(count, 19);
	Read_file(flnm, count, 19, modelInfors);

	id = -1;
	for (i = 0; i < count; i++)
	{
		Ns1 = modelInfors[i][0];    Nr1 = modelInfors[i][1];
		Sx1 = modelInfors[i][8];    Sz1 = modelInfors[i][9];
		Rx1 = modelInfors[i][10];   Rz1 = modelInfors[i][11];
		X1 = modelInfors[i][17];    Z1 = modelInfors[i][18];
		if (Ns == Ns1 && Nr == Nr1)
		{
			id = 1;  break;
		}
	}
	*X = X1;    *Z = Z1; 
	*Rx = Rx1;  *Rz = -Rz1;
	*Sx = Sx1;  *Sz = -Sz1;

	Freedouble2(count, modelInfors);
	return id;
}

void Read_parForwar(int *N, int *Pflag, int *ndiv, int *nstep, int *nsnap, int *sflag, int *stype, double *dt, double *F,
	double *Vp, double *Layx, double *Layz, double *R0, double *Lscale, char *moment)
{
    int N1, Pflag1, ndiv1, nsnap1, nstep1, sflag1, stype1, sf1;
	char cc1[500], file[200], *sf2;
	FILE *fp1;

	sprintf(file, "../Input/parameters_forward.config");
	fp1 = Open_file(file, "r");
    ////////////////////////////////////////////////////////////////////////////
    sf1 = fscanf(fp1, "%d", N);                      sf2 = fgets(cc1, 500, fp1);   // 1
    sf1 = fscanf(fp1, "%lf", dt);                    sf2 = fgets(cc1, 500, fp1);   // 2
    sf1 = fscanf(fp1, "%lf", F);                     sf2 = fgets(cc1, 500, fp1);   // 3
    sf1 = fscanf(fp1, "%d", nstep);                  sf2 = fgets(cc1, 500, fp1);   // 4
    sf1 = fscanf(fp1, "%d", nsnap);                  sf2 = fgets(cc1, 500, fp1);   // 5
    sf1 = fscanf(fp1, "%d", sflag);                  sf2 = fgets(cc1, 500, fp1);   // 6
    sf1 = fscanf(fp1, "%lf", Vp);                    sf2 = fgets(cc1, 500, fp1);   // 7
    sf1 = fscanf(fp1, "%lf", Layx);                  sf2 = fgets(cc1, 500, fp1);   // 8
    sf1 = fscanf(fp1, "%lf", Layz);                  sf2 = fgets(cc1, 500, fp1);   // 9
    sf1 = fscanf(fp1, "%lf", R0);                    sf2 = fgets(cc1, 500, fp1);   // 10
    sf1 = fscanf(fp1, "%d", Pflag);                  sf2 = fgets(cc1, 500, fp1);   // 11
    sf1 = fscanf(fp1, "%lf", Lscale);                sf2 = fgets(cc1, 500, fp1);   // 12
    sf1 = fscanf(fp1, "%d", ndiv);                   sf2 = fgets(cc1, 500, fp1);   // 13
    sf1 = fscanf(fp1, "%s", moment);                 sf2 = fgets(cc1, 500, fp1);   // 14
    sf1 = fscanf(fp1, "%d", stype);                  sf2 = fgets(cc1, 500, fp1);   // 15
    ////////////////////////////////////////////////////////////////////////////
	fclose(fp1);
}

void Read_parInver2d(int *numX, int *numY, int *Ns, int *Nr, int *numStep, int *numCore, double *x1, double *x2, 
    double *y1, double *y2, double *maxPer, double *dampFct, double *soomFct, double *extFct, double *waveLg)
{
    int i, j, ii1, ii2, ii3, num_core, numr, nums, numx, numy, numz, num_stp, sf1;
    double extd, dd1, dd2, deltax, deltay, damping, smooth, per_max, xx1, xx2, yy1, yy2, zz1, zz2;
    char cc1[500], flnm1[100], *sf2;
    FILE *fp1;

    sprintf(flnm1, "../Input/parameters_inversion.config");
    fp1 = Open_file(flnm1, "r");
    ////////////////////////////////////////////////////////////////////////////
    sf1 = fscanf(fp1, "%lf%lf", x1, x2);             sf2 = fgets(cc1, 500, fp1);   // 1
    sf1 = fscanf(fp1, "%lf%lf", y1, y2);             sf2 = fgets(cc1, 500, fp1);   // 2
    sf1 = fscanf(fp1, "%d%d", numX, numY);           sf2 = fgets(cc1, 500, fp1);   // 3
    sf1 = fscanf(fp1, "%d", Ns);                     sf2 = fgets(cc1, 500, fp1);   // 4
    sf1 = fscanf(fp1, "%d", Nr);                     sf2 = fgets(cc1, 500, fp1);   // 5
    sf1 = fscanf(fp1, "%d", numStep);                sf2 = fgets(cc1, 500, fp1);   // 6
    sf1 = fscanf(fp1, "%lf", maxPer);                sf2 = fgets(cc1, 500, fp1);   // 7
    sf1 = fscanf(fp1, "%lf", dampFct);               sf2 = fgets(cc1, 500, fp1);   // 8
    sf1 = fscanf(fp1, "%lf", soomFct);               sf2 = fgets(cc1, 500, fp1);   // 9
    sf1 = fscanf(fp1, "%d", numCore);                sf2 = fgets(cc1, 500, fp1);   // 10
    sf1 = fscanf(fp1, "%lf", waveLg);                sf2 = fgets(cc1, 500, fp1);   // 11
    ////////////////////////////////////////////////////////////////////////////
    fclose(fp1);
}
