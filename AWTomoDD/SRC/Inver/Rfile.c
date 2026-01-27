#include "Rfile.h"
#include "Common.h"

void readParForwar(int *N, int *Pflag, int *ndiv, int *nstep, int *nsnap, int *sflag, int *stype, double *dt, double *F,
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

void readParInver(int *numX, int *numY, int *numZ, int *Ns, int *Nr, int *numStep, int *numCore, double *x1, double *x2, 
    double *y1, double *y2, double *z1, double *z2, double *maxPer, double *dampFct, double *soomFct, double *extFct)
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
    sf1 = fscanf(fp1, "%lf%lf", z1, z2);             sf2 = fgets(cc1, 500, fp1);   // 3
    sf1 = fscanf(fp1, "%d%d%d", numX, numY, numZ);   sf2 = fgets(cc1, 500, fp1);   // 4
    sf1 = fscanf(fp1, "%d", Ns);                     sf2 = fgets(cc1, 500, fp1);   // 5
    sf1 = fscanf(fp1, "%d", Nr);                     sf2 = fgets(cc1, 500, fp1);   // 6
    sf1 = fscanf(fp1, "%d", numStep);                sf2 = fgets(cc1, 500, fp1);   // 7
    sf1 = fscanf(fp1, "%lf", maxPer);                sf2 = fgets(cc1, 500, fp1);   // 8
    sf1 = fscanf(fp1, "%lf", dampFct);               sf2 = fgets(cc1, 500, fp1);   // 9
    sf1 = fscanf(fp1, "%lf", soomFct);               sf2 = fgets(cc1, 500, fp1);   // 10
    sf1 = fscanf(fp1, "%d", numCore);                sf2 = fgets(cc1, 500, fp1);   // 11
    sf1 = fscanf(fp1, "%lf", extFct);                sf2 = fgets(cc1, 500, fp1);   // 12
    ////////////////////////////////////////////////////////////////////////////
    fclose(fp1);
}

void readParInver2d(int *Ns, int *Nr, int *numX, int *numZ, int *numStep, int *numCore, double *x1, double *x2,
    double *z1, double *z2, double *valPer, double *valDam, double *valSmo)
{
	int sf1;
	char cc1[500], *sf2, file[100];
    FILE *fp1;
    ////////////////////////////////////////////////////////////////////////////
    sprintf(file, "../Input/parameters_inversion.conf");
    fp1 = Open_file(file, "r");
    ////////////////////////////////////////////////////////////////////////////
    sf1 = fscanf(fp1, "%lf%lf", x1, x2);             sf2 = fgets(cc1, 500, fp1);   // 1
    sf1 = fscanf(fp1, "%lf%lf", z1, z2);             sf2 = fgets(cc1, 500, fp1);   // 2
    sf1 = fscanf(fp1, "%d%d", numX, numZ);           sf2 = fgets(cc1, 500, fp1);   // 3
    sf1 = fscanf(fp1, "%d", Ns);                     sf2 = fgets(cc1, 500, fp1);   // 4
    sf1 = fscanf(fp1, "%d", Nr);                     sf2 = fgets(cc1, 500, fp1);   // 5
    sf1 = fscanf(fp1, "%d", numStep);                sf2 = fgets(cc1, 500, fp1);   // 6
    sf1 = fscanf(fp1, "%lf", valPer);                sf2 = fgets(cc1, 500, fp1);   // 7
    sf1 = fscanf(fp1, "%lf", valDam);                sf2 = fgets(cc1, 500, fp1);   // 8
    sf1 = fscanf(fp1, "%lf", valSmo);                sf2 = fgets(cc1, 500, fp1);   // 9
    sf1 = fscanf(fp1, "%d", numCore);                sf2 = fgets(cc1, 500, fp1);   // 10
    ////////////////////////////////////////////////////////////////////////////
    fclose(fp1);
}

void VS_StrTrim(char*pStr)
{  
    char *pTmp = pStr;  
     
    while (*pStr != '\0')   
    {  
        if (*pStr != ' ')  
        {  
            *pTmp++ = *pStr;
        }  
        ++pStr;  
    }  
    *pTmp = '\0';  
}
