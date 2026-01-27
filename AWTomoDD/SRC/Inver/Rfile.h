

void readParInver2d(int *Ns, int *Nr, int *numX, int *numZ, int *numStep, int *numCore, double *x1, double *x2,
    double *z1, double *z2, double *valPer, double *valDam, double *valSmo);

void readParForwar(int *N, int *Pflag, int *ndiv, int *nstep, int *nsnap, int *sflag, int *stype, double *dt, double *F,
	double *Vp, double *Layx, double *Layz, double *R0, double *Lscale, char *moment);

void readParInver(int *numX, int *numY, int *numZ, int *Ns, int *Nr, int *numStep, int *numCore, double *x1, double *x2, 
    double *y1, double *y2, double *z1, double *z2, double *maxPer, double *dampFct, double *soomFct, double *extFct);

void VS_StrTrim(char*pStr);
