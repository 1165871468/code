

double Read_par2d(int Nr, int Ns, double *X, double *Z, double *Rx, double *Rz, double *Sx, double *Sz);

void Read_parForwar(int *N, int *Pflag, int *ndiv, int *nstep, int *nsnap, int *sflag, int *stype, double *dt, double *F,
	double *Vp, double *Layx, double *Layz, double *R0, double *Lscale, char *moment);

void Read_parInver(int *numX, int *numY, int *numZ, int *Ns, int *Nr, int *numStep, int *numCore, double *x1, double *x2, 
    double *y1, double *y2, double *z1, double *z2, double *maxPer, double *dampFct, double *soomFct, double *extFct);
    
void VS_StrTrim(char*pStr);