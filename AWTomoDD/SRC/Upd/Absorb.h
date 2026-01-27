#ifndef _Absorb_h
#define _Absorb_h

//////////////////////////////////////////////Matrix.c/////////////////////////////////////
void Stiff_assemble(int nn,long n1,long n2,double add,long *Jtmp,double *Stmp,float *Stiff,long *Istiff,long *Jstiff);


/////////////////////////////////////////////Meshcount.c////////////////////////////////////
void Ini_assemble(int Max,long n1,long n2,long *Jtmp,long *Istiff,long *Jstiff);



/////////////////////////////////////////////File.c//////////////////////////////////////////
void Weight(int Nx,double *W);
void Position(int Nx,double *P);
double diff_fun(int Nx,int i1,int i2);
long *Arraylong1(long Nx);
double *Arraydouble1(long Nx);
double **Arraydouble2(long Nx,long Nz);
int Freedouble2(int Nx,double **Array);
float **Arrayfloat2(long Nx,long Nz);
int Freefloat2(int Nx,float **Array);
float *Arrayfloat1(long Nx);

void Jacobi_for(double *dx,double *dz,double x,double z,double *A11,double *A12,double *A21,double *A22);




#endif
