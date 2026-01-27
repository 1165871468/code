#ifndef _matrix_h
#define _matrix_h

/////////////////////////////////////Matrix.c///////////////////////////////
void Jacobi_for(double *dx,double *dz,double x,double z,double *A11,double *A12,double *A21,double *A22);



void local_mass(short num,double *dx,double *dz,double *W,double *P,double **B);
void local_stiff(short num,long ei,double *W,double *P,double *Ax,double *Az,float **KKK,float **Kxx,float **Kxz,float **Kzx,float **Kzz);
void local_load(short num,double **node_array,long *id,double *W,double *Nx,double *Nz,double *Shape,double **di1,double **di2,
				float ****NT1,float ****NT2,double **Bxx,double **Bzz);
void Stiff_assemble(int nn,long n1,long n2,double add,long *Jtmp,double *Stmp,float *Stiff,long *Istiff,long *Jstiff);


///////////////////////////////////FILE.c///////////////////////////////////////////
double *Arraydouble1(long Nx);
double **Arraydouble2(long Nx,long Nz);
long *Arraylong1(long Nx);
double diff_fun(int Nx,int i1,int i2);
void Weight(int Nx,double *W);
void Position(int Nx,double *P);
float **Arrayfloat2(long Nx,long Nz);
float *Arrayfloat1(long Nx);
void Freefloat2(int Nx,float **Array);
void Freedouble2(int Nx,double **Array);
short *Arrayshort1(long Nx);
double Lagrange_val(short num,short ln,double *P,double lx);
double Lagrange_diff(short num,short ln,double *P,double lx);

#endif
