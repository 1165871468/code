#ifndef _meshcount_h
#define _meshcount_h


/////////////////////////////////File.c///////////////////////////
long *Arraylong1(long Nx);
double *Arraydouble1(long Nx);
void Position(int Nx,double *P);
///////////////////////////////Meshcount.c///////////////////////////////////
void Ini_assemble(int Max,long n1,long n2,long *Jtmp,long *Istiff,long *Jstiff);


#endif
