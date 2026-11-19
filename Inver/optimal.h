#ifndef _optimal_h
#define _optimal_h
///////////////////////////////////////////////////////////////////////
double Total_error(int Ns, int Nr);
void Delta_update(int col,double vt, double *Pcg, double *deltaX1, double *deltaX2);
void Gradient(int Ns, int Nr, int count1, double *Gr);
void Get_gradient(int count1, double *Gr);
void Model_update(int Ns, int Nr, double *error);
void Output_model(int step, int count1, double *deltaX2);
///////////////////////////////////////////////////////////////////////
#endif