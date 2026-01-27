#ifndef _absorb_h
#define _absorb_h
//////////////////////File.c////////////
float *Arrayfloat1(long Nx);

double *Arraydouble1(long Nx);

int **Arrayint2(long Nx,long Nz);

long **Arraylong2(long Nx,long Nz);

void Position(int Nx,double *P);

float **Arrayfloat2(long Nx,long Nz);

float ****Arrayfloat4(int N1,int N2,int N3,int N4);

double **Arraydouble2(long Nx,long Nz);

double ***Arraydouble3(long N1, long N2, long N3);

long *Arraylong1(long Nx);

double diff(int Nx,int Nz,int order,int Ni,int Nj,int d,int style,double Dx,double Dz,float **V,
			double **coeff1,double **coeff2);

double diff1(int Ni,int Nj,int style,int t1,double Dx,double Dz,float **V,double **coeff1);

double diff2(int Ni,int Nj,int style,int t2,double Dx,double Dz,float **V,double **coeff2);

double diff_f(int Nx,int Nz,int order,int Ni,int Nj,int d,int style,double Dx,double Dz,float **V,
			double **coeff1,double **coeff2,float **lamda1,float **mu1,float **rho1);

double diff1_f(int Ni,int Nj,int style,int t1,double Dx,double Dz,float **V,double **coeff1,float **lamda1,float **mu1,float **rho1);

double diff2_f(int Ni,int Nj,int style,int t2,double Dx,double Dz,float **V,double **coeff2,float **lamda1,float **mu1,float **rho1);

void Coeffs(int N1,double **Coeff1,double **Coeff2);

void Coeff_array(int dir,int N1,double *Coeff);

void getvel(int L,int Nx,int Nz,float **lamda,float **mu,float **rho,double **FDM_array_x,double **FDM_array_z);

double Lagrange_val(short num,short ln,double *P,double lx);

double Lagrange_val_1(short num,short ln,double *P,double lx);

double Lagrange_diff(short num,short ln,double *P,double lx);

////////////////////////////////Solver.c//////////////////////////
double Ricker(int nT,double Dt,double F);

void Output(int nm,long node_num,long *node1_array,double **node_array,double **FDM_array_x,double **FDM_array_z,float *Ux,float *Uz,float **Ufx,float **Ufz,int Nxw,int Nzw,int L,double Dx,double Dz,double Z,int order,int Nr_s);

void Output1_diff(int nm1,long node_num,double **node_array,double **FDM_array_x,double **FDM_array_z,float *Uxx,float *Uxz,float *Uzx,float *Uzz,float *Uxtt,float *Uztt,float **Ufxx,float **Ufxz,float **Ufzx,float **Ufzz,float **Ufxtt,float **Ufztt,int Nxw,int Nzw,int L,double len,double Dx,double Dz,double Z,int order,int Nr_s);

void Assemble(short num,long node_num,long elem_num,long *Elem1_count,long *node1_total,long **Elem_array,double *W,double **NT,float *Uxt,float *Uzt,
	float *Ux,float *Uz,float *lam,float *mu,float *Txx,float *Txz,float *Tzx,float *Tzz,float *Rx,float *Rz,float **KKK,float **Kxx,float **Kxz,
	float **Kzx,float **Kzz,double **node_array,float *Uxx,float *Uxz,float *Uzx,float *Uzz,long *node1_array);

void Assemble_pml(short num,long pnode_num,long pelem_num,long *Elem2_count,long **Elem_array,long *node2_total,double *W,float *Rx1,float *Rx2,float *Rx3,
	float *Rz1,float *Rz2,float *Rz3,float *Rxx,float *Rxz,float *Rzx,float *Rzz,float *lam,float *mu,float *Ux,float *Uz,float *Uxt,float *Uzt,float **KKK,
	float **Kxx,float **Kxz,float **Kzx,float **Kzz,double **NT,float *Txx1,float *Txx2,float *Txx3,float *Txx4,float *Txz1,float *Txz2,float *Txz3,float *Txz4,
	float *Tzx1,float *Tzx2,float *Tzx3,float *Tzx4,float *Tzz1,float *Tzz2,float *Tzz3,float *Tzz4);

void Damp_coff(short num,long pnode_num,double X,double Z,double R0,double Vp,double Layx,double Layz,long *node2_array,double **node_array,
	float *Dampx,float *Dampz,float *Dampx1,float *Dampz1);

void Get_record(short num,short Nr,long Ni,long *Rii,long *Elem1_count,long **Elem_array,double *P,double *Rax,double *Raz,float *Ux,float *Uz,float **Syn_rex,float **Syn_rez);

void Output_record(long Nr,long Nt,double dt,float **Rx,float **Rz,int Nr_s);

void Output_record_x(long Nr,long Nt,double dt,float **Rx,float **Rz);

void Output_record_z(long Nr,long Nt,double dt,float **Rx,float **Rz);

long Find_local_position_fdm1(long node_add,double len,double Sx,double Sz,double **node_array);

long Find_local_position_fdm2(long node_add,double len,double Sx,double Sz,double **node_array);

long Find_local_position_fdm3(long node_add,double len,double Sx,double Sz,double **node_array);

long Find_local_position_fdm4(long node_add,double len,double Sx,double Sz,double **node_array);

void Find_local_position_fdm_num(short num, short ss, int Nx, int L, long zf, double Z, long elem_num, double len, double **FDM_array_x,
	double **FDM_array_z, long *Elem1_count, long **Elem_array, double **node_array, double **Fii_fdm, double **Fax_fdm, double **Faz_fdm);
//
////void Find_local_position_fdm_pml(short num,short ss,int Nx,int L,long zf,double Z,long pelem_num,double len,double **FDM_array_x,double **FDM_array_z,long *Elem2_count,long **Elem_array,double **node_array,double **Fii_fdm,double **Fax_fdm,double **Faz_fdm);

void Find_node_sem2fdm_cal(int order, int L, int Nxw, int Nzw, long zf, long node_num, double Dx, double Dz, double Z,
	long *node1_array, double **node_array, double **FDM_array_x, double **FDM_array_z, double **Fii_sem);

void Find_node_sem2fdm_pml(int order, int L, int Nx, int Nxw, long zf, long node_add, long pnode_num, double Dx, double Dz, double Z, 
	long *node2_array, long *node3_total, double **node_array, double **FDM_array_x, double **FDM_array_z, double **Fii_sem_pml);

void Get_record_fdm_num(short num, int Nx, int L, long zf, long *Elem1_count, float *Ux_tmp, float *Uz_tmp, double *P,
	double **Fii_fdm, long **Elem_array, double **Fax_fdm, double **Faz_fdm, float **Ufx_tmp, float **Ufz_tmp);

void Get_record_fdm_pml(short num,int Nx,int L,long zf,double **Fii_fdm,long *Elem2_count,long **Elem_array,double *P,double **Fax_fdm,double **Faz_fdm,float *Ux_tmp,float *Uz_tmp,float **Ufx_tmp,float **Ufz_tmp,double **node_array);

void Trans_node_fdm2sem(int order, int Nxw, int Nzw, int L, long node_num, long *node1_array, float *Ux, float *Uz, float **Ufx, float **Ufz, double Dx,
	double Dz, double Z, double **FDM_array_x, double **FDM_array_z, double **node_array, double **Fii_sem);

void Trans_node_fdm2sem_pml(int order, int Nx, int Nxw, int Nzw, int L, long node_add, long pnode_num, double Dx, double Dz, double Z, long *node2_array, 
	float *Upx1, float *Upx2, float *Upx3, float *Upz1, float *Upz2, float *Upz3, float *Ux, float *Uz, float **Ufx, float **Ufz, double **FDM_array_x, 
	double **FDM_array_z, double **node_array, double **Fii_sem_pml, float ****Ulx, float ****Ulz, float ****Urx, float ****Urz);
//////////////////////////////////////////////////////////////////

void Absorbing(int Nx,int Nz,int L,int order,double Dx,double Dz,double Dt,float **Ux2,float **Uz2,float **Vx,float **Vz,
			   float **lamda,float **mu,float **rho,double **coeff1,double **coeff2,float **BPxlx,float **BPxlz,float **BPzlx,
			   float **BPzlz,float **BPxrx,float **BPxrz,float **BPzrx,float **BPzrz,float **BPxuz,float **BPzuz,float **BPxdz,
			   float **BPzdz,float ****Ulx,float ****Ulz,float ****Urx,float ****Urz,float ****Uux,float ****Uuz,float ****Udx,
			   float ****Udz);
void absorbingl(int Nx,int Nz,int L,int order,double Dx,double Dz,double Dt,double R0,float **Ux2,float **Uz2,
		float **Vx,float **Vz,float **lamda,float **mu,float **rho,double **coeff1,double **coeff2,float **BPxlx,
		float **BPxlz,float **BPzlx,float **BPzlz,float ****Ulx,float ****Ulz);
void absorbingr(int Nx,int Nz,int L,int order,double Dx,double Dz,double Dt,double R0,float **Ux2,float **Uz2,
		float **Vx,float **Vz,float **lamda,float **mu,float **rho,double **coeff1,double **coeff2,float **BPxrx,
		float **BPxrz,float **BPzrx,float **BPzrz,float ****Urx,float ****Urz);
void absorbingu(int Nx,int Nz,int L,int order,double Dx,double Dz,double Dt,double R0,float **Ux2,float **Uz2,
		float **Vx,float **Vz,float **lamda,float **mu,float **rho,double **coeff1,double **coeff2,float **BPxuz,
		float **BPzuz,float ****Uux,float ****Uuz);
void absorbingd(int Nx,int Nz,int L,int order,double Dx,double Dz,double Dt,double R0,float **Ux2,float **Uz2,
		float **Vx,float **Vz,float **lamda,float **mu,float **rho,double **coeff1,double **coeff2,float **BPxdz,
		float **BPzdz,float ****Udx,float ****Udz);
#endif
