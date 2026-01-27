#ifndef _SEM_h
#define _SEM_h
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

//////////////////////////////////////Rfile.c//////////////////////////////////
void Read_tensor(double **Tensor);
void Read_location1(double *Sx,double *Sz,char *source);
void Read_location2(long Nr,double **Rloc,char *source);
void Read_par(short *N,short *Pflag,short *ndiv,short *nstep,short *nsnap,short *sflag,short *stype,long *Nr_s,long *Nr,double *X,double *Z,double *T,double *dt,double *F,
	double *Vp,double *Layx,double *Layz,double *R0,double *Lscale,char *source,char *receiver,char *mesh,char *moment);

/////////////////////////////////////Readmesh.c///////////////////////////////
void Mesh(long *node,long *elem,char *mesh);
void Get_par(long node,long Elem,double **node_array,long **Elem_array,short *Elem_par,char *mesh);
void Count(long node,long Elem,long **Elem_array,short *Count_array,short *n1_point);
void Count_edge(long node,long elem,short max,long **Elem_array,long *Edge,long *Iedge,long *Jedge,long *edge_num);
void Count_node_csr(short num,long elem,long **Elem_array,short *Count_array,short *Count_array_add);
void Resort_node_element(short num,long elem,long **Elem_array);
void Count_material_num(short num,short Pflag,long node_add,long elem,long **Elem_array,short *Elem_par,long *node1_total,long *node_num,long *elem_num);
void Count_material_pml(short num,short Pflag,long node_add,long elem,long **Elem_array,short *Elem_par,long *node1_total,long *node_num,long *elem_num);
void node_sort(long node_add,long node_num,long *node1_total,long *node1_array);
void Elem_computation(short Pflag,long elem,short *Elem_par,long *Elem1_count);
void Elem_pml(short Pflag,long elem,short *Elem_par,long *Elem2_count);
void Collection(long node_add,long *node1_total,long *node2_total,long *node3_total,long *cnode);
void Outside_boundary(long node_add,float X,float Z,double Layx,double Layz,double len,double **node_array,
	long *node4_num,long *node4_total);
/////////////////////////////////////Meshcount.c///////////////////////////////
void re_node_elem(short num,long node,long elem,long edge_num,long **Elem_array,double **node_arrayt,double **node_array,long *Edge,long *Iedge,long *Jedge);
void Elem_in_count(short num,long elem,long elem_num,double dis,double X,double Z,long **Elem_array,double **node_array,long *elem_num1,long *Elem1_total);
void Elem_to_total(long elem,long *Elem1_total,long *Elem1_count);
void node2_count(short num,long node_add,long elem_num1,long *Elem1_count,long **Elem_array,long *node2_total,long *node2_num);
void Stiff_initial(int Node2,short max,long node2_num,long elem_num,long *Elem1_count,long **Elem_array,long *node2_total,long *Istiff,long *Jstiff);
void node_count_num(long node_add,long node_num,long *node2_total,long *node2_array);

/////////////////////////////////////Matrix.c///////////////////////////////////
double pseudo_delt(double x, double dx);

double pseudo_delt_h4(double x, double dx);

double pseudo_delt_filter(double x, double dx);

double pseudo_delt_diff(double x, double dx);

void local_mass(short int num,double *dx,double *dz,double *W,double *P,double **B);

void local_stiff(short num,long ei,double *W,double *P,double *Ax,double *Az,float **KKK,float **Kxx,float **Kxz,float **Kzx,float **Kzz);

void Find_local_position(short num,short nn,long elem_num,double dlen,double Sx,double Sz,long *Elem1_count,
	long **Elem_array,double **node_array,long *ii,double *lx,double *lz);

void Jacobi_for(double *dx,double *dz,double x,double z,double *A11,double *A12,double *A21,double *A22);

void Loadf1(short num,long elem_num,double Sx,double Sz,double **node_array,long *node1_total,long **Elem_array,long *Elem1_count,float *Loadx,float *Loadz);

void Loadf2(short num,long elem_num,double Sx,double Sz,double **node_array,long *node1_total,long **Elem_array,long *Elem1_count,float *Loadx,float *Loadz);

void Loadf(short num,long node_num,long Si,double Sx,double Sz,long *Elem1_count,long *node1_total,long **Elem_array,double **NT,double **node_array,double **M,float *Loadx,float *Loadz);

void Loadf_s(short num,long node_num,long Si,double Sx,double Sz,long *Elem1_count,long *node1_total,long **Elem_array,double **NT,double **node_array,double **M,float **Loadx_s,float **Loadz_s,long Ni);

void Load(short num,long elem_num,long node_num,double Sx,double Sz,long **Elem_array,double **node_array,float ****NT1,
		  float ****NT2,long *node1_array,long *node1_total,float *load);

void Loadf_fdm(int Nxw, int Nzw, double Sx, double Sz, double Dx, double Dz, float **Sx_fdm, float **Sz_fdm, double **FDM_array_x, double **FDM_array_z);

void Loadf_sem1(short num, long elem_num, double len, double Sx, double Sz, double **node_array, long *node1_total, long **Elem_array, long *Elem1_count,
	float *Loadx, float *Loadz, double **Tensor);

void Loadf_sem2(short num, long elem_num, double len, double Sx, double Sz, double **node_array, long *node1_total, long **Elem_array, long *Elem1_count,
	float *Loadx, float *Loadz, double **Tensor);

void Mass_matrix(short num,long elem,float *rho,double **node_array,long **Elem_array,float *Mass,float *Mass_tmp);

void Stiff_matrix(short num,long elem_num,long **Elem_array,long *Elem1_count,double **node_array,float **KKK,float **Kxx,float **Kxz,float **Kzx,float **Kzz);

void Get_Ricker(int Ni, int Nt0, double dt0, double *Ri);
///////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////File.c///////////////////////////////////
char *Arraychar1(long Nx);

long **Arraylong2(long Nx,long Nz);

float **Arrayfloat2(long Nx,long Nz);

short *Arrayshort1(long Nx);

long *Arraylong1(long Nx);

double *Arraydouble1(long Nx);

double **Arraydouble2(long Nx,long Nz);

double ***Arraydouble3(long N1, long N2, long N3);

double ****Arraydouble4(long N1, long N2, long N3, long N4);

double Lagrange_val(short num, short ln, double *P, double lx);

double Lagrange_diff(short num,short ln,double *P,double lx);

float *Arrayfloat1(long Nx);

int count_file(char *FL);

int **Arrayint2(long Nx,long Nz);

int Freeint2(int Nx, int **Array);

int Freefloat2(int Nx, float **Array);

int Freelong2(int Nx, long **Array);

int Freedouble2(int Nx,double **Array);

void Freedouble3(long N1, long N2, double ***Array);

void Freedouble4(long N1, long N2, long N3, double ****Array);

void Average_rank(int rank, int size, int Ns, int *t1, int *t2);

void Getvel(long elem,short int *elem_par,float *rho,float *lamda,float *mu);

void getvel(int L,int Nx,int Nz,float **lamda,float **mu,float **rho,double **FDM_array_x,double **FDM_array_z);

void Change_vel_sem(short num, long elem, long **Elem_array, float *rho, float *lamda,float *mu, double **node_array);

void Change_vel_sem_s1(short num, long elem, long **Elem_array, float *rho, float *lamda,float *mu, double **node_array);

void Change_vel_fdm(int L, int Nx, int Nz, float **lamda, float **mu, float **rho, double **FDM_array_x, double **FDM_array_z);

void Input_vel(long elem, float *rho, float *lamda, float *mu);

void Output_vel(short num, long node, long elem, long **Elem_array, float *rho, float *lamda, float *mu, double **node_array);

void Common_array(short num,double **NT1,double **NT2);

void Weight(int Nx,double *W);

void Position(int Nx,double *P);

double diff(int Nx,int Nz,int order,int Ni,int Nj,int d,int style,double Dx,double Dz,float **V,
			double **coeff1,double **coeff2);

double diff1(int Ni,int Nj,int style,int t1,double Dx,double Dz,float **V,double **coeff1);

double diff2(int Ni,int Nj,int style,int t2,double Dx,double Dz,float **V,double **coeff2);

void Coeffs(int N1,double **Coeff1,double **Coeff2);

void Coeff_array(int dir,int N1,double *Coeff);

void rdfile(char *FL, int count, int n, double **ary);

void Get_kernel(int Ns, int Nr, int Nt, double dt);

void multi_bin(char *fl1, char *fl2, long count, double **arrd2);

void GetRi(int Ns, int Nr, int Nt, int node_num, float **loadsemx, float **loadsemz, float **Rix_sem, float **Riz_sem);
/////////////////////////////////////Absorb.c///////////////////////////////////
void Boundary_node_inf1(short int num,long node_add,long elem_num,long elem,double X,double Z,long **Elem_array,double **node_array,
						long *node_total_pml,long *node_total_pmlc,long *node_numpml_point,long *node_numpmlc_point);
void Boundary_node_inf2(long node,long node_numpml,long node_numpmlc,long *node_total_pml,
						long *node_total_pmlc,long *node_pml,long *node_pmlc);
void node_damp(long node_numpmlc,double Vp,double Lay,double X,double Z,long *node_pmlc,double **node_array,
			   float *dampxz,float *dampxz1);
void PML_initial(int Node2,short int max,long node_numpml,long Elem_num,long Elem,long **Elem_array,long *node_total,
				 long *Istiff,long *Jstiff);
void PML_count(long node_num,long *Istiff_ini,long *Jstiff_ini,long *Istiff);

///////////////////////////////////////////////Solver.c/////////////////////////////////////////////

double Ricker(int nT,double Dt,double F);

void Output(int nm,long node_num,long *node1_array,double **node_array,double **FDM_array_x,double **FDM_array_z,float *Ux,float *Uz,float **Ufx,float **Ufz,int Nxw,int Nzw,int L,double Dx,double Dz,double Z,int order,int Nr_s);

void Output1(int nm, int Ns, int Nr, long node_num, double **node_array, float *Ux,float *Uz, float *Uxx,float *Uxz, float *Uzx,float *Uzz);

void Output1_bs(int nm1, int Ns, int Nr, long node_num, float *Uxx, float *Uxz, float *Uzx, float *Uzz);

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

void Get_record_rec(short num, short Nr, long Ni, double Z, long *Rii, float *Ux, float *Uz, long *Elem1_count, double *P, double *Rax, double *Raz,
	long **Elem_array, float **Ufx, float **Ufz, float **Syn_rex,float **Syn_rez, double **Rec_f);

void Output_record(long Nr,long Nt,double dt,float **Rx,float **Rz,int Nr_s);

void Output_record_x(long Nr,long Nt,double dt,float **Rx,float **Rz);

void Output_record_z(long Nr,long Nt,double dt,float **Rx,float **Rz);

void Find_sf_position_rec(short num, short Nr, int Nxw, int Nzw, double Z, float **Ufx, float **Ufz,
	 double **Rec_f, double **FDM_array_x, double **FDM_array_z, double **Rloc);

long Find_local_position_fdm1(long node_add,double len,double Sx,double Sz,double **node_array);

long Find_local_position_fdm2(long node_add,double len,double Sx,double Sz,double **node_array);

long Find_local_position_fdm3(long node_add,double len,double Sx,double Sz,double **node_array);

long Find_local_position_fdm4(long node_add,double len,double Sx,double Sz,double **node_array);

void Find_local_position_fdm_num(short num, short ss, int Nx, int L, long zf, double Z, long elem_num, double len, double **FDM_array_x,
	double **FDM_array_z, long *Elem1_count, long **Elem_array, double **node_array, double **Fii_fdm, double **Fax_fdm, double **Faz_fdm);
//
////void Find_local_position_fdm_pml(short num,short ss,int Nx,int L,long zf,double Z,long pelem_num,double len,double **FDM_array_x,double **FDM_array_z,long *Elem2_count,long **Elem_array,double **node_array,double **Fii_fdm,double **Fax_fdm,double **Faz_fdm);

void Find_local_position_fdm_pml(short num, int Nx, int L, long zf, double Z, long pelem_num, long pnode_num, double len,
	long *node2_array, long *Elem2_count,long **Elem_array, double **FDM_array_x, double **FDM_array_z, double **node_array,
	double **Fax_fdm, double **Faz_fdm , double ****Fid_s2f_pml);

void Find_local_position_fdm_pml_t(short num, int Nx, int L, long zf, double Z, long pelem_num, long pnode_num, double len,
	long *node2_array, long *Elem2_count, long **Elem_array, double **node_array, double **FDM_array_x, double **FDM_array_z, double ***Fid_s2f);

void Find_node_sem2fdm_cal(int order, int L, int Nxw, int Nzw, long zf, long node_num, double Dx, double Dz, double Z,
	long *node1_array, double **node_array, double **FDM_array_x, double **FDM_array_z, double **Fii_sem);

void Find_node_sem2fdm_pml(int order, int L, int Nx, int Nxw, long zf, long node_add, long pnode_num, double Dx, double Dz, double Z, 
	long *node2_array, long *node3_total, double **node_array, double **FDM_array_x, double **FDM_array_z, double **Fii_sem_pml);

void Get_record_fdm_num(short num, int Nx, int L, long zf, long *Elem1_count, float *Ux_tmp, float *Uz_tmp, double *P,
	double **Fii_fdm, long **Elem_array, double **Fax_fdm, double **Faz_fdm, float **Ufx_tmp, float **Ufz_tmp);

void Get_record_fdm_pml(short num,int Nx,int L,long zf,double **Fii_fdm,long *Elem2_count,long **Elem_array,double *P,double **Fax_fdm,double **Faz_fdm,float *Ux_tmp,float *Uz_tmp,float **Ufx_tmp,float **Ufz_tmp,double **node_array);

void Get_record_fdm_pml_c(short num, int Nx, int L, long zf, long pelem_num, long *Elem2_count, float *Upx1, float *Upx2, float *Upx3, float *Upz1, float *Upz2,
	float *Upz3, float *Ux_tmp, float *Uz_tmp,  double *P, long **Elem_array, float **Ufx_tmp, float **Ufz_tmp, double **Fax_fdm, double **Faz_fdm,
	double **node_array, float ****Ulx, float ****Ulz, float ****Urx, float ****Urz, double ****Fid_s2f_pml);

void Trans_node_fdm2sem(int order, int Nxw, int Nzw, int L, long node_num, long *node1_array, float *Ux, float *Uz, float **Ufx, float **Ufz, double Dx,
	double Dz, double Z, double **FDM_array_x, double **FDM_array_z, double **node_array, double **Fii_sem);

void Trans_node_fdm2sem_pml(int order, int Nx, int Nxw, int Nzw, int L, long node_add, long pnode_num, double Dx, double Dz, double Z, long *node2_array, 
	float *Upx1, float *Upx2, float *Upx3, float *Upz1, float *Upz2, float *Upz3, float *Ux, float *Uz, float **Ufx, float **Ufz, double **FDM_array_x, 
	double **FDM_array_z, double **node_array, double **Fii_sem_pml, float ****Ulx, float ****Ulz, float ****Urx, float ****Urz);

void Update(short ndiv,long elem,double len,short num,short nstep,short nsnap,short sflag,long Ns,long Nr,long Nt,long node_num,long elem_num,long pnode_num,long pelem_num,long node_add,long cnode_num,long bnode_num,double F,double X,
	double Z,double Sx,double Sz,double dt,double R0,double Vp,double Layx,double Layz,long *Rii,long *Elem1_count,long *Elem2_count,long *node1_total,long *node1_array,long *node2_total,long *node2_array,
	long *node3_total,long *node3_array,long *node4_array,long **Elem_array,double **node_array,double **NT,float *lam,float *mu,float *Mass,float *Loadx,float *Loadz,float **KKK,float **Kxx,
	float **Kxz,float **Kzx,float **Kzz,float **PPP,float **Pxx,float **Pxz,float **Pzx,float **Pzz,double *Rax,double *Raz,float **Loadx_s,float **Loadz_s,double **Rloc);
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
