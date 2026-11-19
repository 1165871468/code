#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "SEM.h"
#define PI 3.141592653589793238462643
#define Min 1.0e-10

//////////////////////////////////////assemble matrix/////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
void Assemble(short num,long node_num,long elem_num,long *Elem1_count,long *node1_total,long **Elem_array,double *W,double **NT,float *Uxt,float *Uzt,
	float *Ux,float *Uz,float *lam,float *mu,float *Txx,float *Txz,float *Tzx,float *Tzz,float *Rx,float *Rz,float **KKK,float **Kxx,float **Kxz,
	float **Kzx,float **Kzz,double **node_array,float *Uxx,float *Uxz,float *Uzx,float *Uzz,long *node1_array)
{
	short N2;
	long i,id,id1,i1,i2,i3,i4,j,k1,k2,t1,t2,t3,t4,dd,tt,it;
	float cof1,cof2,cof3,txx,txz,tzx,tzz,sxx,sxz,szx,szz,sxx1,sxz1,szx1,szz1,tmp,tmp1,tmp2,tmp3,tmp4;
	double Sx,Sz,*P;

	N2=num*num;

	P=Arraydouble1(num);
	Position(num-1,P);

	for(i=0;i<node_num;i++)
	{
		Rx[i] = 0.0;   Rz[i] = 0.0;
		Uxx[i] = 0.0;  Uxz[i] = 0.0;
		Uzx[i] = 0.0;  Uzz[i] = 0.0;
	}

	for(i=0;i<elem_num;i++)
	{
		id=Elem1_count[i]-1;
		cof1=lam[id];
		cof2=mu[id];
		cof3=lam[id]+2.0*mu[id];

		for(j=0;j<N2;j++)
		{
			id1=Elem_array[id][j]-1;
			Uxt[j]=Ux[id1];
			Uzt[j]=Uz[id1];
		}

		for(i1=0;i1<num;i1++)
		{
			t1=i1*num;

			for(i2=0;i2<num;i2++)
			{
				t2=t1+i2;

				tmp1=0.0; tmp2=0.0;
				tmp3=0.0; tmp4=0.0;

				for(i3=0;i3<num;i3++)
				{
					id1=i3*num+i2;
					tmp=NT[i3][i1];
					tmp1+=tmp*Uxt[id1];
					tmp2+=tmp*Uzt[id1];

					id1=i1*num+i3;
					tmp=NT[i3][i2];
					tmp3+=tmp*Uxt[id1];
					tmp4+=tmp*Uzt[id1];
				}

				dd=t1+i2;

				sxx=tmp1*Kxx[id][dd]+tmp3*Kzx[id][dd];
				sxz=tmp1*Kxz[id][dd]+tmp3*Kzz[id][dd];
				szx=tmp2*Kxx[id][dd]+tmp4*Kzx[id][dd];
				szz=tmp2*Kxz[id][dd]+tmp4*Kzz[id][dd];

				//////////////////////////difference//////////////////////////
				txx=cof3*sxx+cof1*szz;
				txz=cof2*(sxz+szx);
				tzx=txz;
				tzz=cof3*szz+cof1*sxx;
				Txx[dd]=txx*Kxx[id][dd]+txz*Kxz[id][dd];
				Txz[dd]=txx*Kzx[id][dd]+txz*Kzz[id][dd];
				Tzx[dd]=tzx*Kxx[id][dd]+tzz*Kxz[id][dd];
				Tzz[dd]=tzx*Kzx[id][dd]+tzz*Kzz[id][dd];

				tt = Elem_array[id][dd] - 1;
				tt = node1_total[tt] - 1;

				Uxx[tt] = sxx;  Uxz[tt] = sxz;
				Uzx[tt] = szx;  Uzz[tt] = szz;
			}  
		}

		for(i1=0;i1<num;i1++)
		{
			t1=i1*num;

			for(i2=0;i2<num;i2++)
			{
				tmp1=0.0; tmp3=0.0;
				tmp2=0.0; tmp4=0.0;

				for(i3=0;i3<num;i3++)
				{
					id1=i3*num+i2;
					tmp=W[i3]*KKK[id][id1]*NT[i1][i3];
					tmp1+=tmp*Txx[id1];
					tmp2+=tmp*Tzx[id1];

					id1=i1*num+i3;
					tmp=W[i3]*KKK[id][id1]*NT[i2][i3];
					tmp3+=tmp*Txz[id1];
					tmp4+=tmp*Tzz[id1];
				}

				dd = Elem_array[id][t1+i2]-1;
				dd = node1_total[dd]-1;

				Rx[dd] -= (W[i2]*tmp1+W[i1]*tmp3);
				Rz[dd] -= (W[i2]*tmp2+W[i1]*tmp4);
			}
		}
	}
}

void Assemble_pml(short num,long pnode_num,long pelem_num,long *Elem2_count,long **Elem_array,long *node2_total,double *W,float *Rx1,float *Rx2,float *Rx3,
	float *Rz1,float *Rz2,float *Rz3,float *Rxx,float *Rxz,float *Rzx,float *Rzz,float *lam,float *mu,float *Ux,float *Uz,float *Uxt,float *Uzt,float **KKK,
	float **Kxx,float **Kxz,float **Kzx,float **Kzz,double **NT,float *Txx1,float *Txx2,float *Txx3,float *Txx4,float *Txz1,float *Txz2,float *Txz3,float *Txz4,
	float *Tzx1,float *Tzx2,float *Tzx3,float *Tzx4,float *Tzz1,float *Tzz2,float *Tzz3,float *Tzz4)
{
	short N2;
	long i,j,i1,i2,i3,t1,d1,id,id1,tt,dd;
	double tmp,tmp1,tmp2,tmp3,tmp4,cof1,cof2,cof3,a11,a12,a21,a22,sxx,sxz,szx,szz,tm1,tm2,tmp0;
	double tmpx1,tmpx2,tmpx3,tmpx4,tmpx5,tmpx6,tmpx7,tmpx8,tmpz1,tmpz2,tmpz3,tmpz4,tmpz5,tmpz6,tmpz7,tmpz8;
	
	N2=num*num;
	for(i=0;i<pnode_num;i++)
	{
		Rx1[i]=0.0;  Rz1[i]=0.0;
		Rx2[i]=0.0;  Rz2[i]=0.0;
		Rx3[i]=0.0;  Rz3[i]=0.0;
		Rxx[i]=0.0;  Rzx[i]=0.0;
		Rxz[i]=0.0;  Rzz[i]=0.0;
	}
	
	for(i=0;i<pelem_num;i++)
	{
		id=Elem2_count[i]-1;
		cof1=lam[id];
		cof2=mu[id];
		cof3=lam[id]+2.0*mu[id];

		for(j=0;j<N2;j++)
		{
			id1=Elem_array[id][j]-1;
			Uxt[j]=Ux[id1];
			Uzt[j]=Uz[id1];
		}

		for(i1=0;i1<num;i1++)
		{
			t1=i1*num;
			for(i2=0;i2<num;i2++)
			{
				tmp1=0.0;  tmp2=0.0;
				tmp3=0.0;  tmp4=0.0;
				for(i3=0;i3<num;i3++)
				{
					id1=i3*num+i2;
					tmp=NT[i3][i1];
					tmp1+=tmp*Uxt[id1];
					tmp2+=tmp*Uzt[id1];

					id1=i1*num+i3;
					tmp=NT[i3][i2];
					tmp3+=tmp*Uxt[id1];
					tmp4+=tmp*Uzt[id1];
				}
				d1=t1+i2;
				a11=Kxx[i][d1]; a12=Kzx[i][d1];
				a21=Kxz[i][d1]; a22=Kzz[i][d1];
				
				sxx=tmp1*a11+tmp3*a12;
				sxz=tmp1*a21+tmp3*a22;
				szx=tmp2*a11+tmp4*a12;
				szz=tmp2*a21+tmp4*a22;

				//////////////////////////P/////////////////////////////
				tt=Elem_array[id][d1]-1;
				tt=node2_total[tt]-1;
				Rxx[tt]-=sxx*cof3;
				Rxz[tt]-=sxz*cof2;
				Rzx[tt]-=szx*cof2;
				Rzz[tt]-=szz*cof3;
				/////////////////////Internal-Process/////////////////
				Txx1[d1]=sxx*a11; Txx2[d1]=sxx*a12;
				Txx3[d1]=sxx*a21; Txx4[d1]=sxx*a22;
				
				Txz1[d1]=sxz*a11; Txz2[d1]=sxz*a12;
				Txz3[d1]=sxz*a21; Txz4[d1]=sxz*a22;
				
				Tzx1[d1]=szx*a11; Tzx2[d1]=szx*a12;
				Tzx3[d1]=szx*a21; Tzx4[d1]=szx*a22;
				
				Tzz1[d1]=szz*a11; Tzz2[d1]=szz*a12;
				Tzz3[d1]=szz*a21; Tzz4[d1]=szz*a22;
			}
		}

		for(i1=0;i1<num;i1++)
		{
			t1=i1*num;
			for(i2=0;i2<num;i2++)
			{
				tmpx1=0.0;  tmpx2=0.0;
				tmpx3=0.0;  tmpx4=0.0;
				tmpx5=0.0;  tmpx6=0.0;
				tmpx7=0.0;  tmpx8=0.0;
				
				tmpz1=0.0;  tmpz2=0.0;
				tmpz3=0.0;  tmpz4=0.0;
				tmpz5=0.0;  tmpz6=0.0;
				tmpz7=0.0;  tmpz8=0.0;

				for(i3=0;i3<num;i3++)
				{
					d1=i3*num+i2;
					tmp0=W[i3]*KKK[i][d1];
					tmp=tmp0*NT[i1][i3];
					tmpx1+=Txx1[d1]*tmp;  tmpx2+=Txx3[d1]*tmp;					
					tmpx3+=Txz1[d1]*tmp;  tmpx4+=Txz3[d1]*tmp;
					tmpx5+=Tzx1[d1]*tmp;  tmpx6+=Tzx3[d1]*tmp;
					tmpx7+=Tzz1[d1]*tmp;  tmpx8+=Tzz3[d1]*tmp;
					
					d1=i1*num+i3;
					tmp0=W[i3]*KKK[i][d1];
					tmp=tmp0*NT[i2][i3];
					tmpz1+=Txx2[d1]*tmp;  tmpz2+=Txx4[d1]*tmp;
					tmpz3+=Txz2[d1]*tmp;  tmpz4+=Txz4[d1]*tmp;
					tmpz5+=Tzx2[d1]*tmp;  tmpz6+=Tzx4[d1]*tmp;
					tmpz7+=Tzz2[d1]*tmp;  tmpz8+=Tzz4[d1]*tmp;
				}
				tm1=W[i1];
				tm2=W[i2];
				
				tt=Elem_array[id][t1+i2]-1;	
				tt=node2_total[tt]-1;

				Rx1[tt]-=(tm2*tmpx1+tm1*tmpz1)*cof3;
				Rx2[tt]-=(tm2*tmpx7+tm1*tmpz7)*cof1+(tm2*tmpx6+tm1*tmpz6)*cof2;
				Rx3[tt]-=(tm2*tmpx4+tm1*tmpz4)*cof2;

				Rz1[tt]-=(tm2*tmpx5+tm1*tmpz5)*cof2;
				Rz2[tt]-=(tm2*tmpx2+tm1*tmpz2)*cof1+(tm2*tmpx3+tm1*tmpz3)*cof2;
				Rz3[tt]-=(tm2*tmpx8+tm1*tmpz8)*cof3;
			}   
		}
	}
}

void Collect_treat(long cnode_num,long *node1_total,long *node2_total,long *node3_array,float *Rx1,float *Rx2,float *Rx3,float *Rz1,float *Rz2,float *Rz3,float *Rx,float *Rz)
{
	long i,i1,i2,id;
	for(i=0;i<cnode_num;i++)
	{
		id=node3_array[i]-1;
		i1=node1_total[id]-1;
		i2=node2_total[id]-1;

		Rx[i1]+=Rx1[i2]+Rx2[i2]+Rx3[i2];
		Rz[i1]+=Rz1[i2]+Rz2[i2]+Rz3[i2];
	}	
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////Update///////////////////////////////////////////////////////
void UpdateSyn(short ndiv,long elem,double len,short num,short nstep,short nsnap,short sflag,long Nr,long Ns,long Nt,long node_num,long elem_num,long pnode_num,long pelem_num,long node_add,long cnode_num,long bnode_num,
	double F,double X,double Z,double Sx,double Sz,double dt,double R0,double Vp,double Layx,double Layz,long *Rii,long *Elem1_count,long *Elem2_count,long *node1_total,long *node1_array,long *node2_total,long *node2_array,
	long *node3_total,long *node3_array,long *node4_array,long **Elem_array,double **node_array,double **NT,float *lam,float *mu,float *Mass,float *Mass_tmp,float *Loadx,float *Loadz,float *arrayUusyn,float **KKK,float **Kxx,
	float **Kxz,float **Kzx,float **Kzz,float **PPP,float **Pxx,float **Pxz,float **Pzx,float **Pzz,double *Rax,double *Raz,int Nr_s,float **Loadx_s,float **Loadz_s,double **Rloc)
{
	short N2, ss;

	int M, N, Nx, Nz, Nxw, Nzw, Nx1, Nz1, L, node_fdm, order, order1;
	int t1, t2, t3, t4, nm, nm1;

	long i, i1, i2, i3, i4, id, it, ix, iz, index1, index2, j, k, k1, k2, d1, d2, z, zf, Ni;
	long numIntervalt, numUu;
	long **fx1, **fx2, **fx3, **fx4, **fz1, **fz2, **fz3, **fz4;

	float tmp, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp_x, tmp_z, Loadx_tmp, Loadz_tmp;
	float *Dampx, *Dampz, *Dampx1, *Dampz1, *Plxx, *Plxz, *Plzx, *Plzz, *Rxx, *Rzx, *Rxz, *Rzz, *Upx1, *Upx2, *Upx3;
	float *Utx1, *Utx2, *Utx3, *Rx1, *Rx2, *Rx3, *Upz1, *Upz2, *Upz3, *Utz1, *Utz2, *Utz3, *Rz1, *Rz2, *Rz3;
	float *Ux, *Uz, *Ux_tmp, *Uz_tmp, *Uxt, *Uzt, *Ux1, *Uz1, *Rx, *Rz, *Txx1, *Txx2, *Txx3, *Txx4, *Txz1, *Txz2, *Txz3, *Txz4;
	float *Tzx1, *Tzx2, *Tzx3, *Tzx4, *Tzz1, *Tzz2, *Tzz3, *Tzz4, *Txx, *Txz, *Tzx, *Tzz, **Syn_rex, **Syn_rez;
	float **lamda1, **mu1, **rho1, **Ufx, **Ufz, **Vfx, **Vfz, **tmpx, **tmpz, **Ufx_tmp, **Ufz_tmp, **Sx_fdm, **Sz_fdm;
	float **BPxlx, **BPxlz, **BPzlx, **BPzlz, **BPxrx, **BPxrz, **BPzrx, **BPzrz, **BPxuz, **BPzuz, **BPxdz, **BPzdz;
	float *Uxx, *Uxz, *Uzx, *Uzz, *Uxtt, *Uztt, *Uxx_tmp, *Uxz_tmp, *Uzx_tmp, *Uzz_tmp;
	float **Ufxx, **Ufxz, **Ufzx, **Ufzz, **Ufxtt, **Ufztt, **Syn_uxx, **Syn_uzx, **Syn_uxz, **Syn_uzz;
	float ****Ulx, ****Ulz, ****Urx, ****Urz, ****Uux, ****Uuz, ****Udx, ****Udz;
	float *arrayUu;

	double dt2, cm, dx, dz, dx1, dz1, dxx, dxz, dzx, dzz, dis, xxx, zzz, t, tt, xx, zz;
	double Dx, Dz, X1, X0, Z0, Z1, Ds, S, Ba, Ri, Stress;
	double x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, df1, df2, df3, df4;
	double r1, r2, r3, r4, rr1, rr2, rr3, rr4;
	double tmpf1, tmpf2, tmpf3, tmpf4, tmpf11, tmpf12, tmpf21, tmpf22, valx, valz;
	double  *Sem_ri, *P, *W, *ttx, *ttz, **coeff1, **coeff2;
	double **Fax_fdm, **Faz_fdm, **Fii_fdm, **Fii_sem, **Fii_sem_pml, **Fii, **Fax, **Faz, **FDM_array_x, **FDM_array_z, **Rec_f;
	double ***Fid_s2f, ****Fid_s2f_pml;

	FILE *fp;

	numIntervalt = 5;
	numUu = 4*node_add*(Nt / numIntervalt);

	// printf("%ld %ld %ld\n", Nt,  Nt / numIntervalt,  numUu); sleep(1000);

	nm = 0; nm1 = 0;

	dt2 = dt * dt;
	N2 = num * num;

	P = Arraydouble1(num);
	W = Arraydouble1(num);

	Position(num - 1, P);
	Weight(num - 1, W);

	//////////////////////////PML//////////////////////////////
	Dampx = Arrayfloat1(pnode_num);  Dampx1 = Arrayfloat1(pnode_num);
	Dampz = Arrayfloat1(pnode_num);  Dampz1 = Arrayfloat1(pnode_num);

	Rxx = Arrayfloat1(pnode_num);    Rzx = Arrayfloat1(pnode_num);
	Rxz = Arrayfloat1(pnode_num);    Rzz = Arrayfloat1(pnode_num);

	Plxx = Arrayfloat1(pnode_num);   Plxz = Arrayfloat1(pnode_num);
	Plzx = Arrayfloat1(pnode_num);   Plzz = Arrayfloat1(pnode_num);
	
	Upx1 = Arrayfloat1(pnode_num);   Utx1 = Arrayfloat1(pnode_num);   Rx1 = Arrayfloat1(pnode_num);
	Upx2 = Arrayfloat1(pnode_num);   Utx2 = Arrayfloat1(pnode_num);   Rx2 = Arrayfloat1(pnode_num);
	Upx3 = Arrayfloat1(pnode_num);   Utx3 = Arrayfloat1(pnode_num);   Rx3 = Arrayfloat1(pnode_num);

	Upz1 = Arrayfloat1(pnode_num);   Utz1 = Arrayfloat1(pnode_num);   Rz1 = Arrayfloat1(pnode_num);
	Upz2 = Arrayfloat1(pnode_num);   Utz2 = Arrayfloat1(pnode_num);   Rz2 = Arrayfloat1(pnode_num);
	Upz3 = Arrayfloat1(pnode_num);   Utz3 = Arrayfloat1(pnode_num);   Rz3 = Arrayfloat1(pnode_num);

	Txx1 = Arrayfloat1(N2); Txz1 = Arrayfloat1(N2); Tzx1 = Arrayfloat1(N2); Tzz1 = Arrayfloat1(N2);
	Txx2 = Arrayfloat1(N2); Txz2 = Arrayfloat1(N2); Tzx2 = Arrayfloat1(N2); Tzz2 = Arrayfloat1(N2);
	Txx3 = Arrayfloat1(N2); Txz3 = Arrayfloat1(N2); Tzx3 = Arrayfloat1(N2); Tzz3 = Arrayfloat1(N2);
	Txx4 = Arrayfloat1(N2); Txz4 = Arrayfloat1(N2); Tzx4 = Arrayfloat1(N2); Tzz4 = Arrayfloat1(N2);

	Damp_coff(num,pnode_num,X,Z,R0,Vp,Layx,Layz,node2_array,node_array,Dampx,Dampz,Dampx1,Dampz1);

	////////////////////////// Computation //////////////////////////
	Ux = Arrayfloat1(node_add);         Uxt = Arrayfloat1(N2);  
	Uz = Arrayfloat1(node_add);         Uzt = Arrayfloat1(N2);
	Ux1 = Arrayfloat1(node_num);        Rx = Arrayfloat1(node_num);
	Uz1 = Arrayfloat1(node_num);        Rz = Arrayfloat1(node_num);
	Txx = Arrayfloat1(N2);              Txz = Arrayfloat1(N2);
	Tzx = Arrayfloat1(N2);              Tzz = Arrayfloat1(N2);
	Syn_rex = Arrayfloat2(Nr, Nt);      Syn_rez = Arrayfloat2(Nr, Nt);
	Syn_uxx = Arrayfloat2(Nr, Nt);      Syn_uzx = Arrayfloat2(Nr, Nt);
	Syn_uxz = Arrayfloat2(Nr, Nt);      Syn_uzz = Arrayfloat2(Nr, Nt);
	Sem_ri = Arraydouble1(Nt);
	Ux_tmp = Arrayfloat1(node_add);     Uz_tmp = Arrayfloat1(node_add);

/////////////////////////////////// Partial derivative ///////////////////////////////////
	Uxx = Arrayfloat1(node_add);        Uxz = Arrayfloat1(node_add);
	Uzx = Arrayfloat1(node_add);        Uzz = Arrayfloat1(node_add);
	Uxtt = Arrayfloat1(node_add);       Uztt = Arrayfloat1(node_add);
	Uxx_tmp = Arrayfloat1(node_num);    Uxz_tmp = Arrayfloat1(node_num);
	Uzx_tmp = Arrayfloat1(node_num);    Uzz_tmp = Arrayfloat1(node_num);
////////////////////////////////// Output array --totall //////////////////////////////////
	arrayUu = Arrayfloat1(numUu);

//	Get_Ricker(Ni, Nt, dt, Sem_ri);
////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////// Calculation ////////////////////////////////////////
	for(Ni = 1; Ni <= Nt; Ni++)
	{
		Assemble(num, node_num, elem_num, Elem1_count, node1_total, Elem_array, W, NT, Uxt, Uzt, Ux, Uz, lam, mu, Txx, Txz, Tzx, Tzz, 
			Rx, Rz, KKK, Kxx, Kxz, Kzx, Kzz, node_array, Uxx_tmp, Uxz_tmp, Uzx_tmp, Uzz_tmp, node1_array);

		Assemble_pml(num, pnode_num, pelem_num, Elem2_count, Elem_array, node2_total, W, Rx1, Rx2, Rx3, Rz1, Rz2, Rz3, Rxx, Rxz, Rzx, 
			Rzz, lam, mu, Ux, Uz, Uxt, Uzt, PPP, Pxx, Pxz, Pzx, Pzz, NT, Txx1, Txx2, Txx3, Txx4, Txz1, Txz2, Txz3, Txz4, Tzx1, Tzx2, Tzx3, 
			Tzx4, Tzz1, Tzz2, Tzz3, Tzz4);

		Collect_treat(cnode_num, node1_total, node2_total, node3_array, Rx1, Rx2, Rx3, Rz1, Rz2, Rz3, Rx, Rz);

		Ri = Ricker(Ni - 1, dt, F);

		for(i = 0; i < node_num; i++)
		{
			i1 = node1_array[i]-1;
			tmp = Ux[i1];
			Ux[i1] = 2.0*tmp - Ux1[i] + (Rx[i] + Loadx_s[0][i]*Ri)*dt2*Mass[i1];
			Ux1[i] = tmp;
			Uxtt[i1] = (Ux[i1] - tmp)/dt;

			tmp = Uz[i1];
			Uz[i1] = 2.0*tmp - Uz1[i] + (Rz[i] + Loadz_s[0][i]*Ri)*dt2*Mass[i1];
			Uz1[i] = tmp;
			Uztt[i1] = (Uz[i1] - tmp)/dt;

			Ux_tmp[i1] = Ux1[i];
			Uz_tmp[i1] = Uz1[i];

			Uxx[i1] = Uxx_tmp[i];
			Uxz[i1] = Uxz_tmp[i];
			Uzx[i1] = Uzx_tmp[i];
			Uzz[i1] = Uzz_tmp[i];
		}

		for(i = 0;i < bnode_num; i++)
		{
			i1 = node4_array[i] - 1;
			i1 = node2_total[i1] - 1;
			Rx1[i1] = 0.0;  Rz1[i1] = 0.0; 
			Rx2[i1] = 0.0;  Rz2[i1] = 0.0;
			Rx3[i1] = 0.0;  Rz3[i1] = 0.0;
			Rxx[i1] = 0.0;  Rxz[i1] = 0.0;
			Rzx[i1] = 0.0;  Rzz[i1] = 0.0;
		}

		for(i = 0; i < pnode_num; i++)
		{
			i1 = node2_array[i] - 1;

			if(node3_total[i1] > 0)
			{
				continue;
			}

			cm = Mass[i1];
			dx = Dampx[i]; dx1 = Dampx1[i];
			dz = Dampz[i]; dz1 = Dampz1[i];

			dxx = dx*dx;  dxz = dx*dz;
			dzx = dz*dx;  dzz = dz*dz;
			//////////////////////////////////Ux////////////////////////////////////
			tt = Rx1[i] * cm + Plxx[i];
			tmp = Upx1[i];
			Upx1[i] = ((2.0 - dt2*dxx)*Upx1[i] + (dt*dx - 1.0)*Utx1[i] + dt2*tt) / (1.0 + dt*dx);
			Utx1[i] = tmp;

			tt = Rx2[i] * cm;
			tmp = Upx2[i];
			Upx2[i] = ((2.0 - dt2*dxz)*Upx2[i] + (0.5*dt*(dx + dz) - 1.0)*Utx2[i] + dt2*tt) / (1.0 + 0.5*dt*(dx + dz));
			Utx2[i] = tmp;

			tt = Rx3[i] * cm + Plxz[i];
			tmp = Upx3[i];
			Upx3[i] = ((2.0 - dt2*dzz)*Upx3[i] + (dt*dz - 1.0)*Utx3[i] + dt2*tt) / (1.0 + dt*dz);
			Utx3[i] = tmp;

			Plxx[i] = (1.0 - dt*dx)*Plxx[i] + dt*dx1*Rxx[i]*cm;
			Plxz[i] = (1.0 - dt*dz)*Plxz[i] + dt*dz1*Rxz[i]*cm;
			///////////////////////////////Uz//////////////////////////////////////
			tt = Rz1[i] * cm + Plzx[i];
			tmp = Upz1[i];
			Upz1[i] = ((2.0 - dt2*dxx)*Upz1[i] + (dt*dx - 1.0)*Utz1[i] + dt2*tt) / (1.0 + dt*dx);
			Utz1[i] = tmp;

			tt = Rz2[i] * cm;
			tmp = Upz2[i];
			Upz2[i] = ((2.0 - dt2*dzx) * Upz2[i] + (0.5*dt*(dx + dz) - 1.0)*Utz2[i] + dt2*tt) / (1.0 + 0.5*dt*(dx + dz));
			Utz2[i] = tmp;

			tt = Rz3[i] * cm + Plzz[i];
			tmp = Upz3[i];
			Upz3[i] = ((2.0 - dt2*dzz)*Upz3[i] + (dt*dz - 1.0)*Utz3[i] + dt2*tt) / (1.0 + dt*dz);
			Utz3[i] = tmp;

			Plzx[i] = (1.0 - dt*dx)*Plzx[i] + dt*dx1*Rzx[i]*cm;
			Plzz[i] = (1.0 - dt*dz)*Plzz[i] + dt*dz1*Rzz[i]*cm;
			///////////////////////////////Uz////////////////////////////////////////
			Ux[i1] = Upx1[i] + Upx2[i] + Upx3[i];
			Uz[i1] = Upz1[i] + Upz2[i] + Upz3[i];
			Ux_tmp[i1] = Utx1[i] + Utx2[i] + Utx3[i];
			Uz_tmp[i1] = Utz1[i] + Utz2[i] + Utz3[i];
		}

		if(Ni % nsnap == 0)
		{
			// Output1_syn(nm, Ns, Nr, node_num, node_array, Ux, Uz, Uxx, Uxz, Uzx, Uzz);
			nm++;
		}
////////////////////////////////////////
		if(Ni % numIntervalt == 0)
		{
			nm1++;
			// Output1_bs(nm1, Ns, Nr, node_add, Uxx, Uxz, Uzx, Uzz);
			GetUu_bs(nm1, node_add, Uxx, Uxz, Uzx, Uzz, arrayUu);
		}

		Get_record(num,1 ,Ni, Rii, Elem1_count, Elem_array, P, Rax, Raz, Ux, Uz, Syn_rex, Syn_rez);
	}

	Output_record_syn(Nr, Ns, Nt, dt, Syn_rex, Syn_rez, Nr_s);
	// Output_record(Nr, Ns, Nt, dt, Syn_rex, Syn_rez, Nr_s);

	for (i = 0; i < numUu; i++)
	{
		arrayUusyn[i] = arrayUu[i];
	}

	//Outputotall_bs_syn(nm1, Ns, Nr, numUu, arrayUu);
	///////////////////////FREE ARRAYS/////////////////////////////////////////
	free(P);      free(W);
	free(Dampx);  free(Dampz);
	free(Dampx1); free(Dampz1);
	free(Rxx);    free(Rzx);
	free(Rxz);    free(Rzz);
	free(Plxx);   free(Plzx);
	free(Plxz);   free(Plzz);
	free(Upx1);   free(Utx1);   free(Rx1);
	free(Upx2);   free(Utx2);   free(Rx2);
	free(Upx3);   free(Utx3);   free(Rx3);
	free(Upz1);   free(Utz1);   free(Rz1);
	free(Upz2);   free(Utz2);   free(Rz2);
	free(Upz3);   free(Utz3);   free(Rz3);

	free(Txx1);   free(Txz1);   free(Tzx1);   free(Tzz1);
	free(Txx2);   free(Txz2);   free(Tzx2);   free(Tzz2);
	free(Txx3);   free(Txz3);   free(Tzx3);   free(Tzz3);
	free(Txx4);   free(Txz4);   free(Tzx4);   free(Tzz4);

	free(Ux);     free(Uxt);    free(Ux1);    free(Rx);
	free(Uz);     free(Uzt);    free(Uz1);    free(Rz);
	
	free(Txx);    free(Txz);    free(arrayUu);
	free(Tzx);    free(Tzz);    

	Freefloat2(Nr,Syn_rex);
	Freefloat2(Nr,Syn_rez);
	Freefloat2(Nr,Syn_uxx);
	Freefloat2(Nr,Syn_uzx);
	Freefloat2(Nr,Syn_uxz);
	Freefloat2(Nr,Syn_uzz);
	// printf("ccc\n"), sleep(1000);

	free(Ux_tmp); free(Uz_tmp);
	free(Uxx);    free(Uxz);      free(Uzx);     free(Uzz);      free(Uxtt);
	free(Uztt);   free(Uxx_tmp);  free(Uxz_tmp); free(Uzx_tmp);  free(Uzz_tmp);
	
}


void UpdateAdj(short ndiv,long elem,double len,short num,short nstep,short nsnap,short sflag,long Nr,long Ns,long Nt,long node_num,long elem_num,long pnode_num,long pelem_num,long node_add,long cnode_num,long bnode_num,
	double F,double X,double Z,double Sx,double Sz,double dt,double R0,double Vp,double Layx,double Layz,long *Rii,long *Elem1_count,long *Elem2_count,long *node1_total,long *node1_array,long *node2_total,long *node2_array,
	long *node3_total,long *node3_array,long *node4_array,long **Elem_array,double **node_array,double **NT,float *lam,float *mu,float *Mass,float *Mass_tmp,float *Loadx,float *Loadz,float *arrayUuadj,float **KKK,float **Kxx,
	float **Kxz,float **Kzx,float **Kzz,float **PPP,float **Pxx,float **Pxz,float **Pzx,float **Pzz,double *Rax,double *Raz,int Nr_s,float **Loadx_s,float **Loadz_s,double **Rloc)
{
	short N2, ss;

	int M, N, Nx, Nz, Nxw, Nzw, Nx1, Nz1, L, node_fdm, order, order1;
	int t1, t2, t3, t4, numIntervalt, nm, nm1;

	long i, i1, i2, i3, i4, id, it, ix, iz, index1, index2, j, k, k1, k2, d1, d2, z, zf, numUu, Ni;
	long **fx1, **fx2, **fx3, **fx4, **fz1, **fz2, **fz3, **fz4;

	float tmp, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp_x, tmp_z, Loadx_tmp, Loadz_tmp;
	float *Dampx, *Dampz, *Dampx1, *Dampz1, *Plxx, *Plxz, *Plzx, *Plzz, *Rxx, *Rzx, *Rxz, *Rzz, *Upx1, *Upx2, *Upx3;
	float *Utx1, *Utx2, *Utx3, *Rx1, *Rx2, *Rx3, *Upz1, *Upz2, *Upz3, *Utz1, *Utz2, *Utz3, *Rz1, *Rz2, *Rz3;
	float *Ux, *Uz, *Ux_tmp, *Uz_tmp, *Uxt, *Uzt, *Ux1, *Uz1, *Rx, *Rz, *Txx1, *Txx2, *Txx3, *Txx4, *Txz1, *Txz2, *Txz3, *Txz4;
	float *Tzx1, *Tzx2, *Tzx3, *Tzx4, *Tzz1, *Tzz2, *Tzz3, *Tzz4, *Txx, *Txz, *Tzx, *Tzz, **Syn_rex, **Syn_rez;
	float **lamda1, **mu1, **rho1, **Ufx, **Ufz, **Vfx, **Vfz, **tmpx, **tmpz, **Ufx_tmp, **Ufz_tmp, **Sx_fdm, **Sz_fdm;
	float **BPxlx, **BPxlz, **BPzlx, **BPzlz, **BPxrx, **BPxrz, **BPzrx, **BPzrz, **BPxuz, **BPzuz, **BPxdz, **BPzdz;
	float *Uxx, *Uxz, *Uzx, *Uzz, *Uxtt, *Uztt, *Uxx_tmp, *Uxz_tmp, *Uzx_tmp, *Uzz_tmp;
	float **Ufxx, **Ufxz, **Ufzx, **Ufzz, **Ufxtt, **Ufztt, **Syn_uxx, **Syn_uzx, **Syn_uxz, **Syn_uzz;
	float ****Ulx, ****Ulz, ****Urx, ****Urz, ****Uux, ****Uuz, ****Udx, ****Udz;
	float **Rix_sem, **Riz_sem, *arrayUu;

	double dt2, cm, dx, dz, dx1, dz1, dxx, dxz, dzx, dzz, dis, xxx, zzz, t, tt, xx, zz;
	double Dx, Dz, X1, X0, Z0, Z1, Ds, S, Ba, Ri, Stress;
	double x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, df1, df2, df3, df4;
	double r1, r2, r3, r4, rr1, rr2, rr3, rr4;
	double tmpf1, tmpf2, tmpf3, tmpf4, tmpf11, tmpf12, tmpf21, tmpf22, valx, valz;
	double  *Sem_ri, *P, *W, *ttx, *ttz, **coeff1, **coeff2;
	double **Fax_fdm, **Faz_fdm, **Fii_fdm, **Fii_sem, **Fii_sem_pml, **Fii, **Fax, **Faz, **FDM_array_x, **FDM_array_z, **Rec_f;
	double ***Fid_s2f, ****Fid_s2f_pml;

	FILE *fp;
	/////////////////////////////////////
	nm = 0; nm1 = 0;
	/////////////////////////////////////
	numIntervalt = 5;
	numUu = 4*node_add*(Nt / numIntervalt);
	/////////////////////////////////////

	dt2 = dt*dt;
	N2 = num*num;

	P = Arraydouble1(num);
	W = Arraydouble1(num);

	Position(num-1, P);
	Weight(num-1, W);

	//////////////////////////PML//////////////////////////////
	Dampx = Arrayfloat1(pnode_num);  Dampx1 = Arrayfloat1(pnode_num);
	Dampz = Arrayfloat1(pnode_num);  Dampz1 = Arrayfloat1(pnode_num);

	Rxx = Arrayfloat1(pnode_num);    Rzx = Arrayfloat1(pnode_num);
	Rxz = Arrayfloat1(pnode_num);    Rzz = Arrayfloat1(pnode_num);

	Plxx = Arrayfloat1(pnode_num);   Plxz = Arrayfloat1(pnode_num);
	Plzx = Arrayfloat1(pnode_num);   Plzz = Arrayfloat1(pnode_num);
	
	Upx1 = Arrayfloat1(pnode_num);   Utx1 = Arrayfloat1(pnode_num);   Rx1 = Arrayfloat1(pnode_num);
	Upx2 = Arrayfloat1(pnode_num);   Utx2 = Arrayfloat1(pnode_num);   Rx2 = Arrayfloat1(pnode_num);
	Upx3 = Arrayfloat1(pnode_num);   Utx3 = Arrayfloat1(pnode_num);   Rx3 = Arrayfloat1(pnode_num);

	Upz1 = Arrayfloat1(pnode_num);   Utz1 = Arrayfloat1(pnode_num);   Rz1 = Arrayfloat1(pnode_num);
	Upz2 = Arrayfloat1(pnode_num);   Utz2 = Arrayfloat1(pnode_num);   Rz2 = Arrayfloat1(pnode_num);
	Upz3 = Arrayfloat1(pnode_num);   Utz3 = Arrayfloat1(pnode_num);   Rz3 = Arrayfloat1(pnode_num);

	Txx1 = Arrayfloat1(N2); Txz1 = Arrayfloat1(N2); Tzx1 = Arrayfloat1(N2); Tzz1 = Arrayfloat1(N2);
	Txx2 = Arrayfloat1(N2); Txz2 = Arrayfloat1(N2); Tzx2 = Arrayfloat1(N2); Tzz2 = Arrayfloat1(N2);
	Txx3 = Arrayfloat1(N2); Txz3 = Arrayfloat1(N2); Tzx3 = Arrayfloat1(N2); Tzz3 = Arrayfloat1(N2);
	Txx4 = Arrayfloat1(N2); Txz4 = Arrayfloat1(N2); Tzx4 = Arrayfloat1(N2); Tzz4 = Arrayfloat1(N2);

	Damp_coff(num,pnode_num,X,Z,R0,Vp,Layx,Layz,node2_array,node_array,Dampx,Dampz,Dampx1,Dampz1);

	////////////////////////// Computation //////////////////////////
	Ux = Arrayfloat1(node_add);          Uxt = Arrayfloat1(N2);  
	Uz = Arrayfloat1(node_add);          Uzt = Arrayfloat1(N2);
	Ux1 = Arrayfloat1(node_num);         Rx = Arrayfloat1(node_num);
	Uz1 = Arrayfloat1(node_num);         Rz = Arrayfloat1(node_num);
	Txx = Arrayfloat1(N2);               Txz = Arrayfloat1(N2);
	Tzx = Arrayfloat1(N2);               Tzz = Arrayfloat1(N2);
	Syn_rex = Arrayfloat2(Nr,Nt);        Syn_rez = Arrayfloat2(Nr,Nt);
	Syn_uxx = Arrayfloat2(Nr,Nt);        Syn_uzx = Arrayfloat2(Nr,Nt);
	Syn_uxz = Arrayfloat2(Nr,Nt);        Syn_uzz = Arrayfloat2(Nr,Nt);
	Sem_ri = Arraydouble1(Nt);

	Ux_tmp = Arrayfloat1(node_add);      Uz_tmp = Arrayfloat1(node_add);

/////////////////////////////////// Partial derivative ///////////////////////////////////
	Uxx = Arrayfloat1(node_add);         Uxz = Arrayfloat1(node_add);
	Uzx = Arrayfloat1(node_add);         Uzz = Arrayfloat1(node_add);
	Uxtt = Arrayfloat1(node_add);        Uztt = Arrayfloat1(node_add);
	Uxx_tmp = Arrayfloat1(node_num);     Uxz_tmp = Arrayfloat1(node_num);
	Uzx_tmp = Arrayfloat1(node_num);     Uzz_tmp = Arrayfloat1(node_num);
	Rix_sem = Arrayfloat2(node_num, Nt); Riz_sem = Arrayfloat2(node_num, Nt);
////////////////////////////////// Output array --totall //////////////////////////////////
	arrayUu = Arrayfloat1(numUu);
////////////////////////////////////////////////////////////////////////////////////////////////////////
	GetRi(Ns, Nr, Nt, node_num, Loadx_s, Loadz_s, Rix_sem, Riz_sem);
//	Get_Ricker(Ni, Nt, dt, Sem_ri);

//////////////////////////////////////////////// Calculation ////////////////////////////////////////
	for(Ni = 1; Ni <= Nt; Ni++)
	{
		Assemble(num, node_num, elem_num, Elem1_count, node1_total, Elem_array, W, NT, Uxt, Uzt, Ux, Uz, lam, mu, Txx, Txz, Tzx, Tzz, 
			Rx, Rz, KKK, Kxx, Kxz, Kzx, Kzz, node_array, Uxx_tmp, Uxz_tmp, Uzx_tmp, Uzz_tmp, node1_array);

		Assemble_pml(num, pnode_num, pelem_num, Elem2_count, Elem_array, node2_total, W, Rx1, Rx2, Rx3, Rz1, Rz2, Rz3, Rxx, Rxz, Rzx, 
			Rzz, lam, mu, Ux, Uz, Uxt, Uzt, PPP, Pxx, Pxz, Pzx, Pzz, NT, Txx1, Txx2, Txx3, Txx4, Txz1, Txz2, Txz3, Txz4, Tzx1, Tzx2, Tzx3, 
			Tzx4, Tzz1, Tzz2, Tzz3, Tzz4);

		Collect_treat(cnode_num, node1_total, node2_total, node3_array, Rx1, Rx2, Rx3, Rz1, Rz2, Rz3, Rx, Rz);

		Ri = Ricker(Ni - 1, dt, F);
		// Ri = GetRicker(Ni-1, dt, F);

		for(i = 0; i < node_num; i++)
		{
			i1 = node1_array[i]-1;
			tmp = Ux[i1];
			// Ux[i1] = 2.0*tmp - Ux1[i] + (Rx[i] + Loadx_s[0][i]*Ri)*dt2*Mass[i1];
			Ux[i1] = 2.0*tmp - Ux1[i] + (Rx[i] + Rix_sem[i][Ni])*dt2*Mass[i1];
			Ux1[i] = tmp;
			Uxtt[i1] = (Ux[i1] - tmp)/dt;

			tmp = Uz[i1];
			// Uz[i1] = 2.0*tmp - Uz1[i] + (Rz[i] + Loadz_s[0][i]*Ri)*dt2*Mass[i1];
			Uz[i1] = 2.0*tmp - Uz1[i] + (Rz[i] + Riz_sem[i][Ni])*dt2*Mass[i1];
			Uz1[i] = tmp;
			Uztt[i1] = (Uz[i1] - tmp)/dt;

			Ux_tmp[i1] = Ux1[i];
			Uz_tmp[i1] = Uz1[i];

			Uxx[i1] = Uxx_tmp[i];
			Uxz[i1] = Uxz_tmp[i];
			Uzx[i1] = Uzx_tmp[i];
			Uzz[i1] = Uzz_tmp[i];
		}

		for(i = 0;i < bnode_num; i++)
		{
			i1 = node4_array[i] - 1;
			i1 = node2_total[i1] - 1;
			Rx1[i1] = 0.0;  Rz1[i1] = 0.0; 
			Rx2[i1] = 0.0;  Rz2[i1] = 0.0;
			Rx3[i1] = 0.0;  Rz3[i1] = 0.0;
			Rxx[i1] = 0.0;  Rxz[i1] = 0.0;
			Rzx[i1] = 0.0;  Rzz[i1] = 0.0;
		}

		for(i = 0; i < pnode_num; i++)
		{
			i1 = node2_array[i] - 1;

			if(node3_total[i1] > 0)
			{
				continue;
			}

			cm = Mass[i1];
			dx = Dampx[i]; dx1 = Dampx1[i];
			dz = Dampz[i]; dz1 = Dampz1[i];

			dxx = dx*dx;  dxz = dx*dz;
			dzx = dz*dx;  dzz = dz*dz;
			//////////////////////////////////Ux////////////////////////////////////
			tt = Rx1[i] * cm + Plxx[i];
			tmp = Upx1[i];
			Upx1[i] = ((2.0 - dt2 * dxx) * Upx1[i] + (dt * dx - 1.0) * Utx1[i] + dt2 * tt) / (1.0 + dt * dx);
			Utx1[i] = tmp;

			tt = Rx2[i] * cm;
			tmp = Upx2[i];
			Upx2[i] = ((2.0 - dt2 * dxz) * Upx2[i] + (0.5 * dt * (dx + dz) - 1.0) * Utx2[i] + dt2 * tt) / (1.0 + 0.5 * dt * (dx + dz));
			Utx2[i] = tmp;

			tt = Rx3[i] * cm + Plxz[i];
			tmp = Upx3[i];
			Upx3[i] = ((2.0 - dt2 * dzz) * Upx3[i] + (dt * dz - 1.0) * Utx3[i] + dt2 * tt) / (1.0 + dt * dz);
			Utx3[i] = tmp;

			Plxx[i] = (1.0 - dt * dx) * Plxx[i] + dt * dx1 * Rxx[i] * cm;
			Plxz[i] = (1.0 - dt * dz) * Plxz[i] + dt * dz1 * Rxz[i] * cm;
			///////////////////////////////Uz//////////////////////////////////////
			tt = Rz1[i] * cm + Plzx[i];
			tmp = Upz1[i];
			Upz1[i] = ((2.0 - dt2 * dxx) * Upz1[i] + (dt * dx - 1.0) * Utz1[i] + dt2 * tt) / (1.0 + dt * dx);
			Utz1[i] = tmp;

			tt = Rz2[i] * cm;
			tmp = Upz2[i];
			Upz2[i] = ((2.0 - dt2 * dzx) * Upz2[i] + (0.5 * dt * (dx + dz) - 1.0) * Utz2[i] + dt2 * tt) / (1.0 + 0.5 * dt * (dx + dz));
			Utz2[i] = tmp;

			tt = Rz3[i] * cm + Plzz[i];
			tmp = Upz3[i];
			Upz3[i] = ((2.0 - dt2 * dzz) * Upz3[i] + (dt * dz - 1.0) * Utz3[i] + dt2 * tt) / (1.0 + dt * dz);
			Utz3[i] = tmp;

			Plzx[i] = (1.0 - dt * dx) * Plzx[i] + dt * dx1 * Rzx[i] * cm;
			Plzz[i] = (1.0 - dt * dz) * Plzz[i] + dt * dz1 * Rzz[i] * cm;
			///////////////////////////////Uz////////////////////////////////////////
			Ux[i1] = Upx1[i] + Upx2[i] + Upx3[i];
			Uz[i1] = Upz1[i] + Upz2[i] + Upz3[i];
			Ux_tmp[i1] = Utx1[i] + Utx2[i] + Utx3[i];
			Uz_tmp[i1] = Utz1[i] + Utz2[i] + Utz3[i];
		}

////////////////////////////////////////
		if(Ni % nsnap == 0)
		{
			// Output1_adj(nm, Ns, Nr, node_add, node_array, Ux, Uz, Uxx, Uxz, Uzx, Uzz);
			nm++;
			// printf("Output1_adj: %d\n", nm);
		}

		if(Ni % numIntervalt == 0)
		{
			nm1++;
			// Output1_bs(nm1, Ns, Nr, node_add, Uxx, Uxz, Uzx, Uzz);
			GetUu_bs(nm1, node_add, Uxx, Uxz, Uzx, Uzz, arrayUu);
		}
	}
	
	for (i = 0; i < numUu; i++)
	{
		arrayUuadj[i] = arrayUu[i];
	}

	//Outputotall_bs_adj(nm1, Ns, Nr, numUu, arrayUu);
	///////////////////////FREE ARRAYS/////////////////////////////////////////
	free(P);      free(W);
	free(Dampx);  free(Dampz);
	free(Dampx1); free(Dampz1);
	free(Rxx);    free(Rzx);
	free(Rxz);    free(Rzz);
	free(Plxx);   free(Plzx);
	free(Plxz);   free(Plzz);
	free(Upx1);   free(Utx1);   free(Rx1);
	free(Upx2);   free(Utx2);   free(Rx2);
	free(Upx3);   free(Utx3);   free(Rx3);
	free(Upz1);   free(Utz1);   free(Rz1);
	free(Upz2);   free(Utz2);   free(Rz2);
	free(Upz3);   free(Utz3);   free(Rz3);

	free(Txx1);   free(Txz1);   free(Tzx1);   free(Tzz1);
	free(Txx2);   free(Txz2);   free(Tzx2);   free(Tzz2);
	free(Txx3);   free(Txz3);   free(Tzx3);   free(Tzz3);
	free(Txx4);   free(Txz4);   free(Tzx4);   free(Tzz4);

	free(Ux);     free(Uxt);    free(Ux1);    free(Rx);
	free(Uz);     free(Uzt);    free(Uz1);    free(Rz);
	
	free(Txx);    free(Txz);    free(arrayUu);
	free(Tzx);    free(Tzz);

	Freefloat2(Nr,Syn_rex);
	Freefloat2(Nr,Syn_rez);
	Freefloat2(Nr,Syn_uxx);
	Freefloat2(Nr,Syn_uzx);
	Freefloat2(Nr,Syn_uxz);
	Freefloat2(Nr,Syn_uzz);
	Freefloat2(node_num, Rix_sem);
	Freefloat2(node_num, Riz_sem);
	// printf("ccc\n"), sleep(1000);
	
	free(Ux_tmp); free(Uz_tmp);
	free(Uxx);    free(Uxz);      free(Uzx);     free(Uzz);      free(Uxtt);
	free(Uztt);   free(Uxx_tmp);  free(Uxz_tmp); free(Uzx_tmp);  free(Uzz_tmp);
}

/////////////////////////////ricker wavelet///////////////////////////////////////////////
double GetRicker(int nT, double Dt, double F)
{
	int i, numFs;
	double Ricker;
	FILE *fp;
	fp = fopen("../Input/Ricker_wavelet.txt", "r");
	for (i = 0; i < nT; i++) {
		numFs = fscanf(fp, "%lf", &Ricker);
	}
	fclose(fp);
	return Ricker;
}

double Ricker(int nT, double Dt, double F)
{
	double T0, T, Ft;
	
	T0 = 1.0 / F + 0.1;
	T = Dt*nT;
	Ft = (1.0 - 2.0*(PI*PI*F*F*(T-T0)*(T-T0)))*exp(-PI*PI*F*F*(T-T0)*(T-T0));
	return Ft;
}
//////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////output////////////////////////////////////////////////
void Output(int nm,long node_num,long *node1_array,double **node_array,double **FDM_array_x,double **FDM_array_z,float *Ux,float *Uz,float **Ufx,float **Ufz,int Nxw,int Nzw,int L,double Dx,double Dz,double Z,int order,int Nr_s)
{
	FILE *fp1,*fp2;
	char Ux1[30]={"../Output/"};
	char Uz1[30]={"../Output/"};
	char Ux11[30]={"/Snapshot/Ux"};
	char Uz11[30]={"/Snapshot/Uz"};
	char add1[20]={".dat"};
	char string[20];
	char string1[20];
	long i,i1,j;
	double x,y;

	nm++;
	sprintf(string,"%d",Nr_s);
	strcat(Ux1,string);
	strcat(Uz1,string);
	strcat(Ux1,Ux11);
	strcat(Uz1,Uz11);
	sprintf(string1,"%d",nm);
	strcat(Ux1,string1);
	strcat(Uz1,string1);
	strcat(Ux1,add1);
	strcat(Uz1,add1);
	if((fp1=fopen(Ux1,"w"))==NULL)
	{
		printf("cannot open Ux.dat\n");
		exit(1);
	}
	if((fp2=fopen(Uz1,"w"))==NULL)
	{
		printf("cannot open Uz.dat\n");
		exit(1);
	}
	
	for(i = 0; i < node_num; i++)
	{	
		i1 = node1_array[i] - 1;

		fprintf(fp1,"%e %e %.10e\n", node_array[i1][0], -1.0 * node_array[i1][1], Ux[i1]);
		fprintf(fp2,"%e %e %.10e\n", node_array[i1][0], -1.0 * node_array[i1][1], Uz[i1]);
	}

	for(i = 0; i < Nxw; i++)
	{		
		for(j = 0; j < Nzw; j++)
		{
			if(FDM_array_z[i][j] > Z)
			{
				x = FDM_array_x[i][j];
				y = -1.0 * FDM_array_z[i][j];

				fprintf(fp1,"%e  %e  %.10e\n", x, y, Ufx[i][j]);
				fprintf(fp2,"%e  %e  %.10e\n", x, y, Ufz[i][j]);
			}
		}
	}

	fclose(fp1);
	fclose(fp2);
}

void Output1_syn(int nm, int Ns, int Nr, long node_num, double **node_array, float *Ux, float *Uz,
	float *Uxx, float *Uxz, float *Uzx, float *Uzz)
{
	FILE *fp1,*fp2;
	long i;
	double xx, zz;
	char flnm1[100], flnm2[100];

	nm++;
	sprintf(flnm1,"../Output/Snapshot_syn/wavefiled_s%d_r%d_Ux%d.txt", Ns, Nr, nm);
	sprintf(flnm2,"../Output/Snapshot_syn/wavefiled_s%d_r%d_Uz%d.txt", Ns, Nr, nm);

	if((fp1 = fopen(flnm1, "w")) == NULL)
	{
		printf("cannot open %s\n", flnm1), exit(1);
	}
	if((fp2 = fopen(flnm2, "w")) == NULL)
	{
		printf("cannot open %s\n", flnm2), exit(0);
	}

	for(i = 0; i < node_num; i++)
	{
		xx = node_array[i][0];
		zz = node_array[i][1];
		fprintf(fp1, "%.3f  %.3f  %.10e %.10e %.10e\n", xx, -1.0*zz, Ux[i], Uxx[i], Uxz[i]), fflush(fp1);
		fprintf(fp2, "%.3f  %.3f  %.10e %.10e %.10e\n", xx, -1.0*zz, Uz[i], Uzx[i], Uzz[i]), fflush(fp2);
	}

	fclose(fp1); fclose(fp2);
}

void Output1_adj(int nm, int Ns, int Nr, long node_num, double **node_array, float *Ux, float *Uz,
	float *Uxx, float *Uxz, float *Uzx, float *Uzz)
{
	FILE *fp1,*fp2;
	long i;
	double xx, zz;
	char flnm1[100], flnm2[100];

	nm++;
	sprintf(flnm1,"../Output/Snapshot_adj/wavefiled_s%d_r%d_Ux%d.txt", Ns, Nr, nm);
	sprintf(flnm2,"../Output/Snapshot_adj/wavefiled_s%d_r%d_Uz%d.txt", Ns, Nr, nm);

	if((fp1 = fopen(flnm1, "w")) == NULL) {
		printf("cannot open %s\n", flnm1), exit(1);
	}
	if((fp2 = fopen(flnm2, "w")) == NULL) {
		printf("cannot open %s\n", flnm2), exit(0);
	}

	for(i = 0; i < node_num; i++) {
		xx = node_array[i][0];
		zz = node_array[i][1];
		fprintf(fp1, "%.3f  %.3f  %.10e %.10e %.10e\n", xx, -1.0*zz, Ux[i], Uxx[i], Uxz[i]), fflush(fp1);
		fprintf(fp2, "%.3f  %.3f  %.10e %.10e %.10e\n", xx, -1.0*zz, Uz[i], Uzx[i], Uzz[i]), fflush(fp2);
	}

	fclose(fp1);
	fclose(fp2);
}

void Output1_bs(int nm1, int Ns, int Nr, long node_num, float *Uxx, float *Uxz, float *Uzx, float *Uzz)
{
	long i, j, i1, numUu, id, sizef;
	float arr_f1[4], *arrayUu;
	char file1[100];
	FILE *fp1;

	sprintf(file1,"../Output/Snapshot_bianry_syn/source%d/r%d_Uu%d.bin", Ns, Nr, nm1);

	////////////////// 2 //////////////////
	numUu = 4*node_num;
	sizef = sizeof(float);
	arrayUu = Arrayfloat1(numUu);

	id = 0;
	for(i = 0; i < node_num; i++)
	{	
		arrayUu[id] = Uxx[i];
		id++;
		arrayUu[id] = Uxz[i];
		id++;
		arrayUu[id] = Uzx[i];
		id++;
		arrayUu[id] = Uzz[i];
		id++;
	}

	if((fp1 = fopen(file1, "wb")) == NULL)
	{
		printf("Can not open U_p.dat -- Output1_bs\n"); exit(1);
	}

	fwrite(arrayUu, sizef, numUu, fp1);
	fclose(fp1);

	free(arrayUu);
	///////////////////////////////////////
}

void GetUu_bs(int nm1, long node_num, float *Uxx, float *Uxz, float *Uzx, float *Uzz, float *arrayUu)
{
	long i, id;

	////////////////// 2 //////////////////
	id = (nm1 - 1)*4*node_num;
	for(i = 0; i < node_num; i++)
	{	
		arrayUu[id] = Uxx[i];
		id++;
		arrayUu[id] = Uxz[i];
		id++;
		arrayUu[id] = Uzx[i];
		id++;
		arrayUu[id] = Uzz[i];
		id++;
	}
	///////////////////////////////////////
}

void Outputotall_bs_syn(int nm1, int Ns, int Nr, long numUu, float *arrayUu)
{
	int sizei, sizel, sizef;
	long numNode;
	char file1[100];
	FILE *fp1;

	/// File name ///
	sprintf(file1,"../Output/Snapshot_bianry_syn/Uu_s%d_r%d.bin", Ns, Nr);

	numNode = numUu/(4*nm1);
	/// Output file ///
	sizei = sizeof(int);
	sizel = sizeof(long);
	sizef = sizeof(float);

	if((fp1 = fopen(file1, "wb")) == NULL)
	{
		printf("Can not open U_p.dat -- Output1_bs\n"); exit(1);
	}

	fwrite(&nm1, sizei, 1, fp1);
	fwrite(&numNode, sizel, 1, fp1);
	fwrite(arrayUu, sizef, numUu, fp1);
	fclose(fp1);
	///////////////////////////////////////
}

void Outputotall_bs_adj(int nm1, int Ns, int Nr, long numUu, float *arrayUu)
{
	int sizei, sizel, sizef;
	long numNode;
	char file1[100];
	FILE *fp1;

	/// File name ///
	sprintf(file1,"../Output/Snapshot_bianry_adj/Uu_s%d_r%d.bin", Ns, Nr);

	numNode = numUu/(4*nm1);
	/// Output file ///
	sizei = sizeof(int);
	sizel = sizeof(long);
	sizef = sizeof(float);

	if((fp1 = fopen(file1, "wb")) == NULL)
	{
		printf("Can not open U_p.dat -- Output1_bs\n"); exit(1);
	}

	fwrite(&nm1, sizei, 1, fp1);
	fwrite(&numNode, sizel, 1, fp1);
	fwrite(arrayUu, sizef, numUu, fp1);
	fclose(fp1);
	///////////////////////////////////////
}

void Output1_diff(int nm1,long node_num,double **node_array,double **FDM_array_x,double **FDM_array_z,float *Uxx,float *Uxz,float *Uzx,float *Uzz,float *Uxtt,float *Uztt,float **Ufxx,float **Ufxz,float **Ufzx,float **Ufzz,float **Ufxtt,float **Ufztt,int Nxw,int Nzw,int L,double len,double Dx,double Dz,double Z,int order,int Nr_s)

{
	FILE *fp1,*fp2,*fp3,*fp4,*fp5,*fp6;
	char Ux1[30]={"../Output/"};
	char Ux2[30]={"../Output/"};
	char Uz1[30]={"../Output/"};
	char Uz2[30]={"../Output/"};
	char Uxt[30]={"../Output/"};
	char Uzt[30]={"../Output/"};
	
	char Ux11[30]={"/Diff/Uxx"};
	//char Ux11[30]={"/Diff/Uxx/Uxx"};
	char Ux21[30]={"/Diff/Uxz/Uxz"};
	char Uz11[30]={"/Diff/Uzx/Uzx"};
	char Uz21[30]={"/Diff/Uzz/Uzz"};
	char Uxt1[30]={"/Diff/Uxt/Uxt"};
	char Uzt1[30]={"/Diff/Uzt/Uzt"};

	char add1[20]={".dat"};
	char string[20];
	char string1[20];

	long i,j;
	double x,y;
	nm1++;
	Nr_s++;
	sprintf(string,"%d",Nr_s);
	strcat(Ux1,string);
	strcat(Uz1,string);
	strcat(Ux2,string);
	strcat(Uz2,string);
	strcat(Uxt,string);
	strcat(Uzt,string);
	strcat(Ux1,Ux11);
	strcat(Uz1,Uz11);
	strcat(Ux2,Ux21);
	strcat(Uz2,Uz21);
	strcat(Uxt,Uxt1);
	strcat(Uzt,Uzt1);
	sprintf(string1,"%d",nm1);
	strcat(Ux1,string1);
	strcat(Uz1,string1);
	strcat(Ux2,string1);
	strcat(Uz2,string1);
	strcat(Uxt,string1);
	strcat(Uzt,string1);
	strcat(Ux1,add1);
	strcat(Uz1,add1);
	strcat(Ux2,add1);
	strcat(Uz2,add1);
	strcat(Uxt,add1);
	strcat(Uzt,add1);

	if((fp1=fopen(Ux1,"w"))==NULL)
	{
		printf("cannot open Uxx.dat\n");
		exit(0);
	}

	for(i=0;i<node_num;i++)
	{
		fprintf(fp1,"%e  %e  %.10e %.10e %.10e %.10e %.10e %.10e\n",node_array[i][0],-1.0*node_array[i][1],Uxx[i],Uzx[i],Uxz[i],Uzz[i],Uxtt[i],Uztt[i]);
	}

	for(i=L;i<=Nxw-L;i++)
	{		
		for(j=0;j<Nzw;j++)
		if(FDM_array_z[i][j]>Z)
		{
			x=FDM_array_x[i][j];
			y=-1.0*FDM_array_z[i][j];
			fprintf(fp1,"%e  %e  %.10e %.10e %.10e %.10e %.10e %.10e\n",x,y,Ufxx[i][j],Ufzx[i][j],Ufxz[i][j],Ufzz[i][j],Ufxtt[i][j],Ufztt[i][j]);
		}
	}

	fclose(fp1);
}


void Damp_coff(short num,long pnode_num,double X,double Z,double R0,double Vp,double Layx,double Layz,long *node2_array,double **node_array,
	float *Dampx,float *Dampz,float *Dampx1,float *Dampz1)
{
	long i,id;
	double Lx,Lz,xx,zz;

	for(i=0;i<pnode_num;i++)
	{
		id=node2_array[i]-1;
		xx=node_array[id][0];
		zz=node_array[id][1];

		///////////////face1///////////////
		if(xx<=0.0)
		{
			Lx=xx;
			Lz=0.0;
		}
		///////////////face2//////////////
		if(xx>=X)
		{
			Lx=xx-X;
			Lz=0.0;
		}
		///////////////face3//////////////
		if(zz>=Z)
		{
			Lz=zz-Z;
			if(xx<=0.0)
			{
				Lx=xx;
			}
			else if(xx>=X)
			{
				Lx=xx-X;
			}
			else
			{
				Lx=0.0;
			}
		}

		Dampx[i] = 3.0 * Vp * log10(1.0 / R0) / (2.0 * Layx) * (Lx / Layx) * (Lx / Layx);
		Dampz[i] = 3.0 * Vp * log10(1.0 / R0) / (2.0 * Layz) * (Lz / Layz) * (Lz / Layz);	
		Dampx1[i] = 6.0 * Vp * log10(1.0 / R0) / (2.0 * Layx) * (Lx / Layx) / Layx;
		Dampz1[i] = 6.0 * Vp * log10(1.0 / R0) / (2.0 * Layz) * (Lz / Layz) / Layz;

		// printf("%20.4e%20.4e\n", Dampx[i], Dampx1[i]), sleep(1000);

		// Dampx[i]=3.0*Vp*log10(1.0/R0)/(2.0*Layx)*(Lx/Layx)*(Lx/Layx);
		// Dampz[i]=3.0*Vp*log10(1.0/R0)/(2.0*1)*(Lz/1)*(Lz/1);	
		// Dampx1[i]=6.0*Vp*log10(1.0/R0)/(2.0*Layx)*(Lx/Layx)/Layx;
		// Dampz1[i]=6.0*Vp*log10(1.0/R0)/(2.0*1)*(Lz/1)/1;
	}
}

void Get_record_fdm_pml(short num,int Nx,int L,long zf,double **Fii_fdm,long *Elem2_count,long **Elem_array,double *P,
	double **Fax_fdm,double **Faz_fdm,float *Ux_tmp,float *Uz_tmp,float **Ufx_tmp,float **Ufz_tmp,double **node_array)
{
	long id,id1;
	short i,i1,i2,j,ix,iz;	
	double cfx,cfz,valx,valz;
	double tmp11,tmp21;
	
	for(i = 0; i < Nx + 2*L + 1; i++)
	{
		if(i < L)
		{
			for(j = 0; j <= zf; j++)
			{
				valx = 0.0;  valz = 0.0;

				id = Fii_fdm[i][j];
				id = Elem2_count[id] - 1;

				for(i1 = 1; i1 < num; i1++)
				{
					tmp11=1.0;  tmp21=1.0;

					for(ix = 1; ix < num; ix++)
					{
						if(ix == i1)continue;

						tmp11 *= (P[i1] - P[ix]);
						tmp21 *= (Faz_fdm[i][j] - P[ix]);
					}

					cfx = tmp21 / tmp11;

					for(i2 = 0; i2 < num; i2++)
					{
						cfz = Lagrange_val(num, i2, P, Fax_fdm[i][j]);

						id1 = Elem_array[id][i1*num + i2] - 1;

						valx += cfx*cfz*Ux_tmp[id1];
						valz += cfx*cfz*Uz_tmp[id1];
					}
				}

				Ufx_tmp[i][j] = valx;
				Ufz_tmp[i][j] = valz;
			}
		}

		if(i > Nx + L)
		{
			for(j = 0; j <= zf; j++)
			{
				valx = 0.0;  valz = 0.0;

				id = Fii_fdm[i][j];
				id = Elem2_count[id] - 1;

				for(i1 = 0; i1 < num; i1++)
				{
					cfx = Lagrange_val(num, i1, P, Fax_fdm[i][j]);

					for(i2 = 0; i2 < num - 1; i2++)
					{
						cfz = Lagrange_val(num-1, i2, P, Faz_fdm[i][j]);

						id1 = Elem_array[id][i1*num + i2] - 1;

						valx += cfx*cfz*Ux_tmp[id1];
						valz += cfx*cfz*Uz_tmp[id1];

					}
				}

				Ufx_tmp[i][j] = valx;
				Ufz_tmp[i][j] = valz;

			}
		}


	}
}

void Output_record(long Nr,long Ns,long Nt,double dt,float **Rx,float **Rz,int Nr_s)
{
	long i,j;
	char flnm[100];
	FILE *fp;

	sprintf(flnm, "../Output/Receivers_syn/source%ld/D%ld", Ns, Nr);
	if((fp=fopen(flnm, "w"))==NULL)
	{
		printf("cannot open %s\n", flnm), exit(1);
	}

	for(j=0;j<Nt;j++)
	{
		fprintf(fp,"%e %.10e %.10e\n",(j+1.0)*dt,Rx[0][j],Rz[0][j]), fflush(fp);
	}
	fclose(fp);
}

void Output_record_syn(long Nr,long Ns,long Nt,double dt,float **Rx,float **Rz,int Nr_s)
{
	long i,j;
	char flnm[100];
	FILE *fp;

	sprintf(flnm, "../Output/Receivers_syn/wave_s%ld_r%ld.txt", Ns, Nr);
	if((fp=fopen(flnm, "w"))==NULL) {
		printf("cannot open %s\n", flnm), exit(1);
	}

	for(j=0;j<Nt;j++) {
		fprintf(fp,"%e %.10e %.10e\n",(j+1.0)*dt,Rx[0][j],Rz[0][j]), fflush(fp);
	}
	fclose(fp);
}

void Find_sf_position_rec(short num, short Nr, int Nxw, int Nzw, double Z, float **Ufx, float **Ufz,
	 double **Rec_f, double **FDM_array_x, double **FDM_array_z, double **Rloc)
{
	int i, j, k;
	double fx1, fx2, fz1, fz2, rx, rz, rr, dx;

	dx = FDM_array_x[1][0] - FDM_array_x[0][0];

	for(i = 0; i < Nr; i++)
	{
		rx = Rloc[i][0];
		rz = Rloc[i][1];

		for (j = 0; j < Nxw - 1; j++)
		{
			for (k = 0; k < Nzw - 1; k++)
			{
				fx1 = FDM_array_x[j][k];
				fz1 = FDM_array_z[j][k];

				fx2 = FDM_array_x[j + 1][k + 1];
				fz2 = FDM_array_z[j + 1][k + 1];

				if (rx >= fx1 && rx < fx2 && rz >= fz1 && rz < fz2)
				{
					Rec_f[i][0] = j;
					Rec_f[i][1] = k;
					Rec_f[i][2] = (rx - fx2) / -dx;
					Rec_f[i][3] = (rx - fx1) / dx;
					Rec_f[i][4] = (rz - fz2) / -dx;
					Rec_f[i][5] = (rz - fz1) / dx;


					// printf("%20.4e%20.4e\n", Rec_f[i][2], Rec_f[i][3]);
					// printf("%20.4e%20.4e%20.4e\n", rx, fx1, fx2), sleep(1000);
				}
			}
		}
	}
}

void Get_record(short num,short Nr,long Ni,long *Rii,long *Elem1_count,long **Elem_array,double *P,double *Rax,double *Raz,float *Ux,float *Uz,float **Syn_rex,float **Syn_rez)
{
	long id,id1;
	short i,i1,i2;	
	double cfx,cfz,valx,valz;
	
	for(i=0;i<Nr;i++)
	{
		valx=0.0;
		valz=0.0;
		id=Rii[i];
		id=Elem1_count[id]-1;

		for(i1=0;i1<num;i1++)
		{
			cfx=Lagrange_val(num,i1,P,Rax[i]);
			for(i2=0;i2<num;i2++)
			{
				cfz=Lagrange_val(num,i2,P,Raz[i]);
				id1=Elem_array[id][i1*num+i2]-1;
				valx+=cfx*cfz*Ux[id1];
				valz+=cfx*cfz*Uz[id1];
			}
		}
		Syn_rex[i][Ni-1]=valx;
		Syn_rez[i][Ni-1]=valz;
	}
}


void Get_record_rec(short num, short Nr, long Ni, double Z, long *Rii, float *Ux, float *Uz, long *Elem1_count, double *P, double *Rax, double *Raz,
	long **Elem_array, float **Ufx, float **Ufz, float **Syn_rex,float **Syn_rez, double **Rec_f)
{
	short i, i1, i2, j, k;

	long id, id1;

	double cfx[2], cfz[2], fx[2][2], fz[2][2], coff, valx, valz;

	for(i = 0; i < Nr; i++)
	{
		valx = 0.0;
		valz = 0.0;
		id = Rii[i];

		if(id != 0)
		{
			id = Elem1_count[id] - 1;

			for(i1 = 0; i1 < num; i1++)
			{
				cfx[0] = Lagrange_val(num, i1, P, Rax[i]);

				for(i2 = 0; i2 < num; i2++)
				{
					cfz[0] = Lagrange_val(num, i2, P, Raz[i]);

					coff = cfx[0]*cfz[0];

					id1 = Elem_array[id][i1*num + i2] - 1;
					valx += coff*Ux[id1];
					valz += coff*Uz[id1];
				}
			}
		}
		else
		{
			i1 = round(Rec_f[i][0]);   i2 = round(Rec_f[i][1]);

			cfx[0] = Rec_f[i][2];      cfx[1] = Rec_f[i][3];
			cfz[0] = Rec_f[i][4];      cfz[1] = Rec_f[i][5];

			fx[0][0] = Ufx[i1][i2];    fx[1][0] = Ufx[i1 + 1][i2];
			fx[0][1] = Ufx[i1][i2+ 1]; fx[1][1] = Ufx[i1 + 1][i2 + 1];

			fz[0][0] = Ufz[i1][i2];    fz[1][0] = Ufz[i1 + 1][i2];
			fz[0][1] = Ufz[i1][i2+ 1]; fz[1][1] = Ufz[i1 + 1][i2 + 1];

			for (j = 0; j < 2; j++)
			{
				for (k = 0; k < 2; k++)
				{
					coff = cfx[j]*cfz[k];

					valx += coff*fx[j][k];
					valz += coff*fz[j][k];
				}
			}
		}

		Syn_rex[i][Ni - 1] = valx;
		Syn_rez[i][Ni - 1] = valz;
	}
}

void Find_local_position_fdm_num(short num, short ss, int Nx, int L, long zf, double Z, long elem_num, double len, double **FDM_array_x,
	double **FDM_array_z, long *Elem1_count, long **Elem_array, double **node_array, double **Fii_fdm, double **Fax_fdm, double **Faz_fdm)
{
	short N3, nn1, flag, Cin[4];
	int i, j, k, it, i0, i1, i2, i3, id, flgi;
	double xx, yy, zz, ax, ay, az, a1, a2, a3, a4, a5, a6, a7, a8, tmp0, tmp, Dis,  P1[400];
	double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, px, pz, flgdx, flgdz, semx1, semx2, semz1, semz2;

	Dis = 2.0*len;
	nn1 = 20*20-1;

	Cin[0] = 0;
	Cin[1] = num*(num - 1);
	Cin[2] = num*num - 1;
	Cin[3] = num - 1;

	for(i = 0; i < nn1; i++)
	{
		P1[i] = -1.0 + 2.0 / (nn1 - 1) * i;
	}

	for(i = L; i <= Nx + L + 1; i++)
	{		
		for(j = 0; j <= zf + 1; j++)
		{
			flgi = 0;

			tmp0 = Dis;

			for(k = 0; k < elem_num; k++)
			{
				it = Elem1_count[k]-1;

				i1 = Elem_array[it][Cin[0]]-1;
				x1 = node_array[i1][0];
				z1 = node_array[i1][1];

				i1 = Elem_array[it][Cin[1]]-1;
				x2 = node_array[i1][0];
				z2 = node_array[i1][1];

				i1 = Elem_array[it][Cin[2]]-1;
				x3 = node_array[i1][0];
				z3 = node_array[i1][1];

				i1 = Elem_array[it][Cin[3]]-1;
				x4 = node_array[i1][0];
				z4 = node_array[i1][1];

				semx1 = (x1 + x2 + x3 + x4) / 4;
				semz1 = (z1 + z2 + z3 + z4) / 4;

				xx = semx1 - FDM_array_x[i][j];
				zz = semz1 - FDM_array_z[i][j];

				if(semz1 >= FDM_array_z[0][0])
				{
					tmp = sqrt(xx*xx + zz*zz);

					if(tmp < Dis)
					{
						for(i1 = 0; i1 < nn1; i1++)
						{
							px = P1[i1];

							for(i2 = 0; i2 < nn1; i2++)
							{
								pz = P1[i2];

								a1 = (1.0 - px)*(1.0 - pz) / 4;
								a2 = (1.0 + px)*(1.0 - pz) / 4;
								a3 = (1.0 + px)*(1.0 + pz) / 4;
								a4 = (1.0 - px)*(1.0 + pz) / 4;

								semx2 = a1*x1 + a2*x2 + a3*x3 + a4*x4;
								semz2 = a1*z1 + a2*z2 + a3*z3 + a4*z4;

								xx = semx2 - FDM_array_x[i][j];
								zz = semz2 - FDM_array_z[i][j];

								tmp = sqrt(xx*xx + zz*zz);

								if(tmp < tmp0)
								{
									Fii_fdm[i][j] = k;
									Fax_fdm[i][j] = px;
									Faz_fdm[i][j] = pz;

									flgi = 1;

									tmp0 = tmp;
								}
							}
						}
					}
				}
			}

			if (flgi == 0)
			{
				Fii_fdm[i][j] = -1;
			}
		}
	}
}

//////////////////////////////////////////////// Find nodes from sem to fdm in pml layer //////////////////////////////////////////////////
void Find_local_position_fdm_pml(short num, int Nx, int L, long zf, double Z, long pelem_num, long pnode_num, double len,
	long *node2_array, long *Elem2_count,long **Elem_array, double **FDM_array_x, double **FDM_array_z, double **node_array,
	double **Fax_fdm, double **Faz_fdm , double ****Fid_s2f_pml)
{
	short N3, nn1, flag, Cin[4];
	int i, j, k, l, it, i1, i2, i3, i4, id, id1, id2, NN, flgi1, flgi2;
	double xx1, xx2, zz1, zz2, ax, ay, az, a1, a2, a3, a4, tmp1, tmp2, tmp3, Dis;
	double x1, x2, x3, x4, z1, z2, z3, z4, px, py, pz;
	double semx1, semx2, semz1, semz2, fdmx1, fdmz1, fdmx2, fdmz2, flgd1, P1[400], **arr_e2n;

	Dis = 2.0*len;

	nn1 = 20*20 - 1;

	NN = num*num;

	arr_e2n = Arraydouble2(pelem_num, NN);

	fdmx1 = FDM_array_x[Nx + L][0];
	fdmz1 = FDM_array_z[0][zf];

	Cin[0] = 0;
	Cin[1] = num*(num - 1);
	Cin[2] = num*num - 1;
	Cin[3] = num - 1;

	for(i = 0; i < nn1; i++)
	{
		P1[i] = -1.0 + 2.0 / (nn1 - 1)*i;
	}

	for(i = 0; i < pelem_num; i++)
	{
		i1 = Elem2_count[i] - 1;

		for(j = 0; j < NN; j++)
		{
			i2 = j;
			id1 = Elem_array[i1][i2] - 1;


			// int ii, jj;
			// for (ii = 0; ii < num; ii++)
			// {
			// 	for (jj = 0; jj < num; jj++)
			// 	{
			// 		id = Elem_array[i1][ii*num + jj]-1;

			// 		ax = node_array[id][0];
			// 		az = node_array[id][1];

			// 		printf("%20.4e%20.4e\n", ax, az);
			// 	}
			// 	printf("\n\n");
			// }


			// printf("cc\n"), sleep(1000);

			flgi1 = 0;

			for (k = 0; k < pnode_num; k++)
			{
				id2 = node2_array[k] - 1;

				if (id1 == id2)
				{
					arr_e2n[i][i2] = k;

					flgi1 = 1;
					break;
				}				
			}

			if (flgi1 == 0)
			{
				printf("Sonmething error! cannot find pml nodes!\n"), exit(1);
			}
		}
	}

	for(i = 0; i < Nx + 2*L + 1; i++)
	{		
		if(i < L || i > Nx + L)
		{
			for(j = 0; j <= zf; j++)
			{ 
				fdmx2 = FDM_array_x[i][j];
				fdmz2 = FDM_array_z[i][j];

				flgi2 = 0;
				tmp2 = Dis;

				for(k = 0; k < pelem_num; k++)
				{
					it = Elem2_count[k]-1;

					i1 = Elem_array[it][Cin[0]] - 1;
					x1 = node_array[i1][0];
					z1 = node_array[i1][1];

					i1 = Elem_array[it][Cin[1]] - 1;
					x2 = node_array[i1][0];
					z2 = node_array[i1][1];

					i1 = Elem_array[it][Cin[2]] - 1;
					x3 = node_array[i1][0];
					z3 = node_array[i1][1];

					i1 = Elem_array[it][Cin[3]] - 1;
					x4 = node_array[i1][0];
					z4 = node_array[i1][1];

					semx1 = (x1 + x2 + x3 + x4) / 4;
					semz1 = (z1 + z2 + z3 + z4) / 4;

//					if(semz1 >= FDM_array_z[0][0] && semz1 < Z)
					if(semz1 >= FDM_array_z[0][0])
					{
						xx1 = semx1 - fdmx2;
						zz1 = semz1 - fdmz2;

						tmp1 = sqrt(xx1*xx1 + zz1*zz1);

						if(tmp1 < Dis)
						{
							for(i1 = 0; i1 < nn1; i1++)
							{
								px = P1[i1];

								for(i2 = 0; i2 < nn1; i2++)
								{
									pz = P1[i2];

									a1 = (1.0 - px)*(1.0 - pz) / 4;
									a2 = (1.0 + px)*(1.0 - pz) / 4;
									a3 = (1.0 + px)*(1.0 + pz) / 4;
									a4 = (1.0 - px)*(1.0 + pz) / 4;

									semx2 = a1*x1 + a2*x2 + a3*x3 + a4*x4;
									semz2 = a1*z1 + a2*z2 + a3*z3 + a4*z4;

									xx2 = semx2 - fdmx2;
									zz2 = semz2 - fdmz2;

									tmp3 = sqrt(xx2*xx2 + zz2*zz2);

									if (tmp3 < tmp2)
									{
										flgi2 = 1;

										for (i3 = 0; i3 < NN; i3++)
										{
											Fid_s2f_pml[i][j][k][i3] = arr_e2n[k][i3];
										}

										Fid_s2f_pml[i][j][pelem_num][NN] = k;
										Fid_s2f_pml[i][j][k][NN + 1] = px;
										Fid_s2f_pml[i][j][k][NN + 2] = pz;

										tmp2 = tmp3;
									}
								}
							}
						}
					}
				}

				if (flgi2 == 0)
				{
					Fid_s2f_pml[i][j][pelem_num][NN] = -1;
					// printf("no! %20.4f%20.4f\n", fdmx2, fdmz2);
				}
			}
		}
	}

	Freedouble2(pelem_num, arr_e2n);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Get_record_fdm_num(short num, int Nx, int L, long zf, long *Elem1_count, float *Ux_tmp, float *Uz_tmp, double *P,
	double **Fii_fdm, long **Elem_array, double **Fax_fdm, double **Faz_fdm, float **Ufx_tmp, float **Ufz_tmp)
{
	int i, i1, i2, j;
	long id, id1;
	double cfx, cfz, valx, valz;

	// zf = zf + 4;

	for(i = L; i <= Nx + L + 1; i++)
	{
		for(j = 0; j <= zf + 1; j++)
		{
			valx = 0.0;  valz = 0.0;

			if (Fii_fdm[i][j] >= 0)
			{
				id = Fii_fdm[i][j];
				id = Elem1_count[id] - 1;

				for(i1 = 0 ; i1 < num; i1++)
				{
					cfx = Lagrange_val(num, i1, P, Fax_fdm[i][j]);

					for(i2 = 0; i2 < num; i2++)
					{
						cfz = Lagrange_val(num, i2, P, Faz_fdm[i][j]);

						id1 = Elem_array[id][i1*num + i2] - 1;

						valx += cfx*cfz*Ux_tmp[id1];

						valz += cfx*cfz*Uz_tmp[id1];
					}
				}

				Ufx_tmp[i][j] = valx;
				Ufz_tmp[i][j] = valz;
			}
		}
	}
}

void Get_record_fdm_pml_c(short num, int Nx, int L, long zf, long pelem_num, long *Elem2_count, float *Upx1, float *Upx2, float *Upx3, float *Upz1, float *Upz2,
	float *Upz3, float *Ux_tmp, float *Uz_tmp,  double *P, long **Elem_array, float **Ufx_tmp, float **Ufz_tmp, double **Fax_fdm, double **Faz_fdm,
	double **node_array, float ****Ulx, float ****Ulz, float ****Urx, float ****Urz, double ****Fid_s2f_pml)
{
	int i, i1, i2, j, ix, iz, id1, id2, id3, id4, id5, id6, num1, NN, flgi1;
	double cfx1, cfx2, cfz1, cfz2, valx1, valx2, valx3, valz1, valz2, valz3, vvx, vvz, tmp11, tmp21;
	double plx1, plx2, plx3, plz1, plz2, plz3, semx, semz;

	NN = num*num;

	// zf = zf + 4;

	for(i = 0; i < Nx + 2*L + 1; i++)
	{
		if(i < L)
		{
			for(j = 0; j < zf + 1; j++)
			{
				valx1 = 0.0; valx2 = 0.0; valx3 = 0.0;
				valz1 = 0.0; valz2 = 0.0; valz3 = 0.0;

				vvx = 0.0;   vvz = 0.0;

				id1 = Fid_s2f_pml[i][j][pelem_num][NN];

				if (id1 >= 0)
				{
					cfx1 = Fid_s2f_pml[i][j][id1][NN + 1];
					cfz1 = Fid_s2f_pml[i][j][id1][NN + 2];

					id5 = Elem2_count[id1] - 1;

					for(i1 = 0; i1 < num; i1++)
					{
						cfx2 = Lagrange_val(num, i1, P, cfx1);

						for(i2 = 0; i2 < num; i2++)
						{
							cfz2 = Lagrange_val(num, i2, P, cfz1);

							id2 = i1*num + i2;
							id3 = Fid_s2f_pml[i][j][id1][id2];
							id6 = Elem_array[id5][id2] - 1;

							plx1 = Upx1[id3];          plz1 = Upz1[id3];
							plx2 = Upx2[id3];          plz2 = Upz2[id3];
							plx3 = Upx3[id3];          plz3 = Upz3[id3];

							valx1 += cfx2*cfz2*plx1;   valz1 += cfx2*cfz2*plz1;
							valx2 += cfx2*cfz2*plx2;   valz2 += cfx2*cfz2*plz2;
							valx3 += cfx2*cfz2*plx3;   valz3 += cfx2*cfz2*plz3;

							vvx += cfx2*cfz2*Ux_tmp[id6];
							vvz += cfx2*cfz2*Uz_tmp[id6];
						}
					}

					id4 = i;

					Ulx[1][0][id4][j] = valx1;         Ulz[1][0][id4][j] = valz1;
					Ulx[1][1][id4][j] = valx3;         Ulz[1][1][id4][j] = valz3;
					Ulx[1][2][id4][j] = valx2;         Ulz[1][2][id4][j] = valz2;

					Ufx_tmp[i][j] = Ulx[1][0][id4][j] + Ulx[1][1][id4][j] + Ulx[1][2][id4][j];
					Ufz_tmp[i][j] = Ulz[1][0][id4][j] + Ulz[1][1][id4][j] + Ulz[1][2][id4][j];

					// if (i == L)
					// {
					// 	printf("%20.4e%20.4e\n", Ufx_tmp[i][j], Ufz_tmp[i][j]);
					// }

				}
			}
		}
		else if(i > Nx + L)
		{
			for(j = 0; j <= zf; j++)
			{
				valx1 = 0.0; valx2 = 0.0; valx3 = 0.0;
				valz1 = 0.0; valz2 = 0.0; valz3 = 0.0;

				vvx = 0.0;   vvz = 0.0;

				id1 = Fid_s2f_pml[i][j][pelem_num][NN];

				if (id1 >= 0)
				{
					cfx1 = Fid_s2f_pml[i][j][id1][NN + 1];
					cfz1 = Fid_s2f_pml[i][j][id1][NN + 2];

					id5 = Elem2_count[id1] - 1;
					num1 = num - 1;


					if (i == Nx + L + 1)
					{
						for(i1 = 0; i1 < num1; i1++)
						{
							cfx2 = Lagrange_val(num1, i1, P, cfx1);

							for(i2 = 0; i2 < num; i2++)
							{
								cfz2 = Lagrange_val(num, i2, P, cfz1);

								id2 = i1*num + i2;
								id3 = Fid_s2f_pml[i][j][id1][id2];
								id6 = Elem_array[id5][id2] - 1;

								plx1 = Upx1[id3];          plz1 = Upz1[id3];
								plx2 = Upx2[id3];          plz2 = Upz2[id3];
								plx3 = Upx3[id3];          plz3 = Upz3[id3];

								valx1 += cfx2*cfz2*plx1;   valz1 += cfx2*cfz2*plz1;
								valx2 += cfx2*cfz2*plx2;   valz2 += cfx2*cfz2*plz2;
								valx3 += cfx2*cfz2*plx3;   valz3 += cfx2*cfz2*plz3;

								vvx += cfx2*cfz2*Ux_tmp[id6];
								vvz += cfx2*cfz2*Uz_tmp[id6];
							}
						}
					}
					else
					{
						for(i1 = 0; i1 < num; i1++)
						{
							cfx2 = Lagrange_val(num, i1, P, cfx1);

							for(i2 = 0; i2 < num; i2++)
							{
								cfz2 = Lagrange_val(num, i2, P, cfz1);

								id2 = i1*num + i2;
								id3 = Fid_s2f_pml[i][j][id1][id2];
								id6 = Elem_array[id5][id2] - 1;

								plx1 = Upx1[id3];          plz1 = Upz1[id3];
								plx2 = Upx2[id3];          plz2 = Upz2[id3];
								plx3 = Upx3[id3];          plz3 = Upz3[id3];

								valx1 += cfx2*cfz2*plx1;   valz1 += cfx2*cfz2*plz1;
								valx2 += cfx2*cfz2*plx2;   valz2 += cfx2*cfz2*plz2;
								valx3 += cfx2*cfz2*plx3;   valz3 += cfx2*cfz2*plz3;

								vvx += cfx2*cfz2*Ux_tmp[id6];
								vvz += cfx2*cfz2*Uz_tmp[id6];
							}
						}
					}


					id4 = i - (Nx + L + 1);

					Urx[1][0][id4][j] = valx1;         Urz[1][0][id4][j] = valz1;
					Urx[1][1][id4][j] = valx3;         Urz[1][1][id4][j] = valz3;
					Urx[1][2][id4][j] = valx2;         Urz[1][2][id4][j] = valz2;

					Ufx_tmp[i][j] = vvx;
					Ufz_tmp[i][j] = vvz;
				}
			}
		}
	}
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Find_node_sem2fdm_cal(int order, int L, int Nxw, int Nzw, long zf, long node_num, double Dx, double Dz, double Z,
	long *node1_array, double **node_array, double **FDM_array_x, double **FDM_array_z, double **Fii_sem)
{
	int i, j, k, id1;
	int flgi1, flgi2, flgi3, flgi4;
	double sx, sz, fx1, fx2, fz1, fz2, flgd1;

	// Z = Z - 1500;
	// zf = zf - 5;

	if (fabs(FDM_array_z[0][zf + 1] - Z) >= 0.000001*Dz)
	{
		printf("Error! the trnas layer is not matching in Z direction! \n"), exit(1);
	}

	// Z = Z - 500;
	// zf = zf + 5;

	for(i = 0; i < node_num; i++)
	{
		id1 = node1_array[i] - 1;
		sx = node_array[id1][0];
		sz = node_array[id1][1];

		flgi1 = 0;

		if (fabs(Z - sz) < 0.0001*Dz)
		{
//			for (j = L; j < Nxw - L - 1; j++)
			for (j = L - 1; j < Nxw - L + 1; j++)
			{
				fx1 = FDM_array_x[j][0];
				fx2 = fx1 + Dx;

				if (sx >= fx1 && sx < fx2)
				{
					flgi1 = 1;
					break;
				}
			}
		}

		if (flgi1 == 1)
		{
			Fii_sem[i][0] = j;
			Fii_sem[i][1] = zf + 1;
		}
		else
		{
			Fii_sem[i][0] = -1.0;
			Fii_sem[i][1] = -1.0;
		}
	}
}

void Find_node_sem2fdm_pml(int order, int L, int Nx, int Nxw, long zf, long node_add, long pnode_num, double Dx, double Dz, double Z, 
	long *node2_array, long *node3_total, double **node_array, double **FDM_array_x, double **FDM_array_z, double **Fii_sem_pml)
{
	int i, j, k, id1, id2;
	int flgi1, flgi2, flgi3, flgi4;
	double sx, sz, fx1, fx2, fz1, fz2;

	if (fabs(FDM_array_z[0][zf + 1] - Z) >= 0.000001*Dz)
	{
		printf("Error! the trnas layer is not matching in Z direction! \n"), exit(1);
	}

	for(i = 0; i < pnode_num; i++)
	{
		id1 = node2_array[i] - 1;
		id2 = node3_total[id1];

		sx = node_array[id1][0];
		sz = node_array[id1][1];

		flgi1 = 0;

		if (fabs(Z - sz) < 0.0001*Dz && id2 <= 0)
		{
			for (j = 0; j < Nxw - 1; j++)
			{
				fx1 = FDM_array_x[j][0];
				fx2 = fx1 + Dx;

				if (sx >= fx1 && sx < fx2)
				{
					flgi1 = 1;
					break;
				}
			}
		}

		if (flgi1 == 1)
		{
			
			Fii_sem_pml[i][0] = j;
			Fii_sem_pml[i][1] = zf + 1;
		}
		else
		{
			Fii_sem_pml[i][0] = -1.0;
			Fii_sem_pml[i][1] = -1.0;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Trans_node_fdm2sem(int order, int Nxw, int Nzw, int L, long node_num, long *node1_array, float *Ux, float *Uz, float **Ufx, float **Ufz, double Dx,
	double Dz, double Z, double **FDM_array_x, double **FDM_array_z, double **node_array, double **Fii_sem)
{
	int i, j, k, id1;
	int fx, fz, flgi1, flgi2, flgi3;
	double coff, sx, sz, x1, y1;

	for (i = 0; i < node_num; i++)
	{
		id1 = node1_array[i] - 1;
		sx = node_array[id1][0];
		sz = node_array[id1][1];

		if (Fii_sem[i][0] >= 0)
		{
			fx = Fii_sem[i][0];
			fz = Fii_sem[i][1];

			if (fx < Nxw - 1)
			{
				coff = (sx - FDM_array_x[fx][fz]) / Dx;

				if (fabs(coff) <= 0.000001*Dx)
				{
					Ux[id1] = Ufx[fx][fz];
					Uz[id1] = Ufz[fx][fz];
				}
				else
				{
					Ux[id1] = coff*(Ufx[fx + 1][fz] - Ufx[fx][fz]) + Ufx[fx][fz];
					Uz[id1] = coff*(Ufz[fx + 1][fz] - Ufz[fx][fz]) + Ufz[fx][fz];
				}
			}
		}
	}
}


void Trans_node_fdm2sem_pml(int order, int Nx, int Nxw, int Nzw, int L, long node_add, long pnode_num, double Dx, double Dz, double Z, long *node2_array, 
	float *Upx1, float *Upx2, float *Upx3, float *Upz1, float *Upz2, float *Upz3, float *Ux, float *Uz, float **Ufx, float **Ufz, double **FDM_array_x, 
	double **FDM_array_z, double **node_array, double **Fii_sem_pml, float ****Ulx, float ****Ulz, float ****Urx, float ****Urz)
{
	int i, j, k, n, id1, id2, odr;
	int fx, fz, flgi1, flgi2, flgi3;
	double fplx11, fplx12, fplx13, fplx21, fplx22, fplx23, fplz11, fplz12, fplz13, fplz21, fplz22, fplz23;
	double coff, coff1, sx, sz, x1, y1, *arrfx, *arrfz, **arrfy;

	for(i = 0; i < pnode_num; i++)
	{
		id1 = node2_array[i] - 1;
		sx = node_array[id1][0];
		sz = node_array[id1][1];

		if (Fii_sem_pml[i][0] >= 0)
		{
			fx = Fii_sem_pml[i][0];
			fz = Fii_sem_pml[i][1];

			if (fx < L - 1)
			{
				id2 = fx;

				fplx11 = Ulx[1][0][id2][fz];  fplx21 = Ulx[1][0][id2 + 1][fz];
				fplx12 = Ulx[1][1][id2][fz];  fplx22 = Ulx[1][1][id2 + 1][fz];
				fplx13 = Ulx[1][2][id2][fz];  fplx23 = Ulx[1][2][id2 + 1][fz];

				fplz11 = Ulz[1][0][id2][fz];  fplz21 = Ulz[1][0][id2 + 1][fz];
				fplz12 = Ulz[1][1][id2][fz];  fplz22 = Ulz[1][1][id2 + 1][fz];
				fplz13 = Ulz[1][2][id2][fz];  fplz23 = Ulz[1][2][id2 + 1][fz];

				coff = (sx - FDM_array_x[fx][fz]) / Dx;

				if (fabs(coff) <= 0.000001*Dx)
				{
					Upx1[i] = fplx21;   Upz1[i] = fplz21;
					Upx3[i] = fplx22;   Upz3[i] = fplz22;
					Upx2[i] = fplx23;   Upz2[i] = fplz23;

					// Ux[id1] = Ufx[fx][fz];
					// Uz[id1] = Ufz[fx][fz];
					Ux[id1] = Upx1[i] + Upx2[i] + Upx3[i];
					Uz[id1] = Upz1[i] + Upz2[i] + Upz3[i];
				}
				else
				{
					Upx1[i] = coff*(fplx21 - fplx11) + fplx11;
					Upx3[i] = coff*(fplx22 - fplx12) + fplx12;
					Upx2[i] = coff*(fplx23 - fplx13) + fplx13;

					Upz1[i] = coff*(fplz21 - fplz11) + fplz11;
					Upz3[i] = coff*(fplz22 - fplz12) + fplz12;
					Upz2[i] = coff*(fplz23 - fplz13) + fplz13;

					// Ux[id1] = coff*(Ufx[fx + 1][fz] - Ufx[fx][fz]) + Ufx[fx][fz];
					// Uz[id1] = coff*(Ufz[fx + 1][fz] - Ufz[fx][fz]) + Ufz[fx][fz];
					Ux[id1] = Upx1[i] + Upx2[i] + Upx3[i];
					Uz[id1] = Upz1[i] + Upz2[i] + Upz3[i];
				}
			}
			else if (fx > Nx + L - 1)
			{
				// printf("%20.4e%20.4e%20.4e%20.4e\n", sx, sz, FDM_array_x[fx][fz], FDM_array_z[fx][fz]), sleep(1);
				if (sx <= FDM_array_x[Nx + L + 1][fz])
				{
					id2 = 0;

					Upx1[i] = Urx[1][0][id2][fz];   Upz1[i] = Urz[1][0][id2][fz];
					Upx3[i] = Urx[1][1][id2][fz];   Upz3[i] = Urz[1][1][id2][fz];
					Upx2[i] = Urx[1][2][id2][fz];   Upz2[i] = Urz[1][2][id2][fz];
				}
				else
				{
					id2 = fx - (Nx + L + 1);

					fplx11 = Urx[1][0][id2][fz];  fplx21 = Urx[1][0][id2 + 1][fz];
					fplx12 = Urx[1][1][id2][fz];  fplx22 = Urx[1][1][id2 + 1][fz];
					fplx13 = Urx[1][2][id2][fz];  fplx23 = Urx[1][2][id2 + 1][fz];

					fplz11 = Urz[1][0][id2][fz];  fplz21 = Urz[1][0][id2 + 1][fz];
					fplz12 = Urz[1][1][id2][fz];  fplz22 = Urz[1][1][id2 + 1][fz];
					fplz13 = Urz[1][2][id2][fz];  fplz23 = Urz[1][2][id2 + 1][fz];

					coff = (sx - FDM_array_x[fx][fz]) / Dx;

					if (fabs(coff) <= 0.000001*Dx)
					{
						Upx1[i] = fplx11;   Upz1[i] = fplz11;
						Upx3[i] = fplx12;   Upz3[i] = fplz12;
						Upx2[i] = fplx13;   Upz2[i] = fplz13;

						// Ux[id1] = Ufx[fx][fz];
						// Uz[id1] = Ufz[fx][fz];
						Ux[id1] = Upx1[i] + Upx2[i] + Upx3[i];
						Uz[id1] = Upz1[i] + Upz2[i] + Upz3[i];
					}
					else
					{
						Upx1[i] = coff*(fplx21 - fplx11) + fplx11;
						Upx3[i] = coff*(fplx22 - fplx12) + fplx12;
						Upx2[i] = coff*(fplx23 - fplx13) + fplx13;

						Upz1[i] = coff*(fplz21 - fplz11) + fplz11;
						Upz3[i] = coff*(fplz22 - fplz12) + fplz12;
						Upz2[i] = coff*(fplz23 - fplz13) + fplz13;

						// Ux[id1] = coff*(Ufx[fx + 1][fz] - Ufx[fx][fz]) + Ufx[fx][fz];
						// Uz[id1] = coff*(Ufz[fx + 1][fz] - Ufz[fx][fz]) + Ufz[fx][fz];
						Ux[id1] = Upx1[i] + Upx2[i] + Upx3[i];
						Uz[id1] = Upz1[i] + Upz2[i] + Upz3[i];
					}

				}
			}
		}
	}
}


///////////////////////////////////////////////// Calculate Ux123 and Uz123 in pml domain ////////////////////////////////////////////////
void Absorbing(int Nx,int Nz,int L,int order,double Dx,double Dz,double Dt,float **Ux2,float **Uz2,float **Vx,float **Vz,
			   float **lamda,float **mu,float **rho,double **coeff1,double **coeff2,float **BPxlx,float **BPxlz,float **BPzlx,
			   float **BPzlz,float **BPxrx,float **BPxrz,float **BPzrx,float **BPzrz,float **BPxuz,float **BPzuz,float **BPxdz,
			   float **BPzdz,float ****Ulx,float ****Ulz,float ****Urx,float ****Urz,float ****Uux,float ****Uuz,float ****Udx,
			   float ****Udz)
{
	double R0 = 1.0e-0;

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
	
	absorbingl(Nx,Nz,L,order,Dx,Dz,Dt,R0,Ux2,Uz2,Vx,Vz,lamda,mu,rho,coeff1,coeff2,BPxlx,BPxlz,BPzlx,BPzlz,Ulx,Ulz);

	absorbingr(Nx,Nz,L,order,Dx,Dz,Dt,R0,Ux2,Uz2,Vx,Vz,lamda,mu,rho,coeff1,coeff2,BPxrx,BPxrz,BPzrx,BPzrz,Urx,Urz);

	absorbingu(Nx,Nz,L,order,Dx,Dz,Dt,R0,Ux2,Uz2,Vx,Vz,lamda,mu,rho,coeff1,coeff2,BPxuz,BPzuz,Uux,Uuz);
}

void absorbingl(int Nx, int Nz, int L, int order, double Dx, double Dz, double Dt, double R0, float **Ux2, float **Uz2,
	float **Vx, float **Vz, float **lamda, float **mu, float **rho, double **coeff1, double **coeff2, float **BPxlx,
	float **BPxlz, float **BPzlx, float **BPzlz, float ****Ulx, float ****Ulz)
{
	int i, j, id1, id2, dxx, dzz, Nxw, Nzw, Nxw1, Nzw1, Nzw_L;
	double dix, diz, dix1, diz1, temp, t, Layx, Layz, Lx, Lz, Vp;

	Nxw = Nx + 2*L;  Nzw = Nz + L;
	Nxw1 = Nxw + 1;  Nzw1 = Nzw + 1;  Nzw_L = Nzw - L;

	id1 = Nzw1 / 2;

	Vp = sqrt( (lamda[L][id1] + 2.0*mu[L][id1]) / rho[L][id1] );

	Layx = 6000; Layz = 6000; Vp = 4000;

	for(i = 1; i < L; i++)
	{
		dxx = L - i;

		Lx = dxx*Dx;

		dix = 3.0*Vp*log10(1.0 / R0) / (2.0*Layx)*(Lx / Layx)*(Lx / Layx);

		dix1 = 6.0*Vp*log10(1.0 / R0) / (2.0*Layx)*(Lx / Layx) / Layx;

		// dix = 3.0*log10(1.0 / R0)* Vp/(2.0*L*Dx)*(dxx*Dx/(L*Dx))*(dxx*Dx/(L*Dx));

		// dix1=-6.0*log10(1.0/R0)*Vp/(2.0*L*Dx)*(dxx*Dx/(L*Dx))/(L*Dx);

		for(j = 1; j < Nzw; j++)
		{	
			if(j > Nzw_L)
			{
				dzz = j - Nzw_L;

				diz = 3.0*log10(1.0 / R0)*Vp / (2.0*L*Dz)*(dzz*Dz / (L*Dz))*(dzz*Dz / (L*Dz));

				diz1 = 6.0*log10(1.0 / R0)*Vp / (2.0*L*Dz)*(dzz*Dz / (L*Dz)) / (L*Dz);
			}
			else
			{
				diz = 0.0;
				diz1 = 0.0;
			}

			//======================================================================(1)
			t = diff(Nxw, Nzw, order, i, j, 1, 1, Dx, Dz, Vx, coeff1, coeff2);

			BPxlx[i][j] = (float)(((1.0 - 0.5*Dt*dix)*BPxlx[i][j] - dix1*Dt*t) / (1 + 0.5*Dt*dix));

			temp = Ulx[1][0][i][j];

			t = diff(Nxw, Nzw, order, i, j, 2, 1, Dx, Dz, Vx, coeff1, coeff2);

			Ulx[1][0][i][j] = (float)(((2.0 - Dt*Dt*dix*dix)*Ulx[1][0][i][j] + (Dt*dix - 1.0)*Ulx[0][0][i][j] +
				Dt*Dt*(lamda[i][j] + 2.0*mu[i][j]) / rho[i][j]*(t + BPxlx[i][j])) / (1.0 + Dt*dix));

			Ulx[0][0][i][j] = (float)temp;

			//======================================================================(2)
			t = diff(Nxw, Nzw, order, i, j, 1, 2, Dx, Dz, Vx, coeff1, coeff2);

			BPxlz[i][j] = (float)(((1.0 - 0.5*Dt*diz)*BPxlz[i][j] - diz1*Dt*t) / (1 + 0.5*Dt*diz));

			temp = Ulx[1][1][i][j];

			t = diff(Nxw, Nzw, order, i, j, 2, 2, Dx, Dz, Vx, coeff1, coeff2);

			Ulx[1][1][i][j] = (float)(((2.0 - Dt*Dt*diz*diz)*Ulx[1][1][i][j] + (Dt*diz - 1.0)*Ulx[0][1][i][j] +
				Dt*Dt*mu[i][j]/rho[i][j]*(t + BPxlz[i][j]))/(1.0 + Dt*diz));
			
			Ulx[0][1][i][j] = (float)temp;

			//=======================================================================(3)
			t = diff(Nxw, Nzw, order, i, j, 0, 0, Dx, Dz, Vz, coeff1, coeff2);
			
			temp = Ulx[1][2][i][j];
			
			Ulx[1][2][i][j] = (float)(((2.0 - Dt*Dt*dix*diz)*Ulx[1][2][i][j] + (0.5*Dt*(dix + diz) - 1.0)*Ulx[0][2][i][j] +
				Dt*Dt*(lamda[i][j] + mu[i][j]) / rho[i][j]*t) / (1.0 + 0.5*(dix + diz)*Dt));
			
			Ulx[0][2][i][j] = (float)(temp);

			Ux2[i][j] = Ulx[1][0][i][j] + Ulx[1][1][i][j] + Ulx[1][2][i][j];

			//========================================================================(1)
			t = diff(Nxw, Nzw, order, i, j, 1, 1, Dx, Dz, Vz, coeff1, coeff2);

			BPzlx[i][j] = (float)(((1.0 - 0.5*Dt*dix)*BPzlx[i][j] - dix1*Dt*t) / (1 + 0.5*Dt*dix));

			temp = Ulz[1][0][i][j];
			
			t = diff(Nxw, Nzw, order, i, j, 2, 1, Dx, Dz, Vz, coeff1, coeff2);

			Ulz[1][0][i][j] = (float)(((2.0 - Dt*Dt*dix*dix)*Ulz[1][0][i][j] + (Dt*dix - 1.0)*Ulz[0][0][i][j] +
				Dt*Dt*mu[i][j] / rho[i][j]*(t + BPzlx[i][j])) / (1.0 + Dt*dix));

			Ulz[0][0][i][j] = (float)temp;

			//=========================================================================(2)
			t = diff(Nxw, Nzw, order, i, j, 1, 2, Dx, Dz, Vz, coeff1, coeff2);

			BPzlz[i][j] = (float)(((1.0 - 0.5*Dt*diz)*BPzlz[i][j] - diz1*Dt*t) / (1.0 + 0.5*Dt*diz));

			temp = Ulz[1][1][i][j];

			t = diff(Nxw, Nzw, order, i, j, 2, 2, Dx, Dz, Vz, coeff1, coeff2);

			Ulz[1][1][i][j] = (float)(((2.0 - Dt*Dt*diz*diz)*Ulz[1][1][i][j] + (Dt*diz - 1.0)*Ulz[0][1][i][j] +
				Dt*Dt*(lamda[i][j] + 2.0*mu[i][j]) / rho[i][j]*(t + BPzlz[i][j])) / (1.0 + Dt*diz));
			
			Ulz[0][1][i][j] = (float)temp;

			//=========================================================================(3)
			t = diff(Nxw, Nzw, order, i, j, 0, 0, Dx, Dz, Vx, coeff1, coeff2);

			temp = Ulz[1][2][i][j];

			Ulz[1][2][i][j] = (float)(((2.0 - Dt*Dt*dix*diz)*Ulz[1][2][i][j] + (0.5*Dt*(diz + dix) - 1.0)*Ulz[0][2][i][j] +
				Dt*Dt*(lamda[i][j] + mu[i][j]) / rho[i][j]*t) / (1.0 + 0.5*(dix + diz)*Dt));

			Ulz[0][2][i][j] = (float)temp;

			Uz2[i][j] = Ulz[1][0][i][j] + Ulz[1][1][i][j] + Ulz[1][2][i][j];
		}
	}
}

void absorbingr(int Nx,int Nz,int L,int order,double Dx,double Dz,double Dt,double R0,float **Ux2,float **Uz2,float **Vx,
		float **Vz,float **lamda,float **mu,float **rho,double **coeff1,double **coeff2,float **BPxrx,float **BPxrz,
		float **BPzrx,float **BPzrz,float ****Urx,float ****Urz)
{

	int i, j, id, t1, t2, t3, t4;
	int dxx, dzz, Nxw, Nzw, Nxw1, Nzw1, NxL1, Nzw_L;
	long k, z;
	float tmp1,tmp2, tmp3, tmp4, tmpf1, tmpf2, tmpf3, tmpf4;
	double dix,diz,dix1,diz1,temp,t,Vp;
	double Layx, Layz, Lx, Lz;

	Nxw = Nx+2*L;
	Nzw = Nz+L;
	Nxw1 = Nxw+1;
	Nzw1 = Nzw+1;
	NxL1 = Nx+L+1;
	Nzw_L = Nzw-L;

	// Vp = sqrt((lamda[Nxw][Nzw1/2]+2.0*mu[Nxw][Nzw1/2])/rho[Nxw][Nzw1/2]);
	Layx = 6000; Layz = 6000; Vp = 4000;

	for(i = NxL1; i < Nxw;i++)
	{
		dxx = i-Nx-L;

		Lx = dxx*Dx;

		dix = 3.0*Vp*log10(1.0 / R0) / (2.0*Layx)*(Lx / Layx)*(Lx / Layx);

		dix1 = 6.0*Vp*log10(1.0 / R0) / (2.0*Layx)*(Lx / Layx) / Layx;

		id = dxx - 1;

		for(j = 1; j < Nzw; j++)
		{
			if(j > Nzw_L)
			{
				dzz = j - Nzw_L;

				diz = 3.0*log10(1.0/R0)*Vp/(2.0*L*Dz)*(dzz*Dz/(L*Dz))*(dzz*Dz/(L*Dz));

				diz1 = 6.0*log10(1.0/R0)*Vp/(2.0*L*Dz)*(dzz*Dz/(L*Dz))/(L*Dz);
			}
			else
			{
				diz = 0.0;
				diz1 = 0.0;
			}

			//======================================================================(0)
			t1 = order;  t2 = order;

			if(t1>i)
			{
				t1 = i;
			}
			else if(t1>Nxw-i)
			{
				t1 = Nxw - i;
			}

	
			if(t2>j)
			{
				t2 = j;
			}
			else if(t2>Nzw-j)
			{
				t2 = Nzw-j;
			}

			tmpf1 = 0.0;  tmpf2 = 0.0;
			tmpf3 = 0.0;  tmpf4 = 0.0;

			//======================================================================(1)
			t = diff(Nxw,Nzw,order,i,j,1,1,Dx,Dz,Vx,coeff1,coeff2);

			BPxrx[id][j] = (float)(((1.0-0.5*Dt*dix)*BPxrx[id][j]-dix1*Dt*t)/(1+0.5*Dt*dix));////P中间变量

			temp = Urx[1][0][id][j];

			t = diff(Nxw,Nzw,order,i,j,2,1,Dx,Dz,Vx,coeff1,coeff2);

			Urx[1][0][id][j] = (float)( ( (2.0 - Dt*Dt*dix*dix)*Urx[1][0][id][j]  +(Dt*dix - 1.0) *Urx[0][0][id][j]
				+ Dt*Dt*(lamda[i][j] + 2.0*mu[i][j]) / rho[i][j]*(t + BPxrx[id][j])) / (1.0 + Dt*dix));

			Urx[0][0][id][j] = (float)temp;

			//======================================================================(2)
			t = diff(Nxw,Nzw,order,i,j,1,2,Dx,Dz,Vx,coeff1,coeff2);
			
			BPxrz[id][j] = (float)(((1.0-0.5*Dt*diz)*BPxrz[id][j]-diz1*Dt*t)/(1+0.5*Dt*diz));

			temp = Urx[1][1][id][j];

			t = diff(Nxw,Nzw,order,i,j,2,2,Dx,Dz,Vx,coeff1,coeff2);

			Urx[1][1][id][j] = (float)(((2.0-Dt*Dt*diz*diz)*Urx[1][1][id][j]+(Dt*diz-1.0)*Urx[0][1][id][j]
				+Dt*Dt*(mu[i][j]/rho[i][j])*(t+BPxrz[id][j]))/(1.0+Dt*diz));

			Urx[0][1][id][j] = (float)temp;

			//=======================================================================(3)
			for(k = 1;k <= t1;k++)
			{
				tmp1 = 0.0;  tmp2 = 0.0;

				for(z = 1;z <= t2; z++)
				{
					tmp1+=coeff1[t2-1][z-1]*(lamda[i+k][j+z]+lamda[i+k][j-z])*(Vz[i+k][j+z]-Vz[i+k][j-z])/Dz;
				}

				for(z=1;z<=t2;z++)
				{
					tmp2+=coeff1[t2-1][z-1]*(lamda[i-k][j+z]+lamda[i-k][j-z])*(Vz[i-k][j+z]-Vz[i-k][j-z])/Dz;
				}

				tmpf3+=coeff1[t1-1][k-1]*(tmp1-tmp2)/(2.0*Dx);
			}

			for(k=1;k<=t2;k++)
			{
				tmp1=0.0;
				tmp2=0.0;

				for(z=1;z<=t1;z++)
				{
					tmp1+=coeff1[t1-1][z-1]*(mu[i+z][j+k]+mu[i-z][j+k])*(Vz[i+z][j+k]-Vz[i-z][j+k])/Dx;
				}

				for(z=1;z<=t1;z++)
				{
					tmp2+=coeff1[t1-1][z-1]*(mu[i+z][j-k]+mu[i-z][j-k])*(Vz[i+z][j-k]-Vz[i-z][j-k])/Dx;				
				}

				tmpf4+=coeff1[t2-1][k-1]*(tmp1-tmp2)/(2.0*Dz);
			}

			temp = Urx[1][2][id][j];

			Urx[1][2][id][j]=(float)(((2.0-Dt*Dt*dix*diz)*Urx[1][2][id][j]+(0.5*Dt*(dix+diz)-1.0)*Urx[0][2][id][j]
				+Dt*Dt*(tmpf3+tmpf4)/rho[i][j])/(1.0+0.5*(dix+diz)*Dt));

			Urx[0][2][id][j] = (float)temp;
			
			Ux2[i][j] = Urx[1][0][id][j] + Urx[1][1][id][j] + Urx[1][2][id][j];

			//========================================================================(0)
			tmpf1 = 0.0;  tmpf2 = 0.0;
			tmpf3 = 0.0;  tmpf4 = 0.0;

			//========================================================================(1)
			t = diff(Nxw,Nzw,order,i,j,1,1,Dx,Dz,Vz,coeff1,coeff2);

			BPzrx[id][j] = (float)(((1.0-0.5*Dt*dix)*BPzrx[id][j]-dix1*Dt*t)/(1+0.5*Dt*dix));

			temp = Urz[1][0][id][j];

			t = diff(Nxw,Nzw,order,i,j,2,1,Dx,Dz,Vz,coeff1,coeff2);

			Urz[1][0][id][j] = (float)(((2.0-Dt*Dt*dix*dix)*Urz[1][0][id][j]+(Dt*dix-1.0)*Urz[0][0][id][j]+
				Dt*Dt*(mu[i][j]/rho[i][j])*(t+BPzrx[id][j]))/(1.0+Dt*dix));

			Urz[0][0][id][j] = (float)temp;

			//=========================================================================(2)
			t = diff(Nxw,Nzw,order,i,j,1,2,Dx,Dz,Vz,coeff1,coeff2);

			BPzrz[id][j] = (float)(((1.0-0.5*Dt*diz)*BPzrz[id][j]-diz1*Dt*t)/(1.0+0.5*Dt*diz));

			temp = Urz[1][1][id][j];

			t = diff(Nxw,Nzw,order,i,j,2,2,Dx,Dz,Vz,coeff1,coeff2);

			Urz[1][1][id][j] = (float)(((2.0-Dt*Dt*diz*diz)*Urz[1][1][id][j]+(Dt*diz-1.0)*Urz[0][1][id][j]+
				Dt*Dt*(lamda[i][j]+2.0*mu[i][j])/rho[i][j]*(t+BPzrz[id][j]))/(1.0+Dt*diz));

			Urz[0][1][id][j] = (float)temp;
			//=========================================================================(3)
			for(k=1;k<=t2;k++)
			{
				tmp1=0.0;
				tmp2=0.0;

				for(z=1;z<=t1;z++)
				{
					tmp1+=coeff1[t1-1][z-1]*(Vx[i+z][j+k]-Vx[i-z][j+k])/Dx;
				}
				for(z=1;z<=t1;z++)
				{
					tmp2+=coeff1[t1-1][z-1]*(Vx[i+z][j-k]-Vx[i-z][j-k])/Dx;
				}

				tmpf3+=coeff1[t2-1][k-1]*((lamda[i][j])*tmp1-(lamda[i][j])*tmp2)/(Dz);
			}

			for(k=1;k<=t1;k++)
			{
				tmp1=0.0; tmp2=0.0;

				for(z=1;z<=t2;z++)
				{
					tmp1+=coeff1[t2-1][z-1]*(mu[i+k][j+z]+mu[i+k][j-z])*(Vx[i+k][j+z]-Vx[i+k][j-z])/Dz;
				}

				for(z = 1; z <= t2; z++)
				{
					tmp2 += coeff1[t2-1][z-1]*(mu[i-k][j+z]+mu[i-k][j-z])*(Vx[i-k][j+z]-Vx[i-k][j-z])/Dz;					
				}

				tmpf4 += coeff1[t1-1][k-1]*(tmp1-tmp2)/(2.0*Dx);
			}

			temp = Urz[1][2][id][j];

			Urz[1][2][id][j]=(float)(((2.0-Dt*Dt*dix*diz)*Urz[1][2][id][j]+(0.5*Dt*(diz+dix)-1.0)*Urz[0][2][id][j]
				+Dt*Dt*(tmpf3+tmpf4)/rho[i][j])/(1.0+0.5*(dix+diz)*Dt));

			Urz[0][2][id][j]=(float)temp;

			Uz2[i][j] = Urz[1][0][id][j] + Urz[1][1][id][j] + Urz[1][2][id][j];
		}
	}   
}

void absorbingu(int Nx,int Nz,int L,int order,double Dx,double Dz,double Dt,double R0,float **Ux2,float **Uz2,float **Vx,
		float **Vz,float **lamda,float **mu,float **rho,double **coeff1,double **coeff2,float **BPxuz,float **BPzuz,
		float ****Uux,float ****Uuz)
{
	double diz,diz1,temp,t,Vp;
	int i,j,id1,id2,dzz,Nxw,NxL,NzL1,Nzw;
	
	Nxw=Nx+2*L;
	Nzw=Nz+L;
	NxL=Nx+L;
	NzL1=Nz+1;

	Vp=sqrt((lamda[Nxw/2][NzL1]+2.0*mu[Nxw/2][NzL1])/rho[Nxw/2][NzL1]);
//	Vp = 4000;


	for(i=L;i<=NxL;i++)
	{
		id1=i-L;
		for(j=NzL1;j<Nzw;j++)
		{
			dzz=j-Nz;
			diz=3.0*log10(1.0/R0)*Vp/(2.0*L*Dz)*(dzz*Dz/(L*Dz))*(dzz*Dz/(L*Dz));
			diz1=6.0*log10(1.0/R0)*Vp/(2.0*L*Dz)*(dzz*Dz/(L*Dz))/(L*Dz);
			id2=dzz-1;

			//===========================================================================(1)
			t=diff(Nxw,Nzw,order,i,j,2,1,Dx,Dz,Vx,coeff1,coeff2);
			temp=Uux[1][0][id1][id2];
			Uux[1][0][id1][id2]=(float)(2.0*Uux[1][0][id1][id2]-Uux[0][0][id1][id2]+Dt*Dt*(lamda[i][j]+
				2.0*mu[i][j])/rho[i][j]*t);
			Uux[0][0][id1][id2]=(float)temp;
			//============================================================================(2)
			t=diff(Nxw,Nzw,order,i,j,1,2,Dx,Dz,Vx,coeff1,coeff2);
			BPxuz[id1][id2]=(float)(((1.0-0.5*Dt*diz)*BPxuz[id1][id2]-diz1*Dt*t)/(1.0+0.5*Dt*diz));
			temp=Uux[1][1][id1][id2];
			t=diff(Nxw,Nzw,order,i,j,2,2,Dx,Dz,Vx,coeff1,coeff2);
			Uux[1][1][id1][id2]=(float)(((2.0-Dt*Dt*diz*diz)*Uux[1][1][id1][id2]+
				(Dt*diz-1.0)*Uux[0][1][id1][id2]+Dt*Dt*(mu[i][j]/rho[i][j])*(t+BPxuz[id1][id2]))/(1.0+Dt*diz));
			Uux[0][1][id1][id2]=(float)temp;
			//============================================================================(3)
			t=diff(Nxw,Nzw,order,i,j,0,0,Dx,Dz,Vz,coeff1,coeff2);
			temp=Uux[1][2][id1][id2];
			Uux[1][2][id1][id2]=(float)((2.0*Uux[1][2][id1][id2]+(0.5*Dt*diz-1.0)*Uux[0][2][id1][id2]
				+Dt*Dt*(lamda[i][j]+mu[i][j])/rho[i][j]*t)/(1.0+0.5*diz*Dt));
			Uux[0][2][id1][id2]=(float)temp;

			Ux2[i][j]=Uux[1][0][id1][id2]+Uux[1][1][id1][id2]+Uux[1][2][id1][id2];
			
			//=============================================================================(1)
			t=diff(Nxw,Nzw,order,i,j,2,1,Dx,Dz,Vz,coeff1,coeff2);
			temp=Uuz[1][0][id1][id2];
			Uuz[1][0][id1][id2]=(float)(2.0*Uuz[1][0][id1][id2]-Uuz[0][0][id1][id2]+(mu[i][j]/rho[i][j])*Dt*Dt*t);
			Uuz[0][0][id1][id2]=(float)temp;
			//=============================================================================(2)
			t=diff(Nxw,Nzw,order,i,j,1,2,Dx,Dz,Vz,coeff1,coeff2);
			BPzuz[id1][id2]=(float)(((1.0-0.5*Dt*diz)*BPzuz[id1][id2]-diz1*Dt*t)/(1.0+0.5*Dt*diz));
			temp=Uuz[1][1][id1][id2];
			t=diff(Nxw,Nzw,order,i,j,2,2,Dx,Dz,Vz,coeff1,coeff2);
			Uuz[1][1][id1][id2]=(float)(((2.0-Dt*Dt*diz*diz)*Uuz[1][1][id1][id2]+(Dt*diz-1.0)*
				Uuz[0][1][id1][id2]+Dt*Dt*(lamda[i][j]+2.0*mu[i][j])/rho[i][j]*(t+BPzuz[id1][id2]))/(1.0+Dt*diz));
			Uuz[0][1][id1][id2]=(float)temp;
			//=============================================================================(3)
			t=diff(Nxw,Nzw,order,i,j,0,0,Dx,Dz,Vx,coeff1,coeff2);
			temp=Uuz[1][2][id1][id2];
			Uuz[1][2][id1][id2]=(float)((2.0*Uuz[1][2][id1][id2]+(0.5*Dt*diz-1.0)*Uuz[0][2][id1][id2]
				+Dt*Dt*(lamda[i][j]+mu[i][j])/rho[i][j]*t)/(1.0+0.5*diz*Dt));
			Uuz[0][2][id1][id2]=(float)temp;
			Uz2[i][j]=Uuz[1][0][id1][id2]+Uuz[1][1][id1][id2]+Uuz[1][2][id1][id2];
		}
	}
}
