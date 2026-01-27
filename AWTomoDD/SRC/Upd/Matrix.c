#include "SEM.h"
///////////////////////////////////composite mass matrix////////////////////////////
///////////////////////////////////合成质量矩阵/////////////////////////////
void Mass_matrix(short num,long elem,float *rho,double **node_array,long **Elem_array,float *Mass,float *Mass_tmp)
{
	long i,j,j1,j2,index,*id,N2,n1,n2,n3,n4;
	double add=0.0,**B,*W,*P,*dx,*dz,coeff;

	N2=num*num;
	B=Arraydouble2(num,num);
	id=Arraylong1(N2);
	W=Arraydouble1(num);
	P=Arraydouble1(num);
	dx=Arraydouble1(4);
	dz=Arraydouble1(4);

	Weight(num-1,W);
	Position(num-1,P);

	n1=0;								//初始的四列
	n2=num*(num-1);
	n3=num*num-1;
	n4=num-1;

	for(i=0;i<elem;i++)
	{
		for(j=0;j<N2;j++)
		{
			id[j]=Elem_array[i][j]-1;
		}

		dx[0]=node_array[id[n1]][0];
		dz[0]=node_array[id[n1]][1];
		
		dx[1]=node_array[id[n2]][0];
		dz[1]=node_array[id[n2]][1];
		
		dx[2]=node_array[id[n3]][0];
		dz[2]=node_array[id[n3]][1];
		
		dx[3]=node_array[id[n4]][0];
		dz[3]=node_array[id[n4]][1];

//		printf("%e %e\n",dx[0],dz[0]);
//		printf("%e %e\n",dx[1],dz[1]);
//		printf("%e %e\n",dx[2],dz[2]);
//		printf("%e %e\n",dx[3],dz[3]);

		local_mass(num,dx,dz,W,P,B);				//W、P、B是算雅可比矩阵的一些系数W积分权P积分位置B坐标变化因子
		coeff=rho[i];
		
		for(j1=0;j1<num;j1++)
		{
			for(j2=0;j2<num;j2++)
			{
				index=j1*num+j2;			//0--35
				index=id[index];
				Mass[index]+=(float)(B[j1][j2]*coeff);	//局部计算之后整合到一起
				Mass_tmp[index]+=(float)(B[j1][j2]);
		//if(node_array[index][0]>=0&&node_array[index][0]<=20&&node_array[index][1]>=0&&node_array[index][1]<=20)
		//printf("%f,%f,%f,%f\n",node_array[index][0],node_array[index][1],B[j1][j2],Mass_tmp[index]);
			}
		}
		
	}
	/*for(i=0;i<elem;i++)
	{
		for(j1=0;j1<num;j1++)
		{
			for(j2=0;j2<num;j2++)
			{
				index=j1*num+j2;
				index=Elem_array[i][index]-1;
				printf("Mass %f,%f,%f\n",node_array[index][0],node_array[index][1],Mass_tmp[index]);
			}
		}
		getchar();		
	}*/
	Freedouble2(num,B);
	free(id);
	free(W);
	free(P);
	free(dx);
	free(dz);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////local mass matrix//////////////////////////////////////
void local_mass(short int num,double *dx,double *dz,double *W,double *P,double **B)
{
	double Ja,*A11,*A12,*A21,*A22,a11,a12,a21,a22;
	long i1,i2;
	
	A11=&a11;
	A12=&a12;
	A21=&a21;
	A22=&a22;

	for(i1=0;i1<num;i1++)
	{
		for(i2=0;i2<num;i2++)
		{
			Jacobi_for(dx,dz,P[i1],P[i2],A11,A12,A21,A22);	//雅可比矩阵
			Ja=a11*a22-a12*a21;
			B[i1][i2]=Ja*W[i1]*W[i2];
		}
	}
}
//////////////////////////////////////////Stiff matrix/////////////////////////////////
void Stiff_matrix(short num,long elem_num,long **Elem_array,long *Elem1_count,double **node_array,float **KKK,float **Kxx,float **Kxz,float **Kzx,float **Kzz)
{
	short *Cin;
	long i,it,i1,i2;
	double *W,*Ax,*Az,*P;

	P=Arraydouble1(num);
	W=Arraydouble1(num);
	Cin=Arrayshort1(4);
	Position(num-1,P);
	Weight(num-1,W);
	Ax=Arraydouble1(4);
	Az=Arraydouble1(4);

	Cin[0]=0;
	Cin[1]=num*(num-1);
	Cin[2]=num*num-1;
	Cin[3]=num-1;
	
	for(i=0;i<elem_num;i++)
	{
		it=Elem1_count[i]-1;
		for(i1=0;i1<4;i1++)
		{
			i2=Elem_array[it][Cin[i1]]-1;
			Ax[i1]=node_array[i2][0];
			Az[i1]=node_array[i2][1];
		}
		local_stiff(num,i,W,P,Ax,Az,KKK,Kxx,Kxz,Kzx,Kzz);
	} 
	
	free(P); free(Ax);
	free(W); free(Az);
	free(Cin);
}  
////////////////////////////////////local stiff matrix////////////////////////
///////////////////////////////////////////////////////////////////////////////
void local_stiff(short num,long ei,double *W,double *P,double *Ax,double *Az,float **KKK,float **Kxx,float **Kxz,float **Kzx,float **Kzz)
{
	double px,pz,a11,a12,a21,a22;
	double Nx1,Nx2,Nx3,Nx4,Nz1,Nz2,Nz3,Nz4;
	double b,add,Ja;
	short i1,i2;
	
	for(i1=0;i1<num;i1++)
	{
		px=P[i1];
		for(i2=0;i2<num;i2++)
		{
			pz=P[i2];
			Nx1=-0.25*(1.0-pz);
			Nx2=0.25*(1.0-pz);
			Nx3=0.25*(1.0+pz);
			Nx4=-0.25*(1.0+pz);

			Nz1=-0.25*(1.0-px);
			Nz2=-0.25*(1.0+px);
			Nz3=0.25*(1.0+px);
			Nz4=0.25*(1.0-px);

			a11=Nx1*Ax[0]+Nx2*Ax[1]+Nx3*Ax[2]+Nx4*Ax[3];
			a12=Nx1*Az[0]+Nx2*Az[1]+Nx3*Az[2]+Nx4*Az[3];
			a21=Nz1*Ax[0]+Nz2*Ax[1]+Nz3*Ax[2]+Nz4*Ax[3];
			a22=Nz1*Az[0]+Nz2*Az[1]+Nz3*Az[2]+Nz4*Az[3];

			b=a11*a22-a12*a21;
			Ja=1.0/b;

			KKK[ei][i1*num+i2]=b;
			Kxx[ei][i1*num+i2]=a22*Ja;
			Kzx[ei][i1*num+i2]=-a12*Ja;
			Kxz[ei][i1*num+i2]=-a21*Ja;
			Kzz[ei][i1*num+i2]=a11*Ja;
		}
	}
}
////////////////////////////Stiff_assemble////////////////////////////////////////////////////////
void Stiff_assemble(int nn,long n1,long n2,double add,long *Jtmp,double *Stmp,float *Stiff,long *Istiff,long *Jstiff)
{
	int flag,flag1,flag2;
	long i,k,n11;

	flag=0;
	flag1=0;
	flag2=0;
	n11=n1+1;

	for(i=0;i<nn;i++)
	{
		Jtmp[i]=-1;
		Stmp[i]=0.0;
	}

	for(i=Istiff[n1];i<Istiff[n11];i++)
		if(Jstiff[i]>-1)
		{
			flag=1;
			break;
		}

	if(flag==1)
	{
		for(i=Istiff[n1];i<Istiff[n11];i++)
		{
			if(Jstiff[i]==n2)
			{
				Stiff[i]+=(float)add;
				flag1=1;
				break;
			}
		}	
		if(flag1==0)
		{
			flag2=0;
			for(i=Istiff[n1];i<Istiff[n11];i++)
			{
				if(n2<Jstiff[i])
				{
					flag2=1;
					for(k=i;k<Istiff[n11]-1;k++)
					{
						Stmp[k-i]=Stiff[k];
						Jtmp[k-i]=Jstiff[k];
					}
					Stiff[i]=(float)add;
					Jstiff[i]=n2;
					for(k=i;k<Istiff[n11]-1;k++)
					{
						Stiff[k+1]=(float)Stmp[k-i];
						Jstiff[k+1]=Jtmp[k-i];
					}
					break;
				}
			}		
			if(flag2==0)
			{
				for(i=Istiff[n1];i<Istiff[n11];i++)
					if(Jstiff[i]==-1)
					{	
						Jstiff[i]=n2;
						Stiff[i]+=(float)add;
						break;
					}
			}
		}
	}
	else
	{
		Stiff[Istiff[n1]]=(float)add;
		Jstiff[Istiff[n1]]=n2;
	}
	
}

///////////////////////////////////////local matrix load///////////////////////////////
void local_load(short int num,double **node_array,long *id,double *W,double *Nx,double *Nz,double *Shape,double **di1,double **di2,
				float ****NT1,float ****NT2,double **Bxx,double **Bzz)
{
	double a11,a12,a21,a22,b,b11,b12,b21,b22,add;
	long i1,i2,j,j1,j2,N2;

	N2=num*num;

	for(i1=0;i1<num;i1++)
	{
		for(i2=0;i2<num;i2++)
		{

			a11=0.0;
			a12=0.0;
			a21=0.0;
			a22=0.0;

			for(j1=0;j1<num;j1++)
				for(j2=0;j2<num;j2++)
				{	
					j=j1*num+j2;
					a11+=NT1[i1][i2][j1][j2]*node_array[id[j]][0];
					a12+=NT1[i1][i2][j1][j2]*node_array[id[j]][1];
					a21+=NT2[i1][i2][j1][j2]*node_array[id[j]][0];
					a22+=NT2[i1][i2][j1][j2]*node_array[id[j]][1];
				}


			b=a11*a22-a12*a21;

			b11=a22/b;
			b12=-a12/b;
			b21=-a21/b;
			b22=a11/b;
			
			for(j1=0;j1<num;j1++)
				for(j2=0;j2<num;j2++)
				{
					j=j1*num+j2;						
					Nx[j]=NT1[i1][i2][j1][j2]*b11+NT2[i1][i2][j1][j2]*b12;
					Nz[j]=NT1[i1][i2][j1][j2]*b21+NT2[i1][i2][j1][j2]*b22;			
				}
		
			for(j=0;j<N2;j++)
				Shape[j]=0.0;

			j=i1*num+i2;
			Shape[j]=1.0;

			for(j1=0;j1<N2;j1++)
				for(j2=0;j2<N2;j2++)
				{
					di1[j1][j2]=Shape[j1]*Nx[j2];
					di2[j1][j2]=Shape[j1]*Nz[j2];
				}
					
			add=b*W[i2]*W[i1];
					
			for(j1=0;j1<N2;j1++)
				for(j2=0;j2<N2;j2++)
				{
					 Bxx[j1][j2]+=di1[j1][j2]*add;
					 Bzz[j1][j2]+=di2[j1][j2]*add;
				}
		}
	}

}
/////////////////////////////////////////Jacobi matrix/////////////////////////////////
void Jacobi_for(double *dx,double *dz,double x,double z,double *A11,double *A12,double *A21,double *A22)
{
	double a11,a12,a21,a22,Nx1,Nx2,Nx3,Nx4,Nz1,Nz2,Nz3,Nz4;

	Nx1=-1.0/4.0*(1.0-z);
	Nx2=1.0/4.0*(1.0-z);
	Nx3=1.0/4.0*(1.0+z);
	Nx4=-1.0/4.0*(1.0+z);

	Nz1=-1.0/4.0*(1.0-x);
	Nz2=-1.0/4.0*(1.0+x);
	Nz3=1.0/4.0*(1.0+x);
	Nz4=1.0/4.0*(1.0-x);

	a11=Nx1*dx[0]+Nx2*dx[1]+Nx3*dx[2]+Nx4*dx[3];
	a12=Nx1*dz[0]+Nx2*dz[1]+Nx3*dz[2]+Nx4*dz[3];
	a21=Nz1*dx[0]+Nz2*dx[1]+Nz3*dx[2]+Nz4*dx[3];
	a22=Nz1*dz[0]+Nz2*dz[1]+Nz3*dz[2]+Nz4*dz[3];

	*A11=a11;
	*A12=a12;
	*A21=a21;
	*A22=a22;
}

/////////////////////////////////////////Load Vector//////////////////////////////////////
void Loadf1(short num,long elem_num,double Sx,double Sz,double **node_array,long *node1_total,long **Elem_array,long *Elem1_count,float *Loadx,float *Loadz)
{
	long i,i1,j,j1,j2,*id,N2,n1,n2,n3,n4,in1,in2;
	double add=0.0,**B,*W,*P,*dx,*dz,*Cf,xx,zz,coeff,DD;

	N2=num*num;
	B=Arraydouble2(num,num);
	id=Arraylong1(N2);
	W=Arraydouble1(num);
	P=Arraydouble1(num);
	dx=Arraydouble1(4);
	dz=Arraydouble1(4);
	Cf=Arraydouble1(N2);

	Weight(num-1,W);
	Position(num-1,P);

	n1=0;
	n2=num*(num-1);
	n3=num*num-1;
	n4=num-1;

	DD=200.0;

	for(i=0;i<elem_num;i++)
	{
		i1=Elem1_count[i]-1;
		for(j=0;j<N2;j++)
		{
			id[j]=Elem_array[i1][j]-1;
			xx=node_array[id[j]][0]-Sx;
			zz=node_array[id[j]][1]-Sz;
			Cf[j]=exp(-(xx*xx+zz*zz)/DD);
		}

		dx[0]=node_array[id[n1]][0];
		dz[0]=node_array[id[n1]][1];
		
		dx[1]=node_array[id[n2]][0];
		dz[1]=node_array[id[n2]][1];
		
		dx[2]=node_array[id[n3]][0];
		dz[2]=node_array[id[n3]][1];
		
		dx[3]=node_array[id[n4]][0];
		dz[3]=node_array[id[n4]][1];

		local_mass(num,dx,dz,W,P,B);
		for(j1=0;j1<num;j1++)
		{
			for(j2=0;j2<num;j2++)
			{
				in1=j1*num+j2;
				in2=node1_total[id[in1]]-1;
				Loadz[in2]+=(float)(B[j1][j2]*Cf[in1]);
				Loadx[in2]+=0.0;
			}
		}
	}
	Freedouble2(num,B);
	free(id);
	free(W);
	free(P);
	free(dx);
	free(dz);
	free(Cf);
}


void Loadf2(short num,long elem_num,double Sx,double Sz,double **node_array,long *node1_total,long **Elem_array,long *Elem1_count,float *Loadx,float *Loadz)
{
	long i,i1,j,j1,j2,*id,N2,n1,n2,n3,n4,in1,in2;
	double add=0.0,**B,*W,*P,*dx,*dz,*Cf1,*Cf2,xx,zz,coeff,DD;

	N2=num*num;
	B=Arraydouble2(num,num);
	id=Arraylong1(N2);
	W=Arraydouble1(num);
	P=Arraydouble1(num);
	dx=Arraydouble1(4);
	dz=Arraydouble1(4);
	Cf1=Arraydouble1(N2);
	Cf2=Arraydouble1(N2);

	Weight(num-1,W);
	Position(num-1,P);

	n1=0;
	n2=num*(num-1);
	n3=num*num-1;
	n4=num-1;
	DD=400.0;

	for(i=0;i<elem_num;i++)
	{
		i1=Elem1_count[i]-1;
		for(j=0;j<N2;j++)
		{
			id[j]=Elem_array[i1][j]-1;
			xx=node_array[id[j]][0]-Sx;
			zz=node_array[id[j]][1]-Sz;
			Cf1[j]=-2.0*xx/DD*exp(-(xx*xx+zz*zz)/DD);
			Cf2[j]=-2.0*zz/DD*exp(-(xx*xx+zz*zz)/DD);
		}

		dx[0]=node_array[id[n1]][0];
		dz[0]=node_array[id[n1]][1];
		
		dx[1]=node_array[id[n2]][0];
		dz[1]=node_array[id[n2]][1];
		
		dx[2]=node_array[id[n3]][0];
		dz[2]=node_array[id[n3]][1];
		
		dx[3]=node_array[id[n4]][0];
		dz[3]=node_array[id[n4]][1];

		local_mass(num,dx,dz,W,P,B);

		for(j1=0;j1<num;j1++)
		{
			for(j2=0;j2<num;j2++)
			{
				in1=j1*num+j2;
				in2=node1_total[id[in1]]-1;
				Loadx[in2]+=(float)(B[j1][j2]*Cf1[in1]);
				Loadz[in2]+=(float)(B[j1][j2]*Cf2[in1]);
			}
		}
	}
	Freedouble2(num,B);
	free(id);
	free(W);
	free(P);
	free(dx);
	free(dz);
	free(Cf1);
	free(Cf2);
}

void Loadf_fdm(int Nxw, int Nzw, double Sx, double Sz, double Dx, double Dz, float **Sx_fdm, float **Sz_fdm, double **FDM_array_x, double **FDM_array_z)
{
	int i, j, num1, num2, si, sj, id1, id2;
	double mon, rr_fdm, rr, rr_sx, rr_sz, tmp_sx, tmp_sz, ss[100], par_ss[100];
	double dd, dx, xx;

	rr_fdm = 2 * Dx;
	mon = 1.0e+11;

	for(i = 0; i < Nxw; i++)
	{
		for(j = 0; j < Nzw; j++)
		{
			rr_sx = FDM_array_x[i][j] - Sx;
			rr_sz = FDM_array_z[i][j] - Sz;

			rr = sqrt(rr_sx*rr_sx + rr_sz*rr_sz);

			if (rr <= rr_fdm)
			{
				tmp_sx = FDM_array_x[i][j];
				tmp_sz = FDM_array_z[i][j];

				si = i;
				sj = j;

				rr_fdm = rr;
			}
		}
	}

	num1 = 2;
	num2 = num1*2 + 1;

	dx = Dx;

	for (i = 0; i < num2; i++)
	{
		id1 = i - num1;

		xx = id1*dx;

		dd = (num1 / 2.0)*dx;

		ss[i] = pseudo_delt(xx, dd);

		par_ss[i] = pseudo_delt_diff(xx, dd);
	}


	id1 = 0;
	for (i = si - num1; i <= si + num1; i++)
	{
		id2 = 0;
		for (j = sj - num1; j <= sj + num1; j++)
		{
			Sz_fdm[i][j] = mon*ss[id1]*par_ss[id2];
			id2++;
		}
		id1++;
	}
}

void Loadf_sem1(short num, long elem_num, double len, double Sx, double Sz, double **node_array, long *node1_total, long **Elem_array, long *Elem1_count,
	float *Loadx, float *Loadz, double **Tensor)
{
	long i, j, i1, j1, j2, n1, n2, n3, n4, nn, N2, id1, id2, id3, id4, id[100];
	double xx, zz, fx, fz, coeff, DD, rr, rr_sem, rr_sx, rr_sz, tmp_sx, tmp_sz, ss[100], par_ss[100], cf1[200], cf2[200];
	double *W, *P, *dx, *dz, **B;

	N2 = num * num;

	W = Arraydouble1(num);
	P = Arraydouble1(num);
	dx = Arraydouble1(4);
	dz = Arraydouble1(4);
	B = Arraydouble2(num, num);

	Weight(num - 1, W);
	Position(num - 1, P);

	rr_sem = 2 * len;
	
	n1 = 0;
	n2 = num*(num - 1);
	n3 = num*num - 1;
	n4 = num - 1;

	/*
	for(i = 0; i < elem_num; i++)
	{
		nn = Elem1_count[i] - 1;
		id1 = Elem_array[nn][n1] - 1;
		id2 = Elem_array[nn][n2] - 1;
		id3 = Elem_array[nn][n3] - 1;
		id4 = Elem_array[nn][n4] - 1;

		dx[0] = node_array[id1][0];    dz[0] = node_array[id1][1];
		dx[1] = node_array[id2][0];    dz[1] = node_array[id2][1];
		dx[2] = node_array[id3][0];    dz[2] = node_array[id3][1];
		dx[3] = node_array[id4][0];    dz[3] = node_array[id4][1];

		xx = (dx[0] + dx[1] + dx[2] + dx[3]) / 4.0;
		zz = (dz[0] + dz[1] + dz[2] + dz[3]) / 4.0;

		rr_sx = xx - Sx;
		rr_sz = zz - Sz;

		rr = sqrt(rr_sx*rr_sx + rr_sz*rr_sz);

		if (rr <= rr_sem)
		{
			tmp_sx = xx;
			tmp_sz = zz;

			rr_sem = rr;
		}
	}
	*/
	// Sx = tmp_sx; Sz = tmp_sz;

	// printf("%10.4f%10.4f\n", Sx, Sz); sleep(1000);

	len = len / 1.0;

	for(i = 0; i < elem_num; i++)
	{
		i1 = Elem1_count[i] - 1;

		for(j = 0; j < N2; j++)
		{
			id[j] = Elem_array[i1][j] - 1;

			xx = node_array[id[j]][0] - Sx;
			zz = node_array[id[j]][1] - Sz;

			fx = pseudo_delt_diff(xx, len);
			fz = pseudo_delt(zz, len);

			cf1[j] = fx*fz;

			fx = pseudo_delt(xx, len);
			fz = pseudo_delt_diff(zz, len);

			cf2[j] = fx*fz;
		}

		dx[0] = node_array[id[n1]][0];
		dz[0] = node_array[id[n1]][1];

		dx[1] = node_array[id[n2]][0];
		dz[1] = node_array[id[n2]][1];

		dx[2] = node_array[id[n3]][0];
		dz[2] = node_array[id[n3]][1];

		dx[3] = node_array[id[n4]][0];
		dz[3] = node_array[id[n4]][1];

		local_mass(num, dx, dz, W, P, B);

		for(j1 = 0; j1 < num; j1++)
		{
			for(j2 = 0; j2 < num; j2++)
			{
				id1 = j1*num + j2;
				id2 = node1_total[id[id1]] - 1;

				Loadx[id2] += -1.0*B[j1][j2]*(cf1[id1]*Tensor[0][0] + cf2[id1]*Tensor[0][1]);
				Loadz[id2] += -1.0*B[j1][j2]*(cf1[id1]*Tensor[1][0] + cf2[id1]*Tensor[1][1]);
			}
		}
	}

	free(W);   free(P);
	free(dx);  free(dz);

	Freedouble2(num, B);
}

void Loadf_sem_syn(short num, long elem_num, double len, double Sx, double Sz, double **node_array, long *node1_total, long **Elem_array, long *Elem1_count,
	float *Loadx, float *Loadz, double **Tensor)
{
	long i, j, i1, j1, j2, n1, n2, n3, n4, nn, N2, id1, id2, id3, id4, id[100];
	double xx, zz, fx, fz, coeff, DD, rr, rr_sem, rr_sx, rr_sz, tmp_sx, tmp_sz, ss[100], par_ss[100], cf1[200], cf2[200];
	double *W, *P, *dx, *dz, **B;

	N2 = num * num;

	W = Arraydouble1(num);
	P = Arraydouble1(num);
	dx = Arraydouble1(4);
	dz = Arraydouble1(4);
	B = Arraydouble2(num, num);

	Weight(num - 1, W);
	Position(num - 1, P);

	rr_sem = 2 * len;
	
	n1 = 0;
	n2 = num*(num - 1);
	n3 = num*num - 1;
	n4 = num - 1;

	len = len / (num+1);

	for(i = 0; i < elem_num; i++)
	{
		i1 = Elem1_count[i] - 1;

		for(j = 0; j < N2; j++)
		{
			id[j] = Elem_array[i1][j] - 1;

			xx = node_array[id[j]][0] - Sx;
			zz = node_array[id[j]][1] - Sz;

			// fx = pseudo_delt_filter(xx, len);
			// fz = pseudo_delt_filter(zz, len);
			fx = pseudo_delt_h4(xx, len);
			fz = pseudo_delt_h4(zz, len);

			cf1[j] = fx*fz;

//			fx = pseudo_delt_h4(xx, len);
//			fz = pseudo_delt_h4(zz, len);

			cf2[j] = fx*fz;
		}

		dx[0] = node_array[id[n1]][0];
		dz[0] = node_array[id[n1]][1];

		dx[1] = node_array[id[n2]][0];
		dz[1] = node_array[id[n2]][1];

		dx[2] = node_array[id[n3]][0];
		dz[2] = node_array[id[n3]][1];

		dx[3] = node_array[id[n4]][0];
		dz[3] = node_array[id[n4]][1];

		local_mass(num, dx, dz, W, P, B);

		for(j1 = 0; j1 < num; j1++)
		{
			for(j2 = 0; j2 < num; j2++)
			{
				id1 = j1*num + j2;
				id2 = node1_total[id[id1]] - 1;

				Loadx[id2] += -1.0*B[j1][j2]*(cf1[id1]*Tensor[0][0] + cf2[id1]*Tensor[0][1]);
				Loadz[id2] += -1.0*B[j1][j2]*(cf1[id1]*Tensor[1][0] + cf2[id1]*Tensor[1][1]);
			}
		}
	}

	free(W);   free(P);
	free(dx);  free(dz);

	Freedouble2(num, B);
}

void Loadf_sem_adj(short num, long elem_num, double len, double Sx, double Sz, double **node_array, long *node1_total, long **Elem_array, long *Elem1_count,
	float *Loadx, float *Loadz, double **Tensor)
{
	long i, j, i1, j1, j2, n1, n2, n3, n4, nn, N2, id1, id2, id3, id4, id[100];
	double xx, zz, fx, fz, coeff, DD, rr, rr_sem, ss[100], par_ss[100], cf1[200], cf2[200], mm[2][2];
	double *W, *P, *dx, *dz, **B;


	N2 = num * num;

	W = Arraydouble1(num);
	P = Arraydouble1(num);
	dx = Arraydouble1(4);
	dz = Arraydouble1(4);
	B = Arraydouble2(num, num);

	Weight(num - 1, W);
	Position(num - 1, P);


	for (i = 0; i < 2; i++)
	{
		for (j = 0; j < 2; j++)
		{
			mm[i][j] = Tensor[i][j];
			if (mm[i][j] != 0.0)
			{
				mm[i][j] = -1.0;
			}
		}
	}

	rr_sem = 2*len;
	
	n1 = 0;
	n2 = num*(num - 1);
	n3 = num*num - 1;
	n4 = num - 1;

	len = len / (num+1);

	for(i = 0; i < elem_num; i++)
	{
		i1 = Elem1_count[i] - 1;

		for(j = 0; j < N2; j++)
		{
			id[j] = Elem_array[i1][j] - 1;

			xx = node_array[id[j]][0] - Sx;
			zz = node_array[id[j]][1] - Sz;

			// fx = pseudo_delt_filter(xx, len);
			// fz = pseudo_delt_filter(zz, len);
			fx = pseudo_delt_h4(xx, len);
			fz = pseudo_delt_h4(zz, len);

			cf1[j] = fx*fz;

//			fx = pseudo_delt_h4(xx, len);
//			fz = pseudo_delt_h4(zz, len);

			cf2[j] = fx*fz;
		}

		dx[0] = node_array[id[n1]][0];
		dz[0] = node_array[id[n1]][1];

		dx[1] = node_array[id[n2]][0];
		dz[1] = node_array[id[n2]][1];

		dx[2] = node_array[id[n3]][0];
		dz[2] = node_array[id[n3]][1];

		dx[3] = node_array[id[n4]][0];
		dz[3] = node_array[id[n4]][1];

		local_mass(num, dx, dz, W, P, B);

		for(j1 = 0; j1 < num; j1++)
		{
			for(j2 = 0; j2 < num; j2++)
			{
				id1 = j1*num + j2;
				id2 = node1_total[id[id1]] - 1;

				Loadx[id2] += -1.0*B[j1][j2]*(cf1[id1]*mm[0][0] + cf2[id1]*mm[0][1]);
				Loadz[id2] += -1.0*B[j1][j2]*(cf1[id1]*mm[1][0] + cf2[id1]*mm[1][1]);
			}
		}
	}

	free(W);   free(P);
	free(dx);  free(dz);

	Freedouble2(num, B);
}
//////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////Find position////////////////////////////////////////////
void Find_local_position(short num,short nn,long elem_num,double dlen,double Sx,double Sz,long *Elem1_count,
	long **Elem_array,double **node_array,long *ii,double *lx,double *lz)
{
	short N3,nn1,flag,*Cin;
	int i,it,i0,i1,i2,i3,id,size1;
	double xx,yy,zz,ax,ay,az,a1,a2,a3,a4,a5,a6,a7,a8,tmp0,tmp,*P1,Dis;
	double x1,x2,x3,x4,x5,x6,x7,x8,y1,y2,y3,y4,y5,y6,y7,y8,z1,z2,z3,z4,z5,z6,z7,z8,px,py,pz;

	Dis=10.0*dlen;
	nn1=20*nn-1;						//199将每个元素分成小元素用于计算

	N3=num*num*num;
	P1=Arraydouble1(nn1);
	Cin=Arrayshort1(4);

	Cin[0]=0;
	Cin[1]=num*(num-1);
	Cin[2]=num*num-1;
	Cin[3]=num-1;

	for(i=0;i<nn1;i++)
	{
		P1[i]=-1.0+2.0/(nn1-1)*i;
	}

	ax=-2.0;
	ay=-2.0;
	az=-2.0;
	i0=-1;
	tmp0=1.0e10;

	for(i=0;i<elem_num;i++)
	{
		it=Elem1_count[i]-1;
		i1=Elem_array[it][Cin[0]]-1;
		x1=node_array[i1][0];
		z1=node_array[i1][1];

		i1=Elem_array[it][Cin[1]]-1;
		x2=node_array[i1][0];
		z2=node_array[i1][1];

		i1=Elem_array[it][Cin[2]]-1;
		x3=node_array[i1][0];
		z3=node_array[i1][1];

		i1=Elem_array[it][Cin[3]]-1;
		x4=node_array[i1][0];
		z4=node_array[i1][1];
		
		xx=0.25*(x1+x2+x3+x4)-Sx;
		zz=0.25*(z1+z2+z3+z4)-Sz;
		tmp=sqrt(xx*xx+zz*zz);

		if(tmp<Dis)
		{
			for(i1=1;i1<nn1-1;i1++)
			{
				px=P1[i1];
				for(i2=1;i2<nn1-1;i2++)
				{
					pz=P1[i2];
					a1=(1.0-px)*(1.0-pz)*0.25;
					a2=(1.0+px)*(1.0-pz)*0.25;
					a3=(1.0+px)*(1.0+pz)*0.25;
					a4=(1.0-px)*(1.0+pz)*0.25;
										
					xx=a1*x1+a2*x2+a3*x3+a4*x4-Sx;
					zz=a1*z1+a2*z2+a3*z3+a4*z4-Sz;
					tmp=sqrt(xx*xx+zz*zz);
					
					if(tmp<tmp0)
					{
						*ii=i;
						*lx=px;
						*lz=pz;
						tmp0=tmp;
					}
				}
			}
		}
	}
	free(P1);
	free(Cin);
}


void Loadf(short num,long node_num,long Si,double Sx,double Sz,long *Elem1_count,long *node1_total,long **Elem_array,double **NT,double **node_array,double **M,
	float *Loadx,float *Loadz)
{	
	short *Cin;
	long i,j,id,it,i1,i2;
	double *Vx,*Vz,*Dx,*Dz,*P,*W,*Ax,*Az,px,pz;
	double Nx1,Nx2,Nx3,Nx4,Nz1,Nz2,Nz3,Nz4,b11,b12,b21,b22;
	double b,Ja,a11,a12,a21,a22,mxx,mxz,mzx,mzz,Gxx,Gxz,Gzx,Gzz;

	for(i=0;i<node_num;i++)
	{
		Loadx[i]=0.0;
		Loadz[i]=0.0;
	}

	Vx=Arraydouble1(num);
	Vz=Arraydouble1(num);
	Dx=Arraydouble1(num);
	Dz=Arraydouble1(num);
	W=Arraydouble1(num);
	P=Arraydouble1(num);
	
	Weight(num-1,W);
	Position(num-1,P);

	Ax=Arraydouble1(4);
	Az=Arraydouble1(4);
	Cin=Arrayshort1(4);

	Cin[0]=0;
	Cin[1]=num*(num-1);
	Cin[2]=num*num-1;
	Cin[3]=num-1;
	it=Elem1_count[Si]-1;

	for(i=0;i<4;i++)
	{
		i1=Elem_array[it][Cin[i]]-1;
		Ax[i]=node_array[i1][0];
		Az[i]=node_array[i1][1];
//		printf("1 %f,%f\n",Ax[i],Az[i]);
	}

//	printf("1 %f,%f\n",Sx,Sz);

	for(i=0;i<num;i++)
	{
		Vx[i]=Lagrange_val(num,i,P,Sx);
		Vz[i]=Lagrange_val(num,i,P,Sz);
		Dx[i]=Lagrange_diff(num,i,P,Sx);
		Dz[i]=Lagrange_diff(num,i,P,Sz);
//printf("2 %f,%f,%f,%f\n",Vx[i],Vz[i],Dx[i],Dz[i]);
	}
//////getchar();
	px=Sx;
	pz=Sz;

    	Nx1=-0.25*(1.0-pz);
	Nx2=0.25*(1.0-pz);
	Nx3=0.25*(1.0+pz);
	Nx4=-0.25*(1.0+pz);

	Nz1=-0.25*(1.0-px);
	Nz2=-0.25*(1.0+px);
	Nz3=0.25*(1.0+px);
	Nz4=0.25*(1.0-px);
//printf("3 %f,%f,%f,%f\n",Nx1,Nx2,Nx3,Nx4);
//printf("3 %f,%f,%f,%f\n",Nz1,Nz2,Nz3,Nz4);

	b11=Nx1*Ax[0]+Nx2*Ax[1]+Nx3*Ax[2]+Nx4*Ax[3];
	b12=Nx1*Az[0]+Nx2*Az[1]+Nx3*Az[2]+Nx4*Az[3];
	b21=Nz1*Ax[0]+Nz2*Ax[1]+Nz3*Ax[2]+Nz4*Ax[3];
	b22=Nz1*Az[0]+Nz2*Az[1]+Nz3*Az[2]+Nz4*Az[3];
//printf("4 %f,%f,%f,%f\n",b11,b12,b21,b22);

	b=b11*b22-b12*b21;
	Ja=1.0/b;
	
	a11=b22*Ja;//&1/x1
	a12=-b12*Ja;//&2/x1
	a21=-b21*Ja;//&1/x2
	a22=b11*Ja;//&2/x2
//printf("5 %f,%f,%f,%f,%f\n",a11,a12,a21,a22,Ja);

	
	Gxx=M[0][0]*a11+M[0][1]*a21;//G11
	Gxz=M[0][0]*a12+M[0][1]*a22;//G12
	Gzx=M[1][0]*a11+M[1][1]*a21;//G21
	Gzz=M[1][0]*a12+M[1][1]*a22;//G22
//printf("6 %e,%e,%e,%e\n\n",Gxx,Gxz,Gzx,Gzz);


	for(i1=0;i1<num;i1++)
	{
		for(i2=0;i2<num;i2++)
		{
			id=Elem_array[it][i1*num+i2]-1;
//printf("%d,%d,%f,%f ",i1,i2,node_array[id][0],-1.0*node_array[id][1]);
			id=node1_total[id]-1;
			Loadx[id]=(Gxx*Dx[i1]*Vz[i2]+Gxz*Vx[i1]*Dz[i2]);
			Loadz[id]=(Gzx*Dx[i1]*Vz[i2]+Gzz*Vx[i1]*Dz[i2]);
//printf("%e,%e\n",Loadx[id],Loadz[id]);
		}
	}
//getchar();
	free(Vx);  free(Vz); free(P);
	free(Dx);  free(Dz); free(Cin);
	free(Ax);  free(Az);
}
void Loadf_s(short num,long node_num,long Si,double Sx,double Sz,long *Elem1_count,long *node1_total,long **Elem_array,double **NT,double **node_array,double **M,
	float **Loadx_s,float **Loadz_s,long Ni)
{	
	short *Cin;
	long i,j,id,it,i1,i2;
	double *Vx,*Vz,*Dx,*Dz,*P,*W,*Ax,*Az,px,pz;
	double Nx1,Nx2,Nx3,Nx4,Nz1,Nz2,Nz3,Nz4,b11,b12,b21,b22;
	double b,Ja,a11,a12,a21,a22,mxx,mxz,mzx,mzz,Gxx,Gxz,Gzx,Gzz;

	for(i=0;i<node_num;i++)
	{
		Loadx_s[Ni][i]=0.0;
		Loadz_s[Ni][i]=0.0;
	}

	Vx=Arraydouble1(num);
	Vz=Arraydouble1(num);
	Dx=Arraydouble1(num);
	Dz=Arraydouble1(num);
	W=Arraydouble1(num);
	P=Arraydouble1(num);
	
	Weight(num-1,W);
	Position(num-1,P);


	Ax=Arraydouble1(4);
	Az=Arraydouble1(4);
	Cin=Arrayshort1(4);
//printf("%f,%f\n",Sx,Sz);
	Cin[0]=0;
	Cin[1]=num*(num-1);
	Cin[2]=num*num-1;
	Cin[3]=num-1;
	it=Elem1_count[Si]-1;

	for(i=0;i<4;i++)
	{
		i1=Elem_array[it][Cin[i]]-1;
		Ax[i]=node_array[i1][0];
		Az[i]=node_array[i1][1];
//printf("1 %f,%f\n",Ax[i],Az[i]);
	}
//printf("1 %f,%f\n",Sx,Sz);
	for(i=0;i<num;i++)
	{
		Vx[i]=Lagrange_val(num,i,P,Sx);
		Vz[i]=Lagrange_val(num,i,P,Sz);
		Dx[i]=Lagrange_diff(num,i,P,Sx);
		Dz[i]=Lagrange_diff(num,i,P,Sz);
//printf("2 %f,%f,%f,%f\n",Vx[i],Vz[i],Dx[i],Dz[i]);
	}
//////getchar();
	px=Sx;
	pz=Sz;

    Nx1=-0.25*(1.0-pz);
	Nx2=0.25*(1.0-pz);
	Nx3=0.25*(1.0+pz);
	Nx4=-0.25*(1.0+pz);

	Nz1=-0.25*(1.0-px);
	Nz2=-0.25*(1.0+px);
	Nz3=0.25*(1.0+px);
	Nz4=0.25*(1.0-px);
//printf("3 %f,%f,%f,%f\n",Nx1,Nx2,Nx3,Nx4);
//printf("3 %f,%f,%f,%f\n",Nz1,Nz2,Nz3,Nz4);

	b11=Nx1*Ax[0]+Nx2*Ax[1]+Nx3*Ax[2]+Nx4*Ax[3];
	b12=Nx1*Az[0]+Nx2*Az[1]+Nx3*Az[2]+Nx4*Az[3];
	b21=Nz1*Ax[0]+Nz2*Ax[1]+Nz3*Ax[2]+Nz4*Ax[3];
	b22=Nz1*Az[0]+Nz2*Az[1]+Nz3*Az[2]+Nz4*Az[3];
//printf("4 %f,%f,%f,%f\n",b11,b12,b21,b22);

	b=b11*b22-b12*b21;
	Ja=1.0/b;
	
	a11=b22*Ja;//&1/x1
	a12=-b12*Ja;//&2/x1
	a21=-b21*Ja;//&1/x2
	a22=b11*Ja;//&2/x2
//printf("5 %f,%f,%f,%f,%f\n",a11,a12,a21,a22,Ja);

	
	Gxx=M[0][0]*a11+M[0][1]*a21;//G11
	Gxz=M[0][0]*a12+M[0][1]*a22;//G12
	Gzx=M[1][0]*a11+M[1][1]*a21;//G21
	Gzz=M[1][0]*a12+M[1][1]*a22;//G22
//printf("6 %e,%e,%e,%e\n\n",Gxx,Gxz,Gzx,Gzz);


	for(i1=0;i1<num;i1++)
	{
		for(i2=0;i2<num;i2++)
		{
			id=Elem_array[it][i1*num+i2]-1;
//printf("%d,%d,%f,%f ",i1,i2,node_array[id][0],-1.0*node_array[id][1]);
			id=node1_total[id]-1;
			Loadx_s[Ni][id]=(Gxx*Dx[i1]*Vz[i2]+Gxz*Vx[i1]*Dz[i2]);
			Loadz_s[Ni][id]=(Gzx*Dx[i1]*Vz[i2]+Gzz*Vx[i1]*Dz[i2]);
//printf("%e,%e\n",Loadx[id],Loadz[id]);
		}
	}
//getchar();
	free(Vx);  free(Vz); free(P);
	free(Dx);  free(Dz); free(Cin);
	free(Ax);  free(Az);
}

double pseudo_delt(double x, double dx)
{
	int i, j, k;
	double ff;

	if (x <= -2*dx)
	{
		ff = 0.0;
	}
	else if(x <= -dx)
	{
		ff = 0.5*(x + 2*dx)*(x + 2*dx)/(dx*dx);
	}
	else if(x <= dx)
	{
		ff = -0.5*(x + 2*dx)*(x + 2*dx)/(dx*dx) + 2*(x + 2*dx)/dx - 1.0;
	}
	else if (x <= 2*dx)
	{
		ff = 0.5*(x + 2*dx)*(x + 2*dx)/(dx*dx) - 4*(x + 2*dx)/dx + 8.0;
	}
	else
	{
		ff = 0.0;
	}

	ff = ff*(1 / (2*dx));

	return ff;
}

double pseudo_delt_diff(double x, double dx)
{
	int i, j, k;
	double ff;

	if (x <= -2*dx)
	{
		ff = 0.0;
	}
	else if(x <= -dx)
	{
		ff = (x + 2*dx)/(dx*dx);
	}
	else if(x <= dx)
	{
		ff = -(x + 2*dx)/(dx*dx) + 2/dx;
	}
	else if (x <= 2*dx)
	{
		ff = (x + 2*dx)/(dx*dx) - 4/dx;
	}
	else
	{
		ff = 0.0;
	}

	ff = ff*(1 / (2*dx));

	return ff;
}

double pseudo_delt_h4(double x, double dx)
{
	double x1, x2, xx, ff;

	// x1 = 2.0*dx; x2 = 4.0*dx;

	xx = fabs(x / dx);

	if (xx < 2.0)
	{
		ff = (1.0/32.0)*(16.0 - 4.0*xx - 4.0*xx*xx + xx*xx*xx);
	}
	else if (xx >= 2.0 && xx < 4.0)
	{
		ff = (1.0/96.0)*(48.0 - 44.0*xx + 12.0*xx*xx - xx*xx*xx);
	}
	else
	{
		ff = 0.0;
	}

	ff = ff / dx;

	return ff;
}

double pseudo_delt_filter(double x, double dx)
{
	int i, count1, flgi;
	
	char flnm1[100];
	
	double dd, rr, xx, ff, pi, **delta;

	sprintf(flnm1, "../Input/delta.txt");

	count1 = count_file(flnm1);

	delta = Arraydouble2(count1, 2);

	rdfile(flnm1, count1, 2, delta);
	
	xx = x / dx;

	dd = 10.0 / count1 + 1.0e-10;
	
	flgi = 0;
	for(i = 0; i < count1; i++)
	{
		rr = fabs(xx - delta[i][0]);

		if(rr < dd)
		{
			ff = delta[i][1];
			flgi = 1;
			break;
		}
	}
	
	if(flgi == 0)
	{
		ff = 0.0;
	}
	
	ff = ff / dx;
	
	Freedouble2(count1, delta);
	return ff;
}

void Get_Ricker(int Ni, int Nt0, double dt0, double *Ri)
{
	int i, j, Nt, flgi1;
	char file1[100];
	double t, rr, **ri;

	sprintf(file1, "../Input/st_function");

	Nt = 50000;

	ri = Arraydouble2(Nt, 2);

	rdfile(file1, Nt, 2, ri);

	for (i = 0; i < Nt0; i++)
	{
		Ri[i] = 0.0;
	}

	for (i = 0; i < Nt0; i++)
	{
		t = i*dt0;
		flgi1 = 0;

		for (j = 0; j < Nt; j++)
		{
			rr = fabs(ri[j][0] - t);

			if (rr <= dt0)
			{
				flgi1 = 1;
				break;
			}
		}

		if (flgi1 == 1)
		{
			Ri[i] = ri[j][1];
		}
		else
		{
			Ri[i] = 0.0;
		}
	}
}
