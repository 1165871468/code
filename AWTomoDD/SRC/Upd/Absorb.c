#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include "Absorb.h"
#define min1 1.0e-3
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////node information of PML domain/////////////////////////////
void Boundary_node_inf1(short int num,long node_add,long elem_num,long elem,double X,double Z,long **Elem_array,double **node_array,
						long *node_total_pml,long *node_total_pmlc,long *node_numpml_point,long *node_numpmlc_point)
{
	double dx,dz;
	long i,j,id,node_numpml;
	short int N2;
	
	N2=num*num;
	node_numpml=0;

	for(i=elem_num;i<elem;i++)
		for(j=0;j<N2;j++)
		{
			id=Elem_array[i][j]-1;
			node_total_pml[id]=id+1;
		}
	
	for(i=elem_num;i<elem;i++)
		for(j=0;j<N2;j++)
		{
			id=Elem_array[i][j]-1;
			dx=node_array[id][0];
			dz=node_array[id][1];
			if(dx<-min1||dx>X+min1||dz>Z+min1)
				node_total_pmlc[id]=id+1;
		}		

	for(i=0;i<node_add;i++)
	{
		if(node_total_pml[i]>0)
			node_numpml++;
	}
	*node_numpml_point=node_numpml;
	node_numpml=0;
	for(i=0;i<node_add;i++)
	{
		if(node_total_pmlc[i]>0)
			node_numpml++;
	}
	*node_numpmlc_point=node_numpml;
}
//////////////////////////////////////////node PML to Calculate///////////////////////////////////////
void Boundary_node_inf2(long node,long node_numpml,long node_numpmlc,long *node_total_pml,
						long *node_total_pmlc,long *node_pml,long *node_pmlc)
{
	long i,i1,count;
	count=0;
	for(i=0;i<node;i++)
		if(node_total_pml[i]>0)
		{
			node_pml[count]=node_total_pml[i];
			count++;
		}

	count=0;
	for(i=0;i<node;i++)
		if(node_total_pmlc[i]>0)
		{
			node_pmlc[count]=node_total_pmlc[i];
			count++;
		}

	for(i=0;i<node;i++)
	{
		node_total_pml[i]=0;
		node_total_pmlc[i]=0;
	}

	for(i=0;i<node_numpml;i++)
	{
		i1=node_pml[i]-1;
		node_total_pml[i1]=i+1;
	}

	for(i=0;i<node_numpmlc;i++)
	{
		i1=node_pmlc[i]-1;
		node_total_pmlc[i1]=i+1;
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////nodes damp///////////////////////////////////////////////
void node_damp(long node_numpmlc,double Vp,double Lay,double X,double Z,long *node_pmlc,double **node_array,
			   float *dampxz,float *dampxz1)
{
	long i,i1,id1,id2;
	double Lx,Lz,R0,x1,z1;
	R0=1.0e-6;
	for(i=0;i<node_numpmlc;i++)
	{
		
		i1=node_pmlc[i]-1;

		x1=node_array[i1][0];
		z1=node_array[i1][1];

		id1=2*i;
		id2=id1+1;
			
		if(x1<0&&z1<Z)
		{
			Lx=-x1;
			dampxz[id1]=(float)(3.0*log10(1.0/R0)*Vp/(2.0*Lay)*(Lx/Lay)*(Lx/Lay));
			dampxz[id2]=0.0;
			
			dampxz1[id1]=(float)(-6.0*log10(1.0/R0)*Vp/(2.0*Lay)*(Lx/(Lay))/(Lay));
			dampxz1[id2]=0.0;
		}
		else if(x1<0&&z1>=Z)
		{
			Lx=-x1;
			Lz=z1-Z;
			dampxz[id1]=(float)(3.0*log10(1.0/R0)*Vp/(2.0*Lay)*(Lx/Lay)*(Lx/Lay));
			dampxz[id2]=(float)(3.0*log10(1.0/R0)*Vp/(2.0*Lay)*(Lz/Lay)*(Lz/Lay));

			dampxz1[id1]=(float)(-6.0*log10(1.0/R0)*Vp/(2.0*Lay)*(Lx/(Lay))/(Lay));
			dampxz1[id2]=(float)(6.0*log10(1.0/R0)*Vp/(2.0*Lay)*(Lz/(Lay))/(Lay));
		}
		else if(x1>=0&&x1<=X&&z1>Z)
		{
			Lz=z1-Z;
			dampxz[id1]=0.0;
			dampxz[id2]=(float)(3.0*log10(1.0/R0)*Vp/(2.0*Lay)*(Lz/Lay)*(Lz/Lay));

			dampxz1[id1]=0.0;
			dampxz1[id2]=(float)(6.0*log10(1.0/R0)*Vp/(2.0*Lay)*(Lz/(Lay))/(Lay));
		}
		else if(x1>X&&z1<Z)
		{
			Lx=x1-X;
			dampxz[id1]=(float)(3.0*log10(1.0/R0)*Vp/(2.0*Lay)*(Lx/Lay)*(Lx/Lay));
			dampxz[id2]=0.0;

			dampxz1[id1]=(float)(6.0*log10(1.0/R0)*Vp/(2.0*Lay)*(Lx/(Lay))/(Lay));
			dampxz1[id2]=0.0;
		}
		else if(x1>X&&z1>=Z)
		{
			Lx=x1-X;
			Lz=z1-Z;
			dampxz[id1]=(float)(3.0*log10(1.0/R0)*Vp/(2.0*Lay)*(Lx/Lay)*(Lx/Lay));
			dampxz[id2]=(float)(3.0*log10(1.0/R0)*Vp/(2.0*Lay)*(Lz/Lay)*(Lz/Lay));

			dampxz1[id1]=(float)(6.0*log10(1.0/R0)*Vp/(2.0*Lay)*(Lx/(Lay))/(Lay));
			dampxz1[id2]=(float)(6.0*log10(1.0/R0)*Vp/(2.0*Lay)*(Lz/(Lay))/(Lay));
		}
	} 
}
/////////////////////////////////////////PML initial///////////////////////////////////////////////
void PML_initial(int Node2,short int max,long node_numpml,long Elem_num,long Elem,long **Elem_array,long *node_total,
				 long *Istiff,long *Jstiff)
{
	long i,j,j1,j2,*Jtmp,*nd,*nd1;

	nd=Arraylong1(Node2);
	nd1=Arraylong1(Node2);
	
	for(i=0;i<Istiff[node_numpml];i++)
		Jstiff[i]=-1;

	Jtmp=Arraylong1(max);

	for(i=Elem_num;i<Elem;i++)
	{
		for(j=0;j<Node2;j++)
		{
			nd[j]=Elem_array[i][j]-1;
			nd1[j]=node_total[nd[j]]-1;
		}

		for(j1=0;j1<Node2;j1++)
			for(j2=0;j2<Node2;j2++)
				Ini_assemble(max,nd1[j1],nd1[j2],Jtmp,Istiff,Jstiff);
	}
	free(Jtmp);
	free(nd);
	free(nd1);
}
//////////////////////////////////////////////////////////////////////////////////
void PML_count(long node_num,long *Istiff_ini,long *Jstiff_ini,long *Istiff)
{
	long i,j,i1,count;

	count=0;
	for(i=0;i<node_num;i++)
	{
		i1=i+1;
		for(j=Istiff_ini[i];j<Istiff_ini[i1];j++)
		{
			if(Jstiff_ini[j]>=0)
				count++;
		}

		Istiff[i1]=count;
	}
}
////////////////////////////////////////Kx1 matrix//////////////////////////////////////////////////
void Kx1_matrix(short int num,short int nn,long node_num,long elem_num,long elem,float *lamda,float *mu,long *node3_total,
				long **Elem_array,double **node_array,float ****NT1,float ****NT2,long *IA,long *JA,float *Kx1)
{
	long i,i1,j1,j2,*Jtmp,*id,*id1,NP,N2;
	double *Stmp,*W,*N1,cof,add;
	float **B;

	for(i=0;i<IA[node_num];i++)
		JA[i]=-1;

	N2=num*num;
	id=Arraylong1(N2);
	id1=Arraylong1(N2);
	N1=Arraydouble1(N2);
	B=Arrayfloat2(N2,N2);
	
	Jtmp=Arraylong1(nn);
	Stmp=Arraydouble1(nn);

	W=Arraydouble1(num);

	NP=num-1;
	Weight(NP,W);


	for(i=elem_num;i<elem;i++)
	{
		for(i1=0;i1<N2;i1++)
		{
			id[i1]=Elem_array[i][i1]-1;
			id1[i1]=node3_total[id[i1]]-1;
		}

		local_Kx1(num,node_array,id,W,N1,NT1,NT2,B);

		cof=lamda[i]+2.0*mu[i];
		for(j1=0;j1<N2;j1++)
			for(j2=0;j2<N2;j2++)
			{
				add=cof*B[j1][j2];
				Stiff_assemble(nn,id1[j1],id1[j2],add,Jtmp,Stmp,Kx1,IA,JA);

			}
	}      
	free(Jtmp);
	free(Stmp);
	free(W);
	free(id);
	free(id1);
	free(N1);
	Freefloat2(N2,B);

} 
//////////////////////////////////////////local kx1//////////////////////////////////
void local_Kx1(short int num,double **node_array,long *id,double *W,double *N1,float ****NT1,float ****NT2,float **B)
{
	double a11,a12,a21,a22,b,b11,b12,add;
	long i1,i2,j,j1,j2,N2;

	N2=num*num;

	for(i1=0;i1<N2;i1++)
		for(i2=0;i2<N2;i2++)
			B[i1][i2]=0.0;

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

			for(j1=0;j1<num;j1++)
				for(j2=0;j2<num;j2++)
				{
					j=j1*num+j2;

					N1[j]=NT1[i1][i2][j1][j2]*b11+NT2[i1][i2][j1][j2]*b12;		
				}
			
			add=b*W[i2]*W[i1];	
			for(j1=0;j1<N2;j1++)
				for(j2=0;j2<N2;j2++)
					B[j1][j2]+=(float)(N1[j1]*N1[j2]*add);		
		}
	}

}
////////////////////////////////////////Kx3 matrix//////////////////////////////////////////////////
void Kx3_matrix(short int num,short int nn,long node_numpml,long elem_num,long elem,float *lamda,float *mu,long *node3_total,
				long **Elem_array,double **node_array,float ****NT1,float ****NT2,long *IA,long *JA,float *Kx3)
{
	double *Stmp,*W,*Shape,*N1,cof,add;
	long i,j,k,*Jtmp,*id,*nd,*in,NP;
	float **B;
	short int N2;

	for(i=0;i<IA[node_numpml];i++)
		JA[i]=-1;

	N2=num*num;
	NP=num-1;

	Jtmp=Arraylong1(nn);
	Stmp=Arraydouble1(nn);
	nd=Arraylong1(N2);
	id=Arraylong1(N2);
	in=Arraylong1(N2);
	W=Arraydouble1(N2);
	Shape=Arraydouble1(N2);
	N1=Arraydouble1(N2);
	B=Arrayfloat2(N2,N2);

	Weight(NP,W);

	for(i=elem_num;i<elem;i++)
	{
		for(j=0;j<N2;j++)
		{
			nd[j]=Elem_array[i][j]-1;
			in[j]=node3_total[nd[j]]-1;
		}
		
		local_kx3(num,N2,node_array,nd,W,Shape,N1,NT1,NT2,B);

		cof=lamda[i]+2.0*mu[i];
		for(j=0;j<N2;j++)	
			for(k=0;k<N2;k++)
			{
				add=cof*B[j][k];
				Stiff_assemble(nn,in[j],in[k],add,Jtmp,Stmp,Kx3,IA,JA);
			}
	}
	free(Jtmp);
	free(Stmp);
	free(in);
	free(nd);
	free(id);
	free(N1);
	free(W);
	free(Shape);
	Freefloat2(N2,B);
} 
////////////////////////////////////////////////////////////////////////////////////////////////////
void local_kx3(short int num,short int N2,double **node_array,long *id,double *W,double *Shape,double *N1,
			   float ****NT1,float ****NT2,float **B)
{
	long i1,i2,j,j1,j2;
	double a11,a12,a21,a22,b,b11,b12,add;
	
	for(i1=0;i1<N2;i1++)
		for(i2=0;i2<N2;i2++)
			B[i1][i2]=0.0;
	
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

			for(j1=0;j1<num;j1++)
				for(j2=0;j2<num;j2++)
				{
					j=j1*num+j2;
					N1[j]=NT1[i1][i2][j1][j2]*b11+NT2[i1][i2][j1][j2]*b12;		
				}

			for(j=0;j<N2;j++)
				Shape[j]=0.0;

			j=i1*num+i2;
			Shape[j]=1.0;
					
			add=b*W[i2]*W[i1];
			for(j1=0;j1<N2;j1++)
				for(j2=0;j2<N2;j2++)
					 B[j1][j2]+=(float)(Shape[j1]*N1[j2]*add);	
		}
	}

}
////////////////////////////////////////Kx4 matrix//////////////////////////////////////////////////
void Kx4_matrix(short int num,short int nn,long node_numpml,long elem_num,long elem,float *mu,long *node3_total,long **Elem_array,
				double **node_array,float ****NT1,float ****NT2,long *IA,long *JA,float *Kx4)
{
	double *Stmp,*W,add,coeff,*N1;
	long i,i1,j1,j2,*Jtmp,*id,*id1,NP;
	float **B;
	short N2;

	for(i=0;i<IA[node_numpml];i++)
		JA[i]=-1;

	N2=num*num;
	id=Arraylong1(N2);
	id1=Arraylong1(N2);
	B=Arrayfloat2(N2,N2);
	Jtmp=Arraylong1(nn);
	Stmp=Arraydouble1(nn);
	N1=Arraydouble1(N2);

	W=Arraydouble1(num);

	NP=num-1;
	Weight(NP,W);

	for(i=elem_num;i<elem;i++)
	{

		for(i1=0;i1<N2;i1++)
		{
			id[i1]=Elem_array[i][i1]-1;
			id1[i1]=node3_total[id[i1]]-1;
		}

		local_Kx4(num,N2,node_array,id,N1,W,NT1,NT2,B);

		
		coeff=mu[i];
		for(j1=0;j1<N2;j1++)
			for(j2=0;j2<N2;j2++)
			{
				add=coeff*B[j1][j2];
				Stiff_assemble(nn,id1[j1],id1[j2],add,Jtmp,Stmp,Kx4,IA,JA);
			}
	}      
	free(Jtmp);
	free(Stmp);
	free(id);
	free(id1);
	free(W);
	free(N1);
	Freefloat2(N2,B);

} 
//////////////////////////////////////////local kx4//////////////////////////////////
void local_Kx4(short int num,short int N2,double **node_array,long *id,double *N1,double *W,float ****NT1,float ****NT2,float **B)
{
	double a11,a12,a21,a22,b,b21,b22,add;
	long i1,i2,j,j1,j2;


	for(i1=0;i1<N2;i1++)
		for(i2=0;i2<N2;i2++)
			B[i1][i2]=0.0;
	
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

			b21=-a21/b;
			b22=a11/b;

			for(j1=0;j1<num;j1++)
				for(j2=0;j2<num;j2++)
				{
					j=j1*num+j2;
					N1[j]=NT1[i1][i2][j1][j2]*b21+NT2[i1][i2][j1][j2]*b22;		
				}
				
			add=b*W[i2]*W[i1];
			for(j1=0;j1<N2;j1++)
				for(j2=0;j2<N2;j2++)
					B[j1][j2]+=(float)(N1[j1]*N1[j2]*add);
		}
	}

}

////////////////////////////////////////Kx6 matrix//////////////////////////////////////////////////
void Kx6_matrix(short int num,short int nn,long node_numpml,long elem_num,long elem,float *mu,long *node3_total,
				long **Elem_array,double **node_array,float ****NT1,float ****NT2,long *IA,long *JA,float *Kx6)
{
	long i,j,k,*Jtmp,*nd,*in,N2,NP;
	double *Stmp,*dx,*W,*Shape,*N1,add,coeff;
	float **B;

	for(i=0;i<IA[node_numpml];i++)
		JA[i]=-1;

	N2=num*num;
	NP=num-1;

	Jtmp=Arraylong1(nn);
	Stmp=Arraydouble1(nn);
	nd=Arraylong1(N2);
	in=Arraylong1(N2);
	dx=Arraydouble1(N2);
	W=Arraydouble1(N2);
	Shape=Arraydouble1(N2);
	N1=Arraydouble1(N2);
	B=Arrayfloat2(N2,N2);

	Weight(NP,W);

	for(i=elem_num;i<elem;i++)
	{
		for(j=0;j<N2;j++)
		{
			nd[j]=Elem_array[i][j]-1;
			in[j]=node3_total[nd[j]]-1;
		}

		local_kx6(num,node_array,nd,W,Shape,N1,NT1,NT2,B);

		coeff=mu[i];
		for(j=0;j<N2;j++)
			for(k=0;k<N2;k++)
			{
				add=B[j][k]*coeff;
				Stiff_assemble(nn,in[j],in[k],add,Jtmp,Stmp,Kx6,IA,JA);
			}
	}
	free(Jtmp);
	free(Stmp);
	free(in);
	free(nd);
	free(dx);
	free(N1);
	free(Shape);
	Freefloat2(N2,B);
} 
////////////////////////////////////////////////////////////////////////////////////////////////////
void local_kx6(short int num,double **node_array,long *id,double *W,double *Shape,double *N1,float ****NT1,float ****NT2,float **B)
{
	double a11,a12,a21,a22,b,b21,b22,add;
	long i1,i2,j,j1,j2,N2;
	
	N2=num*num;

	for(i1=0;i1<N2;i1++)
		for(i2=0;i2<N2;i2++)
			B[i1][i2]=0.0;

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

			b21=-a21/b;
			b22=a11/b;

			for(j1=0;j1<num;j1++)
				for(j2=0;j2<num;j2++)
				{
					j=j1*num+j2;	
					N1[j]=NT1[i1][i2][j1][j2]*b21+NT2[i1][i2][j1][j2]*b22;		
				}

			for(j=0;j<N2;j++)
				Shape[j]=0.0;

			j=i1*num+i2;
			Shape[j]=1.0;
			
			add=b*W[i2]*W[i1];		
			for(j1=0;j1<N2;j1++)
				for(j2=0;j2<N2;j2++)
					 B[j1][j2]+=(float)(Shape[j1]*N1[j2]*add);
		}
	}
}
/////////////////////////////////////////////Kx7 ////////////////////////////////////////////////////
void Kx7_matrix(short int num,short int nn,long node_numpml,long elem_num,long elem,float *mu,long *node3_total,
				long **Elem_array,double **node_array,float ****NT1,float ****NT2,long *IA,long *JA,float *Kx7)
{	
	double *Stmp,*W,*N1,add,coeff;
	long i,i1,j1,j2,*Jtmp,*id,*id1,NP,N2;
	float **B;

	for(i=0;i<IA[node_numpml];i++)
		JA[i]=-1;

	N2=num*num;
	id=Arraylong1(N2);
	id1=Arraylong1(N2);
	B=Arrayfloat2(N2,N2);
	Jtmp=Arraylong1(nn);
	Stmp=Arraydouble1(nn);
	W=Arraydouble1(num);
	N1=Arraydouble1(N2);

	NP=num-1;
	Weight(NP,W);

	for(i=elem_num;i<elem;i++)
	{
		for(i1=0;i1<N2;i1++)
		{
			id[i1]=Elem_array[i][i1]-1;
			id1[i1]=node3_total[id[i1]]-1;
		}

		local_Kx7(num,node_array,id,W,N1,NT1,NT2,B);
		coeff=mu[i];

		for(j1=0;j1<N2;j1++)
			for(j2=0;j2<N2;j2++)
			{
				add=B[j1][j2]*coeff;
				Stiff_assemble(nn,id1[j1],id1[j2],add,Jtmp,Stmp,Kx7,IA,JA);
			}
	}      
	free(Jtmp);
	free(Stmp);
	free(W);
	free(id);
	free(id1);
	free(N1);
	Freefloat2(N2,B);

} 
//////////////////////////////////////////local kx4//////////////////////////////////
void local_Kx7(short int num,double **node_array,long *id,double *W,double *N1,float ****NT1,float ****NT2,float **B)
{
	double a11,a12,a21,a22,b,b11,b12,add;
	long i1,i2,j,j1,j2,N2;

	N2=num*num;

	for(i1=0;i1<N2;i1++)
		for(i2=0;i2<N2;i2++)
			B[i1][i2]=0.0;

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

			for(j1=0;j1<num;j1++)
				for(j2=0;j2<num;j2++)
				{
					j=j1*num+j2;
					N1[j]=NT1[i1][i2][j1][j2]*b11+NT2[i1][i2][j1][j2]*b12;		
				}
				
			add=b*W[i2]*W[i1];
			for(j1=0;j1<N2;j1++)
				for(j2=0;j2<N2;j2++)
					B[j1][j2]+=(float)(N1[j1]*N1[j2]*add);
		}
	}

}
////////////////////////////////////////Kx8 matrix//////////////////////////////////////////////////
void Kx8_matrix(short int num,short int nn,long node_numpml,long elem_num,long elem,float *mu,long *node3_total,
				long **Elem_array,double **node_array,float ****NT1,float ****NT2,long *IA,long *JA,float *Kx8)
{
	double *Stmp,*W,add,coeff,*Shape,*N1;
	long i,j,k,*Jtmp,*nd,*in,N2,NP;
	float **B;

	for(i=0;i<IA[node_numpml];i++)
		JA[i]=-1;

	N2=num*num;
	NP=num-1;
	Jtmp=Arraylong1(nn);
	Stmp=Arraydouble1(nn);
	nd=Arraylong1(N2);
	in=Arraylong1(N2);
	W=Arraydouble1(N2);
	B=Arrayfloat2(N2,N2);
	N1=Arraydouble1(N2);
	Shape=Arraydouble1(N2);

	Weight(NP,W);

	for(i=elem_num;i<elem;i++)
	{
		for(j=0;j<N2;j++)
		{
			nd[j]=Elem_array[i][j]-1;
			in[j]=node3_total[nd[j]]-1;	
		}

		local_kx8(num,node_array,nd,W,Shape,N1,NT1,NT2,B);

		coeff=mu[i];
		for(j=0;j<N2;j++)	
			for(k=0;k<N2;k++)
			{
				add=B[j][k]*coeff;
				Stiff_assemble(nn,in[j],in[k],add,Jtmp,Stmp,Kx8,IA,JA);
			}
	}
	free(Jtmp);
	free(Stmp);
	free(in);
	free(nd);
	free(W);
	free(N1);
	free(Shape);
	Freefloat2(N2,B);
} 
////////////////////////////////////////////////////////////////////////////////////////////////////
void local_kx8(short int num,double **node_array,long *id,double *W,double *Shape,double *N1,
			   float ****NT1,float ****NT2,float **B)
{
	double a11,a12,a21,a22,b,b11,b12,add;
	long i1,i2,j,j1,j2,N2;	
	
	N2=num*num;

	for(i1=0;i1<N2;i1++)
		for(i2=0;i2<N2;i2++)
			B[i1][i2]=0.0;
	

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

			for(j1=0;j1<num;j1++)
				for(j2=0;j2<num;j2++)
				{
					j=j1*num+j2;
					N1[j]=NT1[i1][i2][j1][j2]*b11+NT2[i1][i2][j1][j2]*b12;		
				}

			for(j=0;j<N2;j++)
				Shape[j]=0.0;

			j=i1*num+i2;
			Shape[j]=1.0;
			
			add=b*W[i2]*W[i1];
			for(j1=0;j1<N2;j1++)
				for(j2=0;j2<N2;j2++)		
					 B[j1][j2]+=(float)(Shape[j1]*N1[j2]*add);
		}
	}
}
////////////////////////////////////////Kx9 matrix//////////////////////////////////////////////////
void Kx9_matrix(short int num,short int nn,long node_numpml,long elem_num,long elem,float *lamda,float *mu,long *node3_total,
				  long **Elem_array,double **node_array,float ****NT1,float ****NT2,long *IA,long *JA,float *Kx9)
{
	double *Stmp,*W,*N1,cof,add;
	long i,i1,j1,j2,*Jtmp,*id,*id1,NP,N2;
	float **B;

	for(i=0;i<IA[node_numpml];i++)
		JA[i]=-1;

	N2=num*num;
	id=Arraylong1(N2);
	id1=Arraylong1(N2);
	B=Arrayfloat2(N2,N2);
	Jtmp=Arraylong1(nn);
	Stmp=Arraydouble1(nn);
	W=Arraydouble1(num);
	N1=Arraydouble1(N2);

	NP=num-1;
	Weight(NP,W);

	for(i=elem_num;i<elem;i++)
	{
		for(i1=0;i1<N2;i1++)
		{
			id[i1]=Elem_array[i][i1]-1;
			id1[i1]=node3_total[id[i1]]-1;
		}

		local_Kx9(num,node_array,id,W,N1,NT1,NT2,B);

		cof=lamda[i]+2.0*mu[i];
		
		for(j1=0;j1<N2;j1++)
			for(j2=0;j2<N2;j2++)
			{
				add=B[j1][j2]*cof;
				Stiff_assemble(nn,id1[j1],id1[j2],add,Jtmp,Stmp,Kx9,IA,JA);
			}
	}      
	free(Jtmp);
	free(Stmp);
	free(id);
	free(id1);
	free(W);
	free(N1);
	Freefloat2(N2,B);
} 
//////////////////////////////////////////local kx1//////////////////////////////////
void local_Kx9(short int num,double **node_array,long *id,double *W,double *N1,float ****NT1,float ****NT2,float **B)
{
	double a11,a12,a21,a22,b,b21,b22,add;
	long i1,i2,j,j1,j2,N2;

	N2=num*num;

	for(i1=0;i1<N2;i1++)
		for(i2=0;i2<N2;i2++)
			B[i1][i2]=0.0;

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

			b21=-a21/b;
			b22=a11/b;

			for(j1=0;j1<num;j1++)
				for(j2=0;j2<num;j2++)
				{
					j=j1*num+j2;
					N1[j]=NT1[i1][i2][j1][j2]*b21+NT2[i1][i2][j1][j2]*b22;		
				}
	
			add=b*W[i2]*W[i1];
			for(j1=0;j1<N2;j1++)
				for(j2=0;j2<N2;j2++)
					B[j1][j2]+=(float)(N1[j1]*N1[j2]*add);
		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////
void Kx10_matrix(short int num,short int nn,long node_numpml,long elem_num,long elem,float *lamda,float *mu,long *node3_total,
				 long **Elem_array,double **node_array,float ****NT1,float ****NT2,long *IA,long *JA,float *Kx10)
{
	double *W,*Stmp,*Shape,*N1,cof,add;
	long i,j,k,*Jtmp,*nd,*in,N2,NP;
	float **B;
	
	for(i=0;i<IA[node_numpml];i++)
		JA[i]=-1;

	N2=num*num;
	NP=num-1;
	Jtmp=Arraylong1(nn);
	Stmp=Arraydouble1(nn);
	nd=Arraylong1(N2);
	in=Arraylong1(N2);
	W=Arraydouble1(N2);
	B=Arrayfloat2(N2,N2);
	Shape=Arraydouble1(N2);
	N1=Arraydouble1(N2);

	Weight(NP,W);

	for(i=elem_num;i<elem;i++)
	{
		for(j=0;j<N2;j++)
		{
			nd[j]=Elem_array[i][j]-1;
			in[j]=node3_total[nd[j]]-1;
		}
		
		local_kx10(num,node_array,nd,W,Shape,N1,NT1,NT2,B);

		cof=lamda[i]+2.0*mu[i];

		for(j=0;j<N2;j++)
			for(k=0;k<N2;k++)
			{
				add=B[j][k]*cof;
				Stiff_assemble(nn,in[j],in[k],add,Jtmp,Stmp,Kx10,IA,JA);
			}
	}
	free(Jtmp);
	free(Stmp);
	free(in);
	free(nd);
	free(N1);
	free(Shape);
	Freefloat2(N2,B);
} 
////////////////////////////////////////////////////////////////////////////////////////////////////
void local_kx10(short int num,double **node_array,long *id,double *W,double *Shape,double *N1,float ****NT1,
				float ****NT2,float **B)
{
	double a11,a12,a21,a22,b,b21,b22,add;
	long i1,i2,j,j1,j2,N2;
	
	N2=num*num;

	for(i1=0;i1<N2;i1++)
		for(i2=0;i2<N2;i2++)
			B[i1][i2]=0.0;

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

			b21=-a21/b;
			b22=a11/b;

			for(j1=0;j1<num;j1++)
				for(j2=0;j2<num;j2++)
				{
					j=j1*num+j2;
					N1[j]=NT1[i1][i2][j1][j2]*b21+NT2[i1][i2][j1][j2]*b22;		
				}

			for(j=0;j<N2;j++)
				Shape[j]=0.0;

			j=i1*num+i2;
			Shape[j]=1.0;
			
			add=b*W[i2]*W[i1];
			for(j1=0;j1<N2;j1++)
				for(j2=0;j2<N2;j2++)
					 B[j1][j2]+=(float)(Shape[j1]*N1[j2]*add);	
		}
	}
}
//////////////////////////////////////////kx11//////////////////////////////////
void Kx11_matrix(short int num,short int nn,long node_numpml,long elem_num,long elem,float *lamda,float *mu,long *node3_total,
				  long **Elem_array,double **node_array,float ****NT1,float ****NT2,long *IA,long *JA,float *Kx11)
{
	long i,i1,j1,j2,*Jtmp,*id,*id1,NP,N2;
	float **Bxz,**Bzx,add;
	double *Stmp,*W,cf1,cf2,*NN1,*NN2;
	
	for(i=0;i<IA[node_numpml];i++)
		JA[i]=-1;

	N2=num*num;
	id=Arraylong1(N2);
	id1=Arraylong1(N2);
	Bxz=Arrayfloat2(N2,N2);
	Bzx=Arrayfloat2(N2,N2);
	Jtmp=Arraylong1(nn);
	Stmp=Arraydouble1(nn);
	NN1=Arraydouble1(N2);
	NN2=Arraydouble1(N2);

	W=Arraydouble1(num);

	NP=num-1;
	Weight(NP,W);

	for(i=elem_num;i<elem;i++)
	{
		for(i1=0;i1<N2;i1++)
		{
			id[i1]=Elem_array[i][i1]-1;
			id1[i1]=node3_total[id[i1]]-1;
		}

		local_Kx11(num,node_array,id,W,NN1,NN2,NT1,NT2,Bxz,Bzx);
				
		cf1=lamda[i];
		cf2=mu[i];

		for(j1=0;j1<N2;j1++)
			for(j2=0;j2<N2;j2++)
			{
				add=(float)(Bxz[j1][j2]*cf1+Bzx[j1][j2]*cf2);
				Stiff_assemble(nn,id1[j1],id1[j2],add,Jtmp,Stmp,Kx11,IA,JA);
			}
	}      
	free(Jtmp);
	free(Stmp);
	free(id);
	free(id1);
	free(W);
	free(NN1);
	free(NN2);
	Freefloat2(N2,Bxz);
	Freefloat2(N2,Bzx);
} 
//////////////////////////////////////////local kx11//////////////////////////////////
void local_Kx11(short int num,double **node_array,long *id,double *W,double *NN1,double *NN2,float ****NT1,
				float ****NT2,float **Bxz,float **Bzx)
{
	double a11,a12,a21,a22,b,b11,b12,b21,b22,add;
	long i1,i2,j,j1,j2;
	short int N2;

	N2=num*num;

	for(i1=0;i1<N2;i1++)
		for(i2=0;i2<N2;i2++)
		{
			Bxz[i1][i2]=0.0;
			Bzx[i1][i2]=0.0;
		}

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

					NN1[j]=NT1[i1][i2][j1][j2]*b11+NT2[i1][i2][j1][j2]*b12;
					NN2[j]=NT1[i1][i2][j1][j2]*b21+NT2[i1][i2][j1][j2]*b22;
				}
			
			add=b*W[i2]*W[i1];
			for(j1=0;j1<N2;j1++)
				for(j2=0;j2<N2;j2++)
				{
					Bxz[j1][j2]+=(float)(NN1[j1]*NN2[j2]*add);
					Bzx[j1][j2]+=(float)(NN2[j1]*NN1[j2]*add);
				}
		}
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
void Kx13_matrix(short int num,short int nn,long node_numpml,long elem_num,long elem,float *lamda,float *mu,long *node3_total,
				  long **Elem_array,double **node_array,float ****NT1,float ****NT2,long *IA,long *JA,float *Kx13)
{
	long i,i1,j1,j2,*Jtmp,*id,*id1,NP;
	float **Bxz,**Bzx;
	double *Stmp,*W,add,*NN1,*NN2,cf1,cf2;
	short int N2;

	for(i=0;i<IA[node_numpml];i++)
		JA[i]=-1;

	N2=num*num;
	id=Arraylong1(N2);
	id1=Arraylong1(N2);
	Bxz=Arrayfloat2(N2,N2);
	Bzx=Arrayfloat2(N2,N2);
	Jtmp=Arraylong1(nn);
	Stmp=Arraydouble1(nn);
	NN1=Arraydouble1(N2);
	NN2=Arraydouble1(N2);

	W=Arraydouble1(num);

	NP=num-1;
	Weight(NP,W);

	for(i=elem_num;i<elem;i++)
	{
		for(i1=0;i1<N2;i1++)
		{
			id[i1]=Elem_array[i][i1]-1;
			id1[i1]=node3_total[id[i1]]-1;
		}
	
		local_Kx11(num,node_array,id,W,NN1,NN2,NT1,NT2,Bxz,Bzx);

		cf1=lamda[i];
		cf2=mu[i];
		
		for(j1=0;j1<N2;j1++)
			for(j2=0;j2<N2;j2++)
			{
				add=(float)(Bzx[j1][j2]*cf1+Bxz[j1][j2]*cf2);
				Stiff_assemble(nn,id1[j1],id1[j2],add,Jtmp,Stmp,Kx13,IA,JA);
			}
	}      
	free(Jtmp);
	free(Stmp);
	free(W);
	free(id);
	free(id1);
	free(NN1);
	free(NN2);
	Freefloat2(N2,Bxz);
	Freefloat2(N2,Bzx);
} 
