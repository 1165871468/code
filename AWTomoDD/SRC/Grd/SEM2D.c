#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "SEM.h"
#define Min 1.0e-10
#define PI 3.141592653589793238462643

int main(int argc, char *argv[])
{
	clock_t start, finish;

	char *source, *receiver, *mesh, *moment, flnm1[100];
	short N, sflag, stype, nstep, nsnap, max, Pflag, n1, ndiv, d1, d2, d3, d4;
	short *elem_par, *Count_array, *Count_array_add, *Cin;
	int rank, size, s1, s2, count1, numData, numNs, numNr;

	long node, elem, tmp, i, i1, i2, j, k, edge_num, node_add, Sii, Nr, Ns;	
	long *Iedge, *Jedge, *Edge, *Istiff, *Jstiff, *node1_total, *node1_array, node_num, Nt;
	long elem_num, *Elem1_total, *Elem2_total, elem_num1, pnode_num, pelem_num;
	long *Elem1_count, *Elem2_count, *node2_total, *node2_array, *node3_total, bnode_num;
	long node_numpmlc, *node3_array, *node4_array, *node4_total, node_numpml;
	long *IA, *JA, *Istiff_ini, *Jstiff_ini, *IA_ini, *JA_ini, *Rii, cnode_num, **Elem_array;

	float *U, *U1, *Result, *lamda, *mu, *rho, *Mass, *Mass_tmp, *Loadx, *Loadz, *Stiff, *Ux1, *Ux2, *Ux3, *Ux1t, **Loadx_s, **Loadz_s;
	float *Ux2t, *Ux3t, *Vx1, *Vx2, *Vx3, *Uz1, *Uz2, *Uz3, *Uz1t, *Uz2t, *Uz3t, *Vz1, *Vz2, *Vz3;
	float *Kx1, *Kx3, *Kx4, *Kx6, *Kx7, *Kx8, *Kx9, *Kx10, *Kx11, *Kx13;
	float **KKK, **Kxx, **Kxz, **Kzx, **Kzz, **PPP, **Pxx, **Pxz, **Pzx, **Pzz;

	double *dx, *dz, *Rax, *Raz, **NT, **ND, **Tensor, **Rloc, **Sour, **node_arrayt, **node_array, **arrayData;
	double duration, deg, X, Z, Sx, Sz, F, T, dt, Layx, Layz, Vp, Vs, dis, Sax, Saz, len, R0, xx, zz, Vp1, Vs1;

	start=clock();

	mesh = Arraychar1(100);//
	source = Arraychar1(100);//
	receiver = Arraychar1(100);//
	moment = Arraychar1(100);//

	Read_par(&N,&Pflag,&ndiv,&nstep,&nsnap,&sflag,&stype,&Ns,&Nr,&X,&Z,&T,&dt,&F,&Vp,&Layx,&Layz,&R0,&dis,source,receiver,mesh,moment);

	N = N+1;   len = dis;   dis = 0.002*dis;  Nt = (int)((T+Min)/dt);
	node = 0;  elem = 0;    n1 = 0;
	//////////////////////////////////////////////////////////
	sprintf(flnm1, "../Input/arrival_p.txt");
	numData = count_file(flnm1);
	arrayData = Arraydouble2(numData, 3);
	rdfile(flnm1, numData, 3, arrayData);

	numNs = Ns;
	numNr = Nr;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	Average_rank(rank, size, numData, &s1, &s2);
	//////////////////////////////////////////////////////////
	Mesh(&node,&elem,mesh);

	node_arrayt=Arraydouble2(node,2);//
	Elem_array=Arraylong2(elem,N*N);//
	elem_par=Arrayshort1(elem);//

	Get_par(node,elem,node_arrayt,Elem_array,elem_par,mesh);
	Count_array=Arrayshort1(node);
	Count(node,elem,Elem_array,Count_array,&n1);

	Iedge=Arraylong1(node+1);//
	tmp=0;
	for(i=0;i<node;i++)
	{
		tmp+=Count_array[i]+1;				//缓冲
		Iedge[i+1]=tmp;					//所有节点使用总和Iedge[node]=118323
	}
	Jedge=Arraylong1(tmp);//
	Edge=Arraylong1(node);//

	Count_edge(node,elem,n1,Elem_array,Edge,Iedge,Jedge,&edge_num);	//数边edge_num=47549并计算插值点

	node_add=node+edge_num*(N-2)+elem*(N-2)*(N-2);			//591611插值后点的个数
	node_array=Arraydouble2(node_add,2);//

	re_node_elem(N,node,elem,edge_num,Elem_array,node_arrayt,node_array,Edge,Iedge,Jedge);//插值并更新节点信息
	
	Count_array_add=Arrayshort1(node_add);//
	Count_node_csr(N,elem,Elem_array,Count_array,Count_array_add);	//计算csr节点；csr是一种矩阵的高效存储方式

	Resort_node_element(N,elem,Elem_array);				//插值后重新记录元素单元

	lamda = Arrayfloat1(elem);//
	mu = Arrayfloat1(elem);//
	rho = Arrayfloat1(elem);//

//	Output_vel(N, node, elem, Elem_array, rho, lamda, mu, node_array);  printf("Vel yes!\n"); sleep(1000);
	Input_vel(elem, rho, lamda, mu);

	ND=Arraydouble2(N,N);//
	NT=Arraydouble2(N,N);//
	Common_array(N,ND,NT);
	//////////////////////////////////Mass matrix///////////////////////
	Mass = Arrayfloat1(node_add);
	Mass_tmp = Arrayfloat1(node_add);
	Mass_matrix(N,elem,rho,node_array,Elem_array,Mass,Mass_tmp);

	for(i = 0; i < node_add; i++)
	{
		Mass[i] = 1.0 / Mass[i];
	}

	/////////////////////////////////Stiff matrix/////////////////////////
	node1_total=Arraylong1(node_add);//
	node2_total=Arraylong1(node_add);//

	Count_material_num(N,Pflag,node_add,elem,Elem_array,elem_par,node1_total,&node_num,&elem_num);

	node1_array=Arraylong1(node_num);//
	Elem1_count=Arraylong1(elem_num);//
	node_sort(node_add,node_num,node1_total,node1_array);
	Elem_computation(Pflag,elem,elem_par,Elem1_count);

	KKK=Arrayfloat2(elem_num,N*N);//
	Kxx=Arrayfloat2(elem_num,N*N);//
	Kxz=Arrayfloat2(elem_num,N*N);//
	Kzx=Arrayfloat2(elem_num,N*N);//
	Kzz=Arrayfloat2(elem_num,N*N);//
	Stiff_matrix(N,elem_num,Elem_array,Elem1_count,node_array,KKK,Kxx,Kxz,Kzx,Kzz);

	//////////////////////////////////////PML matrix////////////////////////////////
	Count_material_pml(N,Pflag,node_add,elem,Elem_array,elem_par,node2_total,&pnode_num,&pelem_num);

	node2_array=Arraylong1(pnode_num);//
	Elem2_count=Arraylong1(pelem_num);//

	node_sort(node_add,pnode_num,node2_total,node2_array);
	Elem_pml(Pflag,elem,elem_par,Elem2_count);

	PPP=Arrayfloat2(pelem_num,N*N);//
	Pxx=Arrayfloat2(pelem_num,N*N);//
	Pxz=Arrayfloat2(pelem_num,N*N);//
	Pzx=Arrayfloat2(pelem_num,N*N);//
	Pzz=Arrayfloat2(pelem_num,N*N);//

	Stiff_matrix(N,pelem_num,Elem_array,Elem2_count,node_array,PPP,Pxx,Pxz,Pzx,Pzz);

	//////////////////////////////////////Load matrix////////////////////////////
	Loadx=Arrayfloat1(node_num);//
	Loadz=Arrayfloat1(node_num);//
	Loadx_s=Arrayfloat2(Nr, node_num);//
	Loadz_s=Arrayfloat2(Nr, node_num);//
	Tensor=Arraydouble2(2,2);//
	Sour=Arraydouble2(Ns,2);//

	Read_tensor(Tensor);
	Read_location2(Ns,Sour,source);
	////////////////////////////////////Find-Station//////////////////////////////
	Rii=Arraylong1(Nr);//
	Rax=Arraydouble1(Nr);//
	Raz=Arraydouble1(Nr);//
	Rloc=Arraydouble2(Nr,2);//

	Read_location2(Nr,Rloc,receiver);

	//////////////////////////////////Collection//////////////////////////////////
	node3_total=Arraylong1(node_add);//
	node4_total=Arraylong1(node_add);//

	Collection(node_add,node1_total,node2_total,node3_total,&cnode_num);

	node3_array=Arraylong1(cnode_num);//

	node_sort(node_add,cnode_num,node3_total,node3_array);

	////////////////////////////////Outside-Boundary//////////////////////////////
	Outside_boundary(node_add,X,Z,Layx,Layz,dis,node_array,&bnode_num,node4_total);

	node4_array=Arraylong1(bnode_num);//

	node_sort(node_add,bnode_num,node4_total,node4_array);

	for(i = s1; i < s2; i++) {
		if (rank == 0) {
			deg = (i + 1.0)/s2*100.0;
			printf("%10.1f%%\n", deg);
		}

		Ns = arrayData[i][0] - 1;
		Nr = arrayData[i][1] - 1;
		numNs = Nr;

		Sii = -2; Sax = -2.0; Saz = -2.0;
		Find_local_position(N,ndiv,elem_num,len,Rloc[numNs][0],Rloc[numNs][1],Elem1_count,Elem_array,node_array,&Sii,&Sax,&Saz);
		Sx = Rloc[numNs][0]; Sz = Rloc[numNs][1];

		for (j = 0; j < node_num; j++) {
			Loadx[j] = 0.0, Loadz[j] = 0.0;
		}

		if (Sii != -2) {
			Loadf_sem2(N, elem_num, len, Sx, Sz, node_array, node1_total, Elem_array, Elem1_count, Loadx, Loadz, Tensor);

			for(j = 0; j < node_num; j++) {
				Loadx_s[numNs][j] += Loadx[j];
				Loadz_s[numNs][j] += Loadz[j];
			}
		}
		// printf("cc0 %ld %ld %d\n", Ns, Nr, numNs), sleep(1000);
		Update(ndiv,elem,len,N,nstep,nsnap,sflag,Ns,Nr,Nt,node_num,elem_num,pnode_num,pelem_num,node_add,cnode_num,bnode_num,F,X,Z,Sx,Sz,dt,R0,Vp,Layx,Layz,Rii,Elem1_count,Elem2_count,
		node1_total,node1_array,node2_total,node2_array,node3_total,node3_array,node4_array,Elem_array,node_array,NT,lamda,mu,Mass,Loadx,Loadz,KKK,
		Kxx,Kxz,Kzx,Kzz,PPP,Pxx,Pxz,Pzx,Pzz,Rax,Raz,Loadx_s,Loadz_s,Rloc);
	}

	//////////////////////////////////// Free /////////////////////////////////////
	free(mesh);        free(source);        free(receiver);        free(moment);
	free(Edge);        free(Iedge);         free(Jedge);           free(elem_par);
	free(Count_array); free(lamda);         free(mu);              free(rho);
	free(Mass);        free(Mass_tmp);      free(node1_total);     free(node2_total);
	free(node3_total); free(node4_total);   free(node1_array);     free(node2_array);
	free(node3_array); free(node4_array);   free(Elem1_count);     free(Elem2_count);
	free(Loadx);       free(Loadz);         free(Rii);             free(Rax);
	free(Raz);         

	Freelong2(elem, Elem_array);

	Freefloat2(elem_num, KKK);              Freefloat2(elem_num, Kxx);
	Freefloat2(elem_num, Kxz);              Freefloat2(elem_num, Kzx);
	Freefloat2(elem_num, Kzz);              Freefloat2(pelem_num, PPP);
	Freefloat2(pelem_num, Pxx);             Freefloat2(pelem_num, Pxz);
	Freefloat2(pelem_num, Pzx);             Freefloat2(pelem_num, Pzz);
	Freefloat2(Nr, Loadx_s);                Freefloat2(Nr, Loadz_s);

	Freedouble2(node_add, node_array);      Freedouble2(node, node_arrayt);
	Freedouble2(N, ND);                     Freedouble2(N, NT);
	Freedouble2(2, Tensor);                 Freedouble2(Ns, Sour);
	Freedouble2(Nr, Rloc);                  Freedouble2(numData, arrayData);


	MPI_Finalize();

	finish=clock();
	duration=(double)(finish-start)/CLOCKS_PER_SEC;
	// printf("t=%f seconds\n",duration);
	// printf("Yes Grident\n");
	return 0;
}
