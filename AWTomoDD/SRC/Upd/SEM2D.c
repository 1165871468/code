#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "SEM.h"
#include "Common.h"
#define Min 1.0e-10
#define PI 3.141592653589793238462643

int main(int argc,char *argv[])
{
	clock_t start,finish;

	int rank, size, s1, s2, nr, ns, numData, numDt, numIntervalt, flg, fs;
	short N,sflag,stype,nstep,nsnap,max,Pflag,*elem_par,*Count_array,n1,ndiv,*Count_array_add;
	short d1,d2,d3,d4,*Cin;
	long node,elem,tmp,i,i1,i2,j,k,edge_num,node_add,Rii,Sii,Nr,Ns,Riis[2],**Elem_array;	
	long *Iedge,*Jedge,*Edge,*Istiff,*Jstiff,*node1_total,*node1_array,node_num,Nt,numUu;
	long elem_num,*Elem1_total,*Elem2_total,elem_num1,pnode_num,pelem_num;
	long *Elem1_count,*Elem2_count,*node2_total,*node2_array,*node3_total,bnode_num;
	long node_numpmlc,*node3_array,*node4_array,*node4_total,node_numpml;
	long *IA,*JA,*Istiff_ini,*Jstiff_ini,*IA_ini,*JA_ini,cnode_num;
	float *U,*U1,*Result,*lamda,*mu,*rho,*Mass,*Mass_tmp,*Loadx,*Loadz,*Stiff,*Ux1,*Ux2,*Ux3,*Ux1t,**Loadx_s,**Loadz_s;
	float *Ux2t,*Ux3t,*Vx1,*Vx2,*Vx3,*Uz1,*Uz2,*Uz3,*Uz1t,*Uz2t,*Uz3t,*Vz1,*Vz2,*Vz3;
	float *Kx1,*Kx3,*Kx4,*Kx6,*Kx7,*Kx8,*Kx9,*Kx10,*Kx11,*Kx13;
	float **KKK,**Kxx,**Kxz,**Kzx,**Kzz,**PPP,**Pxx,**Pxz,**Pzx,**Pzz;
	float *arrayUadj, *arrayUsyn;
	double duration, X, Z, Rx, Rz, Sx, Sz, F, T, dt, Layx, Layz, Vp, Vs, dis, Rax, Raz, Sax, Saz, len, R0;
	double Raxs[2], Razs[2], *dx, *dz, **timepDatas, **Rlocs, **NT, **ND, **Tensor, **Sour, **node_arrayt, **node_array;
	double TmpD1, TmpD2;
	char *source,*receiver,*mesh,*moment,ff1[100];
	start = clock();
	////// Locd model parameters ///
	mesh = Arraychar1(100);
	source = Arraychar1(100);
	receiver = Arraychar1(100);
	moment = Arraychar1(100);

	numIntervalt = 5;
    Read_parForwar(&N, &Pflag, &ndiv, &nstep, &nsnap, &sflag, &stype, &dt, &F,
    &Vp, &Layx, &Layz, &R0, &dis, &X, &Z, moment);

	sprintf(mesh, "../Input/Mesh_test.txt");

	N = N + 1;              len = dis;         dis = 0.002*dis;
	Layx = Layx*len;        Layz = Layz*len;
	//////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////
	numData = count_file("../Input/arrival_p.txt");
	timepDatas = Arraydouble2(numData, 3);
	rdfile("../Input/arrival_p.txt", numData, 3, timepDatas);

	/// Mpi start! ! ! ///
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	Average_rank(rank, size, numData, &s1, &s2);
	//////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////
	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == 0) {
		printf("\nFoward siesmic wavefileds: \n");

		FILE *fp = fopen("../Output/Error/Error.txt", "w"); fclose(fp);
	}

	for (k = s1; k < s2; k++) {

		Ns = timepDatas[k][0];
		Nr = timepDatas[k][1];
		T = timepDatas[k][2];

		loadRecAndSor(Nr, Ns, &Rx, &Rz, &Sx, &Sz);
		//////////////////////////////////////////////////////////
		Nt = round(T / dt);
		
		node = 0;            elem = 0;          n1 = 0;
		Mesh(&node, &elem, mesh);
		Count_array = Arrayshort1(node);
		elem_par = Arrayshort1(elem);
		Iedge = Arraylong1(node+1);
		Edge = Arraylong1(node);
		Elem_array = Arraylong2(elem, N*N);
		node_arrayt = Arraydouble2(node, 2);
		Get_par(node,elem,node_arrayt,Elem_array,elem_par,mesh);
		Count(node,elem,Elem_array,Count_array,&n1);
		tmp = 0;

		for(i = 0; i < node; i++) {
			tmp += Count_array[i] + 1;
			Iedge[i+1] = tmp;
		}

		Jedge = Arraylong1(tmp);
		Count_edge(node,elem,n1,Elem_array,Edge,Iedge,Jedge,&edge_num);
		node_add = node+edge_num*(N-2)+elem*(N-2)*(N-2);
		node_array = Arraydouble2(node_add,2);
		re_node_elem(N,node,elem,edge_num,Elem_array,node_arrayt,node_array,Edge,Iedge,Jedge);
		Count_array_add = Arrayshort1(node_add);
		Count_node_csr(N, elem, Elem_array, Count_array,Count_array_add);
		Resort_node_element(N,elem,Elem_array);
		lamda = Arrayfloat1(elem);
		mu = Arrayfloat1(elem);
		rho = Arrayfloat1(elem);

		// Output_vel1(N, node, elem, Elem_array, rho, lamda, mu, node_array);  printf("Vel yes!\n"); sleep(1000);
		Input_vel(elem, rho, lamda, mu);

		ND = Arraydouble2(N,N);
		NT = Arraydouble2(N,N);
		Common_array(N,ND,NT);
		//////////////////////////////////Mass matrix///////////////////////
		Mass = Arrayfloat1(node_add);
		Mass_tmp = Arrayfloat1(node_add);
		Mass_matrix(N,elem,rho,node_array,Elem_array,Mass,Mass_tmp);
		/////////////////////////////////Stiff matrix/////////////////////////
		node1_total=Arraylong1(node_add);
		Count_material_num(N,Pflag,node_add,elem,Elem_array,elem_par,node1_total,&node_num,&elem_num);
		node1_array=Arraylong1(node_num);
		node_sort(node_add,node_num,node1_total,node1_array);
		Elem1_count=Arraylong1(elem_num);
		Elem_computation(Pflag,elem,elem_par,Elem1_count);
		KKK=Arrayfloat2(elem_num,N*N);
		Kxx=Arrayfloat2(elem_num,N*N);
		Kxz=Arrayfloat2(elem_num,N*N);
		Kzx=Arrayfloat2(elem_num,N*N);
		Kzz=Arrayfloat2(elem_num,N*N);
		Stiff_matrix(N,elem_num,Elem_array,Elem1_count,node_array,KKK,Kxx,Kxz,Kzx,Kzz);
		////// Locd PML matrix ///
		node2_total=Arraylong1(node_add);
		Count_material_pml(N,Pflag,node_add,elem,Elem_array,elem_par,node2_total,&pnode_num,&pelem_num);
		node2_array=Arraylong1(pnode_num);
		Elem2_count=Arraylong1(pelem_num);
		node_sort(node_add,pnode_num,node2_total,node2_array);
		Elem_pml(Pflag,elem,elem_par,Elem2_count);
		PPP=Arrayfloat2(pelem_num,N*N);
		Pxx=Arrayfloat2(pelem_num,N*N);
		Pxz=Arrayfloat2(pelem_num,N*N);
		Pzx=Arrayfloat2(pelem_num,N*N);
		Pzz=Arrayfloat2(pelem_num,N*N);
		Stiff_matrix(N,pelem_num,Elem_array,Elem2_count,node_array,PPP,Pxx,Pxz,Pzx,Pzz);
		for(i = 0; i < node_add; i++)
		{
			Mass[i] = 1.0 / Mass[i];
			Mass_tmp[i] = 1.0 / Mass_tmp[i];
		}
		////// Collection ///
		node3_total=Arraylong1(node_add);
		Collection(node_add,node1_total,node2_total,node3_total,&cnode_num);
		node3_array=Arraylong1(cnode_num);
		node_sort(node_add,cnode_num,node3_total,node3_array);
		////// Load Outside-Boundary ///
		node4_total=Arraylong1(node_add);
		Outside_boundary(node_add,X,Z,Layx,Layz,dis,node_array,&bnode_num,node4_total);
		node4_array=Arraylong1(bnode_num);
		node_sort(node_add,bnode_num,node4_total,node4_array);
		////// Load locations of Source and Reciver ///
		Loadx = Arrayfloat1(node_num);
		Loadz = Arrayfloat1(node_num);
		Loadx_s = Arrayfloat2(1, node_num);
		Loadz_s = Arrayfloat2(1, node_num);
		Rlocs = Arraydouble2(1, 2);
		Tensor = Arraydouble2(2,2);
		Read_tensor(Tensor);
		Sii = -2;   Sax = -2.0;  Sax = -2.0;
		Rii = -2;   Rax = -2.0;  Rax = -2.0;
		Find_local_position(N,ndiv,elem_num,len,Sx,Sz,Elem1_count,Elem_array,node_array,&Sii,&Sax,&Saz);
		Find_local_position(N,ndiv,elem_num,len,Rx,Rz,Elem1_count,Elem_array,node_array,&Rii,&Rax,&Raz);
		if (Sii == -2)
		{
			printf("Can not find the source location in the model! -- %.1f %.1f\n", Sx, Sz); return 1;
		}
		if (Rii == -2)
		{
			printf("Can not find the reciver location in the model! -- %.1f %.1f\n", Rx, Rz); return 1;
		}
		////// Load Source and Reciver location ///
		Riis[0] = Rii;
		Raxs[0] = Rax;
		Razs[0] = Raz;

		Loadf_sem_syn(N, elem_num, len, Sx, Sz, node_array, node1_total, Elem_array, Elem1_count, Loadx, Loadz, Tensor);
		for(j = 0; j < node_num; j++)
		{
			Loadx_s[0][j] = Loadx[j];
			Loadz_s[0][j] = Loadz[j];
		}

		///
		numDt = Nt / numIntervalt;
		numUu = 4*node_add*(Nt / numIntervalt);
		arrayUadj = Arrayfloat1(numUu);
		arrayUsyn = Arrayfloat1(numUu);
		///
		////// Calculate the synthetic wavefiled ///
		UpdateSyn(ndiv,elem,len,N,nstep,nsnap,sflag,Nr,Ns,Nt,node_num,elem_num,pnode_num,pelem_num,node_add,cnode_num,bnode_num,F,X,Z,Sx,Sz,dt,R0,Vp,Layx,Layz,Riis,Elem1_count,Elem2_count,
		node1_total,node1_array,node2_total,node2_array,node3_total,node3_array,node4_array,Elem_array,node_array,NT,lamda,mu,Mass,Mass_tmp,Loadx,Loadz,arrayUsyn,KKK,
		Kxx,Kxz,Kzx,Kzz,PPP,Pxx,Pxz,Pzx,Pzz,Raxs,Razs,i,Loadx_s,Loadz_s,Rlocs);
		//////
		////// Calculate the adjoint Source ///
		// getAdjSource(Ns, Nr, dt, T);
		GetSourceTimeAdj(Ns, Nr, dt, T);
		//////
		////// Calculate the adjoint wavefiled ///
		Sii = -2;    Sax = -2.0;  Sax = -2.0;
		Rii = -2;    Rax = -2.0;  Rax = -2.0;
		TmpD1 = Sx;  TmpD2 = Sz;
		Sx = Rx;     Sz = Rz;
		Rx = TmpD1;  Rz = TmpD2;
		Find_local_position(N,ndiv,elem_num,len,Sx,Sz,Elem1_count,Elem_array,node_array,&Sii,&Sax,&Saz);
		Find_local_position(N,ndiv,elem_num,len,Rx,Rz,Elem1_count,Elem_array,node_array,&Rii,&Rax,&Raz);
		if (Sii == -2.0)
		{
			printf("Can not find the source location in the model! -- %.1f %.1f\n", Sx, Sz); return 1;
		}
		if (Rii == -2.0)
		{
			printf("Can not find the reciver location in the model! -- %.1f %.1f\n", Rx, Rz); return 1;
		}

		for(j = 0; j < node_num; j++)
		{
			Loadx[j] = 0.0;
			Loadz[j] = 0.0;
		}
		////// Load Source and Reciver location ///
		Riis[0] = Rii;
		Raxs[0] = Rax;
		Razs[0] = Raz;
		///
		Loadf_sem_adj(N, elem_num, len, Sx, Sz, node_array, node1_total, Elem_array, Elem1_count, Loadx, Loadz, Tensor);

		for(j = 0; j < node_num; j++)
		{
			Loadx_s[0][j] = Loadx[j];
			Loadz_s[0][j] = Loadz[j];
		}

		UpdateAdj(ndiv,elem,len,N,nstep,nsnap,sflag,Nr,Ns,Nt,node_num,elem_num,pnode_num,pelem_num,node_add,cnode_num,bnode_num,F,X,Z,Sx,Sz,dt,R0,Vp,Layx,Layz,Riis,Elem1_count,Elem2_count,
		node1_total,node1_array,node2_total,node2_array,node3_total,node3_array,node4_array,Elem_array,node_array,NT,lamda,mu,Mass,Mass_tmp,Loadx,Loadz,arrayUadj,KKK,
		Kxx,Kxz,Kzx,Kzz,PPP,Pxx,Pxz,Pzx,Pzz,Raxs,Razs,i,Loadx_s,Loadz_s,Rlocs);

		getKernel(Ns, Nr, numIntervalt, numDt, node_add, numUu, dt, arrayUsyn, arrayUadj);
		//////
		////// Free arrays ///
		free(Edge);        free(Iedge);         free(Jedge);           free(elem_par);
		free(Count_array); free(lamda);         free(mu);              free(rho);
		free(Mass);        free(Mass_tmp);      free(node1_total);     free(node2_total);
		free(node3_total); free(node4_total);   free(node1_array);     free(node2_array);
		free(node3_array); free(node4_array);   free(Elem1_count);     free(Elem2_count);
		free(Loadx);       free(Loadz);         free(arrayUadj);       free(arrayUsyn);
		///
		Freelong2(elem, Elem_array);
		Freefloat2(elem_num, KKK);              Freefloat2(elem_num, Kxx);
		Freefloat2(elem_num, Kxz);              Freefloat2(elem_num, Kzx);
		Freefloat2(elem_num, Kzz);              Freefloat2(pelem_num, PPP);
		Freefloat2(pelem_num, Pxx);             Freefloat2(pelem_num, Pxz);
		Freefloat2(pelem_num, Pzx);             Freefloat2(pelem_num, Pzz);
		Freefloat2(1, Loadx_s);                 Freefloat2(1, Loadz_s);
		///
		Freedouble2(node_add, node_array);      Freedouble2(node, node_arrayt);
		Freedouble2(N, ND);                     Freedouble2(N, NT);
		Freedouble2(2, Tensor);                 Freedouble2(1, Rlocs);
		//////
		if (rank == 0 && size > 1)
		{
			printf("%5s %10.1f%% ...\n", " ", 100.0*((k + 1.0)/s2));
		}
        else if (size == 1)
        {
            printf("%5s %10.1f%% ...\n", " ", 100.0*((k + 1.0)/s2));
        }
	}
	/// Free arrays ///
	free(mesh);        free(source);
	free(receiver);    free(moment);
	Freedouble2(numData, timepDatas);
	/// Mpi end ///
	MPI_Finalize();
	//////
	return 0;
}
