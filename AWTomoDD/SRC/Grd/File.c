#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include "SEM.h"

///////////////////////////////// Load source time function ///////////////////////////
void GetRi(int Ns, int Nr, int Nt, int node_num, float **loadsemx, float **loadsemz, float **Rix_sem, float **Riz_sem)
{
	int  dt, fs, length1, length2, flgi1, flgi2, flgi3;
	int i, j, k, l;
	char file1[100];
	double a[3], b[2];
	FILE *fp1;

	sprintf(file1, "../Output/Source_adj/Source_times%d_r%d.txt", Ns+1, Nr+1);

	length1 = count_file(file1);

	fp1 = fopen(file1, "r");
	if(NULL == fp1)
	{
		printf("Open file failed! --- U%d_df.txt\n", Nr+1), exit(1);
	}

	if (length1 != Nt)
	{
		printf("Receiver data is error! -- %d Vs %d\n", length1, Nt), exit(1);
	}

	for (j = 0; j < Nt; j++)
	{
		a[0] = -1.0, a[1] = -1.0, a[2] = -1.0;

		fs = fscanf(fp1, "%lf%lf%lf ", &a[0], &a[1], &a[2]);

		for (k = 0; k < node_num; k++)
		{
			if(loadsemx[Nr][k] != 0.0 && a[1] != 0.0)
			{
				Rix_sem[k][j] += loadsemx[Nr][k]*a[1];
			}

			if(loadsemz[Nr][k] != 0.0 && a[1] != 0.0)
			{
				Riz_sem[k][j] += loadsemz[Nr][k]*a[1];
			}
		}
	}
	fclose(fp1);
}

void Input_vel(long elem, float *rho, float *lamda, float *mu)
{
	int i, j, fs;
	char flnm1[100];
	double xx, zz, Vp, Vs, P;
	FILE *fp;

	sprintf(flnm1, "../Output/Model/Forver/vel.dat");

	if ((fp = fopen(flnm1, "r")) == NULL)
	{
		printf("Read filed! -- %s\n", flnm1), exit(1);
	}

	for (i = 0; i < elem; i++)
	{
		fs = fscanf(fp, "%lf%lf%lf%lf%lf", &xx, &zz, &Vp, &Vs, &P);

		rho[i] = P;
		mu[i] = Vs*Vs*P;
		lamda[i] = Vp*Vp*P - 2.0*mu[i];
	}

	fclose(fp);
}

void Output_vel(short num, long node, long elem, long **Elem_array, float *rho, float *lamda, float *mu, double **node_array) {
	int i, j, k, i1, i2, i3, i4, Cin[4];
	float sx, sz, fx, fz, xx1, xx2, zz1, zz2, sx1, sx2, sx3, sx4, sz1, sz2, sz3, sz4, Vp1, Vp2, Vs1, Vs2, P1, P2, per;

	FILE *fp;

	per = 0.00;

	xx1 = 46000; xx2 = 54000;
	zz1 = 16000; zz2 = 24000;

	Vp1 = 4000.0; Vs1 = 2300.0; P1 = 2600.0;

	fp = fopen("../Initial/Forver/vel_initial.dat", "w");

	if (fp == NULL)
	{
		printf("can't open this file --- lam.txt\n"), exit(1);
	}
	
	for(i = 0; i < elem; i++)
	{
		Cin[0] = 0;
		Cin[1] = num * (num - 1);
		Cin[2] = num * num - 1;
		Cin[3] = num - 1;

		i1 = Elem_array[i][Cin[0]] - 1, i2 = Elem_array[i][Cin[1]] - 1;
		i3 = Elem_array[i][Cin[2]] - 1, i4 = Elem_array[i][Cin[3]] - 1;

		sx1 = node_array[i1][0],   sz1 = node_array[i1][1];
		sx2 = node_array[i2][0],   sz2 = node_array[i2][1];
		sx3 = node_array[i3][0],   sz3 = node_array[i3][1];
		sx4 = node_array[i4][0],   sz4 = node_array[i4][1];

		sx = (sx1 + sx2 + sx3 + sx4) / 4.0;
		sz = (sz1 + sz2 + sz3 + sz4) / 4.0;

		if (sx >= xx1 && sx <= xx2 && sz >= zz1 && sz <= zz2)
		{
			Vp2 = Vp1*(1.0 + per);
			Vs2 = Vs1;
			P2 = P1;
		}
		else
		{
			Vp2 = Vp1;
			Vs2 = Vs1;
			P2 = P1;
		}

		fprintf(fp, "%.1f%10.1f%10.1f%10.1f%10.1f\n", sx, sz, Vp2, Vs2, P2), fflush(fp);
	}

	fclose(fp);
}

void Getvel(long elem,short int *elem_par,float *rho,float *lamda,float *mu)
{
	double Vp, Vs, P;
	long i;
	
	for(i=0;i<elem;i++)
	{

		Vp = 5000.0;
		Vs = 2800.0;
		P = 2700.0;

		rho[i] = P;
		mu[i] = Vs*Vs*rho[i];
		lamda[i] = Vp*Vp*rho[i] - 2.0*mu[i];
	}
}

void Change_vel_sem(short num, long elem, long **Elem_array, float *rho, float *lamda,float *mu, double **node_array)
{
	int i, j, k, i1, i2, i3, i4, num_dp, count, fs, Cin[4];
	char ch;
	double Vp, Vs, rr, dis, flgfx1, flgfz1, flgfx2, flgfz2, flgfx3, flgfz3, flgfx4, flgfz4, flgf1, flgf2, tmp;
	double xx, zz, dd, vvp, vvs, layer, num_layer[16], **vv;

	FILE *fp;

	fp = fopen("../Input/vel_model.txt", "r");

	if(fp == NULL)
	{
		printf("cannot open this file -- vel_model.txt\n"); exit(1);
	}

	count = 0;
	while ((ch = fgetc(fp)) != EOF)
	{
		if(ch == '\n')
		{
			count++;
		}
	}
	rewind(fp);

	vv = Arraydouble2(count, 4);

	for (i = 0; i < count; i++)
	{
		fs = fscanf(fp, "%lf%lf%lf%lf", &xx, &zz, &vvp, &vvs);

		vv[i][0] = xx;  vv[i][1] = zz;
		vv[i][2] = vvp; vv[i][3] = vvs;

		if (i < 16)
		{
			num_layer[i] = zz;
		}
	}
	
	for(i = 0; i < elem; i++)
	{
		Cin[0] = 0;
		Cin[1] = num * (num - 1);
		Cin[2] = num * num - 1;
		Cin[3] = num - 1;

		i1 = Elem_array[i][Cin[0]] - 1, i2 = Elem_array[i][Cin[1]] - 1;
		i3 = Elem_array[i][Cin[2]] - 1, i4 = Elem_array[i][Cin[3]] - 1;

		flgfx1 = node_array[i1][0], flgfz1 = node_array[i1][1];
		flgfx2 = node_array[i2][0], flgfz2 = node_array[i2][1];
		flgfx3 = node_array[i3][0], flgfz3 = node_array[i3][1];
		flgfx4 = node_array[i4][0], flgfz4 = node_array[i4][1];

		flgf1 = (flgfx1 + flgfx2 + flgfx3 + flgfx4) / 4;
		flgf2 = (flgfz1 + flgfz2 + flgfz3 + flgfz4) / 4;

		// tmp = flgf2;

		for (j = 0; j < 15; j++)
		{
			if (flgf2 > num_layer[15])
			{
				flgf2 = num_layer[15];
				break;
			}
			else if(flgf2 < num_layer[0])
			{
				flgf2 = num_layer[0];
				break;
			}
			else if (flgf2 <= num_layer[j + 1] && flgf2 >= num_layer[j])
			{
				flgf2 = num_layer[j];
				break;
			}
		}

		rr = 4000;

		// for (j = 0; j < count; j++)
		// {
		// 	xx = vv[j][0];  zz = vv[j][1];
		// 	vvp = vv[j][2]; vvs =vv[j][3];

		// 	dis = sqrt((xx - flgf1) * (xx - flgf1) + (zz - flgf2) * (zz - flgf2));

		// 	if (dis <= rr)
		// 	{
		// 		Vp = vvp; Vs = vvs;

		// 		rr = dis;
		// 	}
		// }

		Vp = vv[j][2]; Vs = vv[j][3];

		rho[i] = 2600.0;

		mu[i] = Vs * Vs * rho[i];

		lamda[i] = Vp * Vp * rho[i] - 2.0 * mu[i];
	}

	Freedouble2(count, vv);
	fclose(fp);
}

void Change_vel_sem_s1(short num, long elem, long **Elem_array, float *rho, float *lamda,float *mu, double **node_array)
{
	int i, j, k, i1, i2, i3, i4, num_dp, count, fs, Cin[4];
	char ch;
	double Vp, Vs, P, rr, dis, flgfx1, flgfz1, flgfx2, flgfz2, flgfx3, flgfz3, flgfx4, flgfz4, flgf1, flgf2, tmp;
	double xx, yy, dd, vvp, vvs, layer, num_layer[16], **vv;
	double x1, x2, y1, y2, fk, FF1;

	x1 = 0.0;      x2 = 2000.0;
	y1 = 1200.0;   y2 = 1500.0;
	
	for(i = 0; i < elem; i++)
	{
		Cin[0] = 0;
		Cin[1] = num * (num - 1);
		Cin[2] = num * num - 1;
		Cin[3] = num - 1;

		i1 = Elem_array[i][Cin[0]] - 1, i2 = Elem_array[i][Cin[1]] - 1;
		i3 = Elem_array[i][Cin[2]] - 1, i4 = Elem_array[i][Cin[3]] - 1;

		flgfx1 = node_array[i1][0], flgfz1 = node_array[i1][1];
		flgfx2 = node_array[i2][0], flgfz2 = node_array[i2][1];
		flgfx3 = node_array[i3][0], flgfz3 = node_array[i3][1];
		flgfx4 = node_array[i4][0], flgfz4 = node_array[i4][1];

		xx = (flgfx1 + flgfx2 + flgfx3 + flgfx4) / 4.0;
		yy = (flgfz1 + flgfz2 + flgfz3 + flgfz4) / 4.0;

		fk = (y2 - y1) / (x2 - x1);
//		fk = 0.0; y1 = 1200.0;

		FF1 = fk*(xx - x1) + y1;

		if (yy <= FF1)
		{
			Vp = 3000.0;
			Vs = 2000.0;
			P = 2000.0;
		}
		else
		{
			Vp = 4000.0;
			Vs = 2500.0;
			P = 2500.0;
		}

		Vp = 5000.0;
		Vs = 2800.0;
		P = 2600.0;


		rho[i] = P;
		mu[i] = Vs * Vs * rho[i];
		lamda[i] = Vp * Vp * rho[i] - 2.0 * mu[i];
	}
}

void Change_vel_fdm(int L, int Nx, int Nz, float **lamda, float **mu, float **rho, double **FDM_array_x, double **FDM_array_z)
{
	int i, j, k, count, fs;

	char ch;

	double Vp, Vs, rr, dis, flgf1, flgf2;

	double xx, zz, dd, vvp, vvs, layer, num_layer[16], **vv;

	FILE *fp;

	fp = fopen("../Input/vel_model.txt", "r");

	if(fp == NULL)
	{
		printf("cannot open this file -- vel_model.txt\n"); exit(1);
	}

	count = 0;
	while ((ch = fgetc(fp)) != EOF)
	{
		if(ch == '\n')
		{
			count++;
		}
	}
	rewind(fp);

	vv = Arraydouble2(count, 4);

	for (i = 0; i < count; i++)
	{
		fs = fscanf(fp, "%lf%lf%lf%lf", &xx, &zz, &vvp, &vvs);

		vv[i][0] = xx;  vv[i][1] = zz;
		vv[i][2] = vvp; vv[i][3] = vvs;

		if (i < 16)
		{
			num_layer[i] = zz;
		}
	}

	for (i = 0; i <= Nx + 2 * L; i++)
	{
		for (j = 0; j <= Nz + L; j++)
		{
			flgf1 = FDM_array_x[i][j];
			flgf2 = FDM_array_z[i][j];

			for (k = 0; k < 15; k++)
			{
				if (flgf2 > num_layer[15])
				{
					flgf2 = num_layer[15];
					break;
				}
				else if(flgf2 < num_layer[0])
				{
					flgf2 = num_layer[0];
					break;
				}
				else if (flgf2 <= num_layer[k + 1] && flgf2 >= num_layer[k])
				{
					flgf2 = num_layer[k];
					break;
				}
			}

			rr = 4000;

			// for (k = 0; k < count; k++)
			// {
			// 	xx = vv[k][0];  zz = vv[k][1];
			// 	vvp = vv[k][2]; vvs =vv[k][3];

			// 	dis = sqrt((xx - flgf1) * (xx - flgf1) + (zz - flgf2) * (zz - flgf2));

			// 	if (dis <= rr)
			// 	{
			// 		Vp = vvp; Vs = vvs;

			// 		rr = dis;
			// 	}
			// }

			Vp = vv[k][2]; Vs = vv[k][3];

			rho[i][j] = 2600.0;

			mu[i][j] = rho[i][j] * Vs * Vs;

			lamda[i][j] = rho[i][j] * Vp * Vp - 2.0 * mu[i][j];
		}
	}

	Freedouble2(count, vv);
	fclose(fp);
}

////////////////////////////////////////////fdm define lamda, mu and rho////////////////////////////////////////////////
void getvel(int L, int Nx, int Nz, float **lamda, float **mu, float **rho, double **FDM_array_x, double **FDM_array_z)
{
	int i, j,index1, index2, startx, endx, startz, endz;

	double Vp, Vs;

	startx = L;  endx = L + Nx;
	startz = 0;  endz = Nz;

	for (i = 0; i <= Nx + 2 * L; i++)
	{
		for (j = 0; j <= Nz + L; j++)
		{
			Vp=6000.0;
			Vs=3500.0;

//			Vp=3200.0;
//			Vs=1800.0;

			rho[i][j]=2600.0;

			mu[i][j] = rho[i][j] * Vs * Vs;

			lamda[i][j] = rho[i][j] * Vp * Vp - 2.0 * mu[i][j];
		}
	}
}

////////////////////////////////differential coefficients/////微分系数///////////
void Coeffs(int N1,double **Coeff1,double **Coeff2)
{
	double *Arraydouble1(long Nx);
	void Coeff_array(int dir,int N1,double *Coeff);
	double *Coeff;
	int i,j;

	Coeff=Arraydouble1(N1);

	for(i=0;i<N1;i++)
	{
		Coeff_array(1,i+1,Coeff);
		for(j=0;j<N1;j++)	
			Coeff1[i][j]=Coeff[j];

		Coeff_array(2,i+1,Coeff);
		for(j=0;j<N1;j++)
			Coeff2[i][j]=Coeff[j];
	}
}

void Coeff_array(int dir,int N1,double *Coeff)
{
	if(dir==1)
	{
		if(N1==1)
			Coeff[0]=0.5;
		else if(N1==2)
		{
			Coeff[0]=6.666666666666666300e-001;
			Coeff[1]=-8.333333333333332900e-002;
		}
		else if(N1==3)
		{
			Coeff[0]=7.500000000000000000e-001;
			Coeff[1]=-1.499999999999999900e-001;
			Coeff[2]=1.666666666666666600e-002;	
		}
		else if(N1==4)
		{
			Coeff[0]=7.999999999999998200e-001;
			Coeff[1]=-1.999999999999999800e-001;
			Coeff[2]=3.809523809523809900e-002;
			Coeff[3]=-3.571428571428571800e-003;
		}
		else if(N1==5)
		{
			Coeff[0]=8.789874e-001;
			Coeff[1]=-2.951612e-001;
			Coeff[2]=9.706293e-002;
			Coeff[3]=-2.398844e-002;
			Coeff[4]=3.259713e-003;	
		}
		else if(N1==6)
		{
			Coeff[0]=8.571428571428574300e-001;
			Coeff[1]=-2.678571428571429600e-001;
			Coeff[2]=7.936507936507945800e-002;
			Coeff[3]=-1.785714285714290100e-002;
			Coeff[4]=2.597402597402607400e-003;
			Coeff[5]=-1.803751803751813000e-004;
		}
		else if(N1==7)
		{
			Coeff[0]=8.750000000000007800e-001;
			Coeff[1]=-2.916666666666679100e-001;
			Coeff[2]=9.722222222222309800e-002;
			Coeff[3]=-2.651515151515190800e-002;
			Coeff[4]=5.303030303030411800e-003;
			Coeff[5]=-6.798756798756967900e-004;
			Coeff[6]=4.162504162504277000e-005;
		}
		else if(N1==8)
		{
			Coeff[0]=8.888888888888929500e-001;
			Coeff[1]=-3.111111111111163300e-001;
			Coeff[2]=1.131313131313170000e-001;
			Coeff[3]=-3.535353535353724000e-002;
			Coeff[4]=8.702408702409316500e-003;
			Coeff[5]=-1.554001554001683900e-003;
			Coeff[6]=1.776001776001941300e-004;
			Coeff[7]=-9.712509712510691700e-006;
		}
	}
	else if(dir==2)
	{
		if(N1==1)
			Coeff[0]=1.0;
		else if(N1==2)
		{
			Coeff[0]=1.333333333333333300e+000;
			Coeff[1]=-8.333333333333332900e-002;
		}
		else if(N1==3)
		{
			Coeff[0]=1.500000000000000000e+000;
			Coeff[1]=-1.499999999999999900e-001;
			Coeff[2]=1.111111111111111200e-002;
		}
		else if(N1==4)
		{
			Coeff[0]=1.599999999999999600e+000;
			Coeff[1]=-1.999999999999999800e-001;
			Coeff[2]=2.539682539682539700e-002;
			Coeff[3]=-1.785714285714285900e-003;
		}
		else if(N1==5)
		{
			Coeff[0]=1.666666666666666500e+000;
			Coeff[1]=-2.380952380952380500e-001;
			Coeff[2]=3.968253968253967300e-002;
			Coeff[3]=-4.960317460317458300e-003;
			Coeff[4]=3.174603174603173000e-004;
		}
		else if(N1==6)
		{
			Coeff[0]=1.714285714285714000e+000;
			Coeff[1]=-2.678571428571426800e-001;
			Coeff[2]=5.291005291005284500e-002;
			Coeff[3]=-8.928571428571409000e-003;
			Coeff[4]=1.038961038961035500e-003;
			Coeff[5]=-6.012506012505984900e-005;
		}
		else if(N1==7)
		{
			Coeff[0]=1.749999999999998400e+000;
			Coeff[1]=-2.916666666666656300e-001;
			Coeff[2]=6.481481481481432500e-002;
			Coeff[3]=-1.325757575757559300e-002;
			Coeff[4]=2.121212121212084600e-003;
			Coeff[5]=-2.266252266252217300e-004;
			Coeff[6]=1.189286903572587900e-005;
		}
		else if(N1==8)
		{
			Coeff[0]=1.777777777777772400e+000;
			Coeff[1]=-3.111111111111075100e-001;
			Coeff[2]=7.542087542087357000e-002;
			Coeff[3]=-1.767676767676699300e-002;
			Coeff[4]=3.480963480963301000e-003;
			Coeff[5]=-5.180005180004861600e-004;
			Coeff[6]=5.074290788576158800e-005;
			Coeff[7]=-2.428127428127255300e-006;
		}
	}
	
}
double diff(int Nx,int Nz,int order,int Ni,int Nj,int d,int style,double Dx,double Dz,float **V,
			double **coeff1,double **coeff2)
{
	double diff1(int Ni,int Nj,int style,int t1,double Dx,double Dz,float **V,double **coeff1);
	double diff2(int Ni,int Nj,int style,int t2,double Dx,double Dz,float **V,double **coeff2);
	double t=0.0,tmp1,tmp2;
	int i,t1,t2; 
	
	t1=order;
	t2=order;

	if(t1>Ni)
		t1=Ni;
	else if(t1>Nx-Ni)
		t1=Nx-Ni;
	
	if(t2>Nj)
		t2=Nj;
	else if(t2>Nz-Nj)
		t2=Nz-Nj;

    if(d==1&&style==1)							//垂直对x
		t=diff1(Ni,Nj,1,t1,Dx,Dz,V,coeff1);
	else if(d==1&&style==2)						//垂直对z
		t=diff1(Ni,Nj,2,t2,Dx,Dz,V,coeff1);
	else if(d==2&&style==1)						//水平对x
		t=diff2(Ni,Nj,1,t1,Dx,Dz,V,coeff2);
	else if(d==2&&style==2)						//水平对z
		t=diff2(Ni,Nj,2,t2,Dx,Dz,V,coeff2);

	else if(d==0&&style==0)
	{
		for(i=1;i<=t1;i++)
		{
			tmp1=diff1(Ni+i,Nj,2,t2,Dx,Dz,V,coeff1);	//垂直对zx
			tmp2=diff1(Ni-i,Nj,2,t2,Dx,Dz,V,coeff1);
			t+=coeff1[t1-1][i-1]*(tmp1-tmp2)/Dx;
		}
	}
   return t;
}

double diff1(int Ni,int Nj,int style,int t1,double Dx,double Dz,float **V,double **coeff1)
{
	int i;
	double t;
	t=0.0;

	if(style==1)
		for(i=1;i<=t1;i++)
			t+=coeff1[t1-1][i-1]*(V[Ni+i][Nj]-V[Ni-i][Nj])/Dx;

	else if(style==2)
		for(i=1;i<=t1;i++)
			t+=coeff1[t1-1][i-1]*(V[Ni][Nj+i]-V[Ni][Nj-i])/Dz;

	return t;
}


double diff2(int Ni,int Nj,int style,int t2,double Dx,double Dz,float **V,double **coeff2)	
{
	int i;
	double t=0.0;

	if(style==1)
	{
		for(i=1;i<=t2;i++)
			{
			//if(((Ni>=0&&Ni<=5)||(Ni>=395&&Ni<=400))&&(Nj==0||Nj==1||Nj==2))
			//printf("V-1=%e,V0=%e,V1%e\n",V[Ni+i][Nj],V[Ni-i][Nj],V[Ni][Nj]);
			t+=coeff2[t2-1][i-1]*(V[Ni+i][Nj]+V[Ni-i][Nj]-2.0*V[Ni][Nj])/(Dx*Dx);
			//if(((Ni>=0&&Ni<=5)||(Ni>=395&&Ni<=400))&&(Nj==0||Nj==1||Nj==2))
			//printf("i=%d,j=%d,t=%e,coeff2=%e\n",Ni,Nj,t,coeff2[t2-1][i-1]);
			}
//		printf("test-diff,Ni=%d,Nj=%d\n",Ni,Nj);
//		getchar();
	}
	else if(style==2)
		for(i=1;i<=t2;i++)
			t+=coeff2[t2-1][i-1]*(V[Ni][Nj+i]+V[Ni][Nj-i]-2.0*V[Ni][Nj])/(Dz*Dz);
	return t;
}


////////////////////////////////////////////////////fdm/////////////////////////////////////////////////////////////////////////////////////////////////

double diff_f(int Nx,int Nz,int order,int Ni,int Nj,int d,int style,double Dx,double Dz,float **V,
			double **coeff1,double **coeff2,float **lamda1,float **mu1,float **rho1)
{
	double diff1_f(int Ni,int Nj,int style,int t1,double Dx,double Dz,float **V,double **coeff1,float **lamda1,float **mu1,float **rho1);
	double diff2_f(int Ni,int Nj,int style,int t2,double Dx,double Dz,float **V,double **coeff2,float **lamda1,float **mu1,float **rho1);
	double t=0.0,tmp1,tmp2;
	int i,t1,t2; 
	
	t1=order;
	t2=order;

	if(t1>Ni)
		t1=Ni;
	else if(t1>Nx-Ni)
		t1=Nx-Ni;
	
	if(t2>Nj)
		t2=Nj;
	else if(t2>Nz-Nj)
		t2=Nz-Nj;

    if(d==1&&style==1)							//垂直对x
		t=diff1_f(Ni,Nj,1,t1,Dx,Dz,V,coeff1,lamda1,mu1,rho1);
	else if(d==1&&style==2)						//垂直对z
		t=diff1_f(Ni,Nj,2,t2,Dx,Dz,V,coeff1,lamda1,mu1,rho1);
	else if(d==2&&style==1)						//水平对xx
		t=diff2_f(Ni,Nj,1,t1,Dx,Dz,V,coeff2,lamda1,mu1,rho1);
	else if(d==2&&style==2)						//水平对zz
		t=diff2_f(Ni,Nj,2,t2,Dx,Dz,V,coeff2,lamda1,mu1,rho1);

	else if(d==0&&style==0)
	{
		for(i=1;i<=t1;i++)
		{
			tmp1=diff1_f(Ni+i,Nj,2,t2,Dx,Dz,V,coeff1,lamda1,mu1,rho1);	//垂直对zx
			tmp2=diff1_f(Ni-i,Nj,2,t2,Dx,Dz,V,coeff1,lamda1,mu1,rho1);
			t+=coeff1[t1-1][i-1]*(tmp1-tmp2)/Dx;
		}
	}
   return t;
}
double diff1_f(int Ni,int Nj,int style,int t1,double Dx,double Dz,float **V,double **coeff1,float **lamda1,float **mu1,float **rho1)
{
	int i;
	double t;
	t=0.0;

	if(style==1)
		for(i=1;i<=t1;i++)
			t+=coeff1[t1-1][i-1]*(V[Ni+i][Nj]-V[Ni-i][Nj])/Dx;

	else if(style==2)
		for(i=1;i<=t1;i++)
			t+=coeff1[t1-1][i-1]*(V[Ni][Nj+i]-V[Ni][Nj-i])/Dz;

	return t;
}


double diff2_f(int Ni,int Nj,int style,int t2,double Dx,double Dz,float **V,double **coeff2,float **lamda1,float **mu1,float **rho1)	
{
	int i;
	double t=0.0;

	if(style==1)
	{
		for(i=1;i<=t2;i++)
			{
			t+=coeff2[t2-1][i-1]*(V[Ni+i][Nj]+V[Ni-i][Nj]-2.0*V[Ni][Nj])/(Dx*Dx);
			}
	}
	else if(style==2)
		for(i=1;i<=t2;i++)
			t+=coeff2[t2-1][i-1]*(V[Ni][Nj+i]+V[Ni][Nj-i]-2.0*V[Ni][Nj])/(Dz*Dz);
	return t;
}

/////////////////////////////////allocate arrays//////////////////////////////
//////////////////////////////////分配数组//////////////////////////////////
double *Arraydouble1(long Nx)
{
	long i;
	double *Array;
	Array=(double *)malloc(Nx*sizeof(double));
	for(i=0;i<Nx;i++)
		Array[i]=0.0;
	return Array;
}

double **Arraydouble2(long Nx,long Nz)
{
	long i,j;
	double **Array;
	Array=(double **)malloc(Nx*sizeof(double *));
	for(i=0;i<Nx;i++)
		*(Array+i)=(double *)malloc(Nz*sizeof(double));
	for(i=0;i<Nx;i++)
		for(j=0;j<Nz;j++)
			Array[i][j]=0.0;
	return Array;
}

double ***Arraydouble3(long N1, long N2, long N3)
{
	long i, j, k;
	double ***Array;

	Array = (double ***)malloc(N1 * sizeof(double **));

	for(i = 0; i < N1; i++)
	{
		Array[i] = (double **)malloc(N2 * sizeof(double *));

		for(j = 0; j < N2; j++)
		{
			Array[i][j] = (double *)malloc(N3 * sizeof(double));			
		}		
	}

	for(i = 0; i < N1; i++)
	{
		for(j = 0; j < N2; j++)
		{
			for(k = 0; k < N3; k++)
			{
				Array[i][j][k] = 0.0;	
			}
		}
	}

	return Array;
}

double ****Arraydouble4(long N1, long N2, long N3, long N4)
{
	long i, j, k, l;
	double ****Array;

	Array = (double ****)malloc(N1 * sizeof(double ***));

	for(i = 0; i < N1; i++)
	{
		Array[i] = (double ***)malloc(N2 * sizeof(double **));

		for(j = 0; j < N2; j++)
		{
			Array[i][j] = (double **)malloc(N3 * sizeof(double *));

			for (k = 0; k < N3; k++)
			{
				Array[i][j][k] = (double *)malloc(N4 * sizeof(double));
			}
		}
	}

	for(i = 0; i < N1; i++)
	{
		for(j = 0; j < N2; j++)
		{
			for(k = 0; k < N3; k++)
			{
				for(l = 0; l < N4; l++)
				{
					Array[i][j][k][l] = 0.0;
				}
			}
		}
	}

	return Array;
}


int Freedouble2(int Nx,double **Array)
{
	int i;

	for(i = 0; i < Nx; i++)
		free(Array[i]);

	free(Array);
	return 0;
}

void Freedouble3(long N1, long N2, double ***Array)
{
	long i, j;

	for(i = 0; i < N1; i++)
	{
		for(j=0;j<N2;j++)
		{
			free(Array[i][j]);
		}
		free(Array[i]);
	}

	free(Array);
}

void Freedouble4(long N1, long N2, long N3, double ****Array)
{
	long i, j, k;

	for(i = 0; i < N1; i++)
	{
		for(j = 0; j < N2; j++)
		{
			for(k = 0; k < N3; k++)
			{
				free(Array[i][j][k]);
			}
			free(Array[i][j]);
		}
		free(Array[i]);
	}

	free(Array);
}

float **Arrayfloat2(long Nx,long Nz)
{
	long i,j;
	float **Array;
	Array=(float **)malloc(Nx*sizeof(float *));
	for(i=0;i<Nx;i++)
		*(Array+i)=(float *)malloc(Nz*sizeof(float));
	for(i=0;i<Nx;i++)
		for(j=0;j<Nz;j++)
			Array[i][j]=0.0;
	return Array;
}

float ****Arrayfloat4(long N1,long N2,long N3,long N4)
{
	long i,j,k,l;
	float ****Array;
	Array=(float ****)malloc(N1*sizeof(float ***));
	for(i=0;i<N1;i++)
	{
		Array[i]=(float ***)malloc(N2*sizeof(float **));
		for(j=0;j<N2;j++)
		{
			Array[i][j]=(float **)malloc(N3*sizeof(float *));
			for(k=0;k<N3;k++)
			{
				Array[i][j][k]=(float *)malloc(N4*sizeof(float));
			}
		}
	}

	for(i = 0; i < N1; i++)
	{
		for(j = 0; j < N2; j++)
		{
			for(k = 0; k < N3; k++)
			{
				for(l = 0; l < N4; l++)
				{
					Array[i][j][k][l] = 0.0;						
				}
			}
		}
	}

	return Array;
}


void Freefloat4(long N1,long N2,long N3,float ****Array)
{
	long i,j,k;
	for(i=0;i<N1;i++)
	{
		for(j=0;j<N2;j++)
		{
			for(k=0;k<N3;k++)
			{
				free(Array[i][j][k]);
			}
			free(Array[i][j]);
		}
		free(Array[i]);
	}
	free(Array);
}

long **Arraylong2(long Nx,long Nz)
{
	long i,j;
	long **Pt;
	Pt=(long **)malloc(Nx*sizeof(long *));
	for(i=0;i<Nx;i++)
		*(Pt+i)=(long *)malloc(Nz*sizeof(long));
	for(i=0;i<Nx;i++)
		for(j=0;j<Nz;j++)
			Pt[i][j]=0;
	return Pt;
}

int **Arrayint2(long Nx,long Nz)
{
	long i,j;
	int **Array;
	Array=(int **)malloc(Nx*sizeof(int *));
	for(i=0;i<Nx;i++)
		*(Array+i)=(int *)malloc(Nz*sizeof(int));
	for(i=0;i<Nx;i++)
		for(j=0;j<Nz;j++)
			Array[i][j]=0;
	return Array;
}



float *Arrayfloat1(long Nx)
{
	long i;
	float *Array;
	Array=(float *)malloc(Nx*sizeof(float));
	for(i=0;i<Nx;i++)
		Array[i]=0.0;
	return Array;
}

long *Arraylong1(long Nx)
{
	long i;
	long *Array;
	Array=(long *)malloc(Nx*sizeof(long));
	for(i=0;i<Nx;i++)
		Array[i]=0;
	return Array;
}

int *Arrayint1(long Nx)
{
	long i;
	int *Array;
	Array=(int *)malloc(Nx*sizeof(int));
	for(i=0;i<Nx;i++)
		Array[i]=0;
	return Array;
}


short *Arrayshort1(long Nx)
{
	long i;
	short *Array;
	Array=(short *)malloc(Nx*sizeof(short));
	for(i=0;i<Nx;i++)
		Array[i]=0;
	return Array;
}

int Freeint2(int Nx, int **Array)
{
	int i;
	for(i=0;i<Nx;i++)
		free(Array[i]);
	free(Array);
	return 0;
}

int Freefloat2(int Nx, float **Array)
{
	int i;
	for(i=0;i<Nx;i++)
		free(Array[i]);
	free(Array);
	return 0;
}

int Freelong2(int Nx,long **Array)
{
	int i;
	for(i=0;i<Nx;i++)
		free(Array[i]);
	free(Array);
	return 0;
}

char *Arraychar1(long Nx)
{
	long i;
	char *Array;

	Array=(char *)malloc(Nx*sizeof(char));
	Array[0]='\0';

	return Array;
}
///////////////////////////GLL weight and GLL position////////////////////////////////////
void Weight(int Nx,double *W)
{
	if(Nx==1)
	{
		W[0]=1.0;
		W[1]=1.0;
	}
	else if(Nx==2)
	{
		W[0]=0.33333333333333;
		W[1]=1.333333333333333;
		W[2]=0.33333333333333;
	}
	else if(Nx==3)
	{
		W[0]=0.16666666666667;
		W[1]=0.83333333333333;
		W[2]=0.83333333333333;
		W[3]=0.16666666666667;
	}
	else if(Nx==4)
	{
		W[0]=0.10000000000000;
		W[1]=0.54444444444444;
		W[2]=0.71111111111112;
		W[3]=0.54444444444444;
		W[4]=0.10000000000000;
	}
	else if(Nx==5)
	{
		W[0]=0.06666666666667;
		W[1]=0.3784749562;
		W[2]=0.5548583766;
		W[3]=0.5548583766;
		W[4]=0.3784749562;
		W[5]=0.06666666666667;
	}
	else if(Nx==6)
	{
		W[0]=0.04761904761905;
		W[1]=0.2768260470;
		W[2]=0.4317453812;
		W[3]=0.4876190476;
		W[4]=0.43174538121;
		W[5]=0.2768260470;
		W[6]=0.04761904761905;
	}
	else if(Nx==7)
	{
		W[0]=0.03571428571429;
		W[1]=0.2107042271;
		W[2]=0.3411226922;
		W[3]=0.4124587946;
		W[4]=0.4124587946;
		W[5]=0.3411226922;
		W[6]=0.2107042271;
		W[7]=0.03571428571429;
	}
	else if(Nx==8)
	{
		W[0]=0.02777777777778;
		W[1]=0.1654953612;
		W[2]=0.2745387130;
		W[3]=0.3464285109;
		W[4]=0.3715192744;
		W[5]=0.3464285109;
		W[6]=0.2745387130;
		W[7]=0.1654953612;
		W[8]=0.02777777777778;
	}
	else if(Nx==9)
	{
		W[0]=0.02222222222222;
		W[1]=0.1333059906;
		W[2]=0.2248893424;
		W[3]=0.2920426833;
		W[4]=0.3275397613;	
		W[5]=0.3275397613;
		W[6]=0.2920426833;
		W[7]=0.2248893424;
		W[8]=0.1333059906;
		W[9]=0.02222222222222;
	}
}
void Position(int Nx,double *P)
{
	if(Nx==1)
	{
		P[0]=-1.0;
		P[1]=1.0;
	}
	else if(Nx==2)
	{
		P[0]=-1.0;
		P[1]=0.0;
		P[2]=1.0;
	}
	else if(Nx==3)
	{
		P[0]=-1.0;
		P[1]=-0.44721359549996;
		P[2]=0.44721359549996;
		P[3]=1.0;
	}
	else if(Nx==4)
	{
		P[0]=-1.0;
		P[1]=-0.65465367070798;
		P[2]=0.0;
		P[3]=0.65465367070798;
		P[4]=1.0;
	}
	else if(Nx==5)
	{
		P[0]=-1.0;
		P[1]=-0.76505532392946;
		P[2]=-0.28523151648065;
		P[3]=0.28523151648065;
		P[4]=0.76505532392946;
		P[5]=1.0;
	}
	else if(Nx==6)
	{
		P[0]=-1.0;
		P[1]=-0.83022389627857;
		P[2]=-0.46884879347071;
		P[3]=0.0;
		P[4]=0.46884879347071;
		P[5]=0.83022389627857;
		P[6]=1.0;
	}
	else if(Nx==7)
	{
		P[0]=-1.0;
		P[1]=-0.8717401487;
		P[2]=-0.5917001814;
		P[3]=-0.2092992172;
		P[4]=0.2092992172;
		P[5]=0.5917001814;
		P[6]=0.8717401487;
		P[7]=1.0;
	}
	else if(Nx==8)
	{
		P[0]=-1.0;
		P[1]=-0.8997579958;
		P[2]=-0.6771862797;
		P[3]=-0.3631174632;
		P[4]=0.0;
		P[5]=0.3631174632;
		P[6]=0.6771862797;
		P[7]=0.8997579958;
		P[8]=1.0;
	}
	else if(Nx==9)
	{
		P[0]=-1.0;
		P[1]=-0.919533908166459;
		P[2]=-0.738773865105505;
		P[3]=-0.477924949810444;
		P[4]=-0.165278957666387;	
		P[5]=0.165278957666387;
		P[6]=0.477924949810444;
		P[7]=0.738773865105505;
		P[8]=0.919533908166459;
		P[9]=1.0;
	}
}
////////////////////////////////////difference of lagrange function////拉格朗日的差分/////////////////
double diff_fun(int Nx,int i1,int i2)
{
	double Val,*V,*P;
	void Val_legendre(int Nx,double *V);
	V=Arraydouble1(Nx+1);
	P=Arraydouble1(Nx+1);
	Val_legendre(Nx,V);
	Position(Nx,P);

	if(i1==0&&i2==0)
		Val=-1.0/4.0*Nx*(Nx+1.0);
	else if(i1==Nx&&i2==Nx)
		Val=1.0/4.0*Nx*(Nx+1.0);
	else if(i1>=1&&i1<=Nx-1&&i1==i2)
		Val=0.0;
	else
		Val=V[i2]/(V[i1]*(P[i2]-P[i1]));

	free(V);
	free(P);
	return Val;
}
/////////////////////////////////////////Value of Legendre/////////////////////////////////////
void Val_legendre(int Nx,double *V)
{
	if(Nx==1)
	{
		V[0]=1.0;
		V[1]=1.0;	
	}
	else if(Nx==2)
	{
		V[0]=1.0;
		V[1]=-0.50;
		V[2]=1.0;
	}
	else if(Nx==3)
	{
		V[0]=-1.000000000;
		V[1]=0.4472135955;
		V[2]=-0.4472135955;
		V[3]=1.000000000;
	}
	else if(Nx==4)
	{
		V[0]=1.000000000;
		V[1]=-0.4285714285;
		V[2]=0.3750000000;
		V[3]=-0.4285714285;
		V[4]=1.000000000;
	}
	else if(Nx==5)
	{
		V[0]=-1.000000000;
		V[1]=0.4196969342;
		V[2]=-0.3466277252;
		V[3]=0.3466277252;
		V[4]=-0.4196969342;
		V[5]=1.0;
	}
	else if(Nx==6)
	{
		V[0]=1.000000000;
		V[1]=-0.4147504607;
		V[2]=0.3321058323;
		V[3]=-0.3125000000;
		V[4]=0.3321058323;
		V[5]=-0.4147504607;
		V[6]=1.0;
	}
	else if(Nx==7)
	{
		V[0]=-1.0;
		V[1]=0.4117033114;
		V[2]=-0.3235680853;
		V[3]=0.2942596406;
		V[4]=-0.2942596406;
		V[5]=0.3235680853;
		V[6]=-0.4117033114;
		V[7]=1.0;
	}
	else if(Nx==8)
	{
		V[0]=1.000000000;
		V[1]=-0.4096904467;
		V[2]=0.3180877564;
		V[3]=-.2831665412;
		V[4]=0.2734375000;
		V[5]=-0.2831665412;
		V[6]=0.3180877564;
		V[7]=-0.4096904467;
		V[8]=1.000000000;
	}
	else if(Nx==9)
	{
		V[0]=-1.0;
		V[1]=0.4082901567;
		V[2]=-0.3143469898;
		V[3]=0.2758484691;
		V[4]=-0.2604724104;
		V[5]=0.2604724104;
		V[6]=-0.2758484691;
		V[7]=0.3143469898;
		V[8]=-0.4082901567;
		V[9]=1.0;
	}

}
	
void Common_array(short num,double **ND,double **NT)
{
	short i1,i2,NP;
	NP=num-1;
//	printf("NP=%d\n",NP);
	for(i1=0;i1<num;i1++)
		for(i2=0;i2<num;i2++)
		{
			if(i1!=i2)
				ND[i1][i2]=0.0;
			else
				ND[i1][i2]=1.0;
		}

	for(i1=0;i1<num;i1++)
		for(i2=0;i2<num;i2++)
		{
			NT[i1][i2]=diff_fun(NP,i1,i2);
			//printf("%d,%d,%f\n",i1,i2,NT[i1][i2]);
		}
	//	getchar();
}


double Lagrange_val(short num, short ln, double *P, double lx)
{
	short i;
	double tmp1, tmp2;

	tmp1 = 1.0; tmp2 = 1.0;

	for(i = 0; i < num; i++)
	{
		if(i == ln) continue;

		tmp1 *= (P[ln] - P[i]);
		tmp2 *= (lx - P[i]);
	}

	tmp1 = tmp2 / tmp1;

	return tmp1;
}

double Lagrange_val_1(short num,short ln,double *P,double lx)
{
	short i;
	double tmp1,tmp2;

	tmp1=1.0;
	tmp2=1.0;
	for(i=1;i<num;i++)
	{
		if(i==ln)continue;
		tmp1*=(P[ln]-P[i]);
		tmp2*=(lx-P[i]);
	}
	tmp1=tmp2/tmp1;
	return tmp1;
}

double Lagrange_diff(short num,short ln,double *P,double lx)
{
	double tmp,tmp1,tmp2,dx;
	if(fabs(lx+1)<=1.0e-10)
	{
		dx=1.0e-10;
		tmp1=Lagrange_val(num,ln,P,lx);
		tmp2=Lagrange_val(num,ln,P,dx+lx);
		tmp=(tmp2-tmp1)/(dx);
//printf("11 %f,%f,%f\n",lx,tmp1,tmp2);
	}
	else if(fabs(lx-1)<=1.0e-10)
	{
		dx=1.0e-10;
		tmp1=Lagrange_val(num,ln,P,lx-dx);
		tmp2=Lagrange_val(num,ln,P,lx);
		tmp=(tmp2-tmp1)/(dx);
//printf("22 %f,%f,%f\n",lx,tmp1,tmp2);
	}
	else
	{	
		dx=1.0e-10;
		tmp1=Lagrange_val(num,ln,P,-dx+lx);
		tmp2=Lagrange_val(num,ln,P,dx+lx);	
		tmp=(tmp2-tmp1)/(2.0*dx);
//printf("33 %e,%e,%e,%e\n",lx,tmp1,tmp2,tmp);
	}
	return tmp;
}

int count_file(char *FL)
{
	int i, sf, count;

	char ch, FF[100];

	double tmp1, tmp2, tmp3, tmp4, tmp5;

	FILE *fp;

	fp = fopen(FL,"r");

	if(fp == NULL)
	{
		printf("cannot open this file -- %s\n", FL); exit(1);
	}

	count = 0;
	while ((ch = fgetc(fp)) != EOF)
	{
		if(ch == '\n')
		{
			count++;
		}
	}

	rewind(fp);
	fclose(fp);

	return count;
}

void rdfile(char *FL, int count, int n, double **ary)
{
	int i, sf;

	char ch, FF[100], *ff;

	double tmp1, tmp2, tmp3, tmp4, tmp5;

	FILE *fp;

	fp = fopen(FL,"r");

	if (fp == NULL)
	{
		printf("can't open this file --- %s !\n", FL), exit(1);		
	}

	count = 0;
	while ((ch = fgetc(fp)) != EOF)
	{
		if(ch == '\n')
		{
			count++;
		}
	}

	rewind(fp);
	
	if (n == 1)
	{
		for ( i = 0; i < count; i++)
		{
			sf = fscanf(fp, "%lf", &tmp1);

			ary[i][0] = tmp1;
		}
	}

	if (n == 2)
	{
		for ( i = 0; i < count; i++)
		{
			sf = fscanf(fp, "%lf%lf", &tmp1, &tmp2);

			ary[i][0] = tmp1;
			ary[i][1] = tmp2;
		}
	}

	if (n == 3)
	{
		for ( i = 0; i < count; i++)
		{
			sf = fscanf(fp, "%lf%lf%lf", &tmp1, &tmp2, &tmp3);
			ary[i][0] = tmp1;
			ary[i][1] = tmp2;
			ary[i][2] = tmp3;
		}
	}

	if (n == 4)
	{
		for ( i = 0; i < count; i++)
		{
			sf = fscanf(fp, "%lf%lf%lf%lf", &tmp1, &tmp2, &tmp3, &tmp4);
			ary[i][0]=tmp1;
			ary[i][1]=tmp2;
			ary[i][2]=tmp3;
			ary[i][3]=tmp4;
		}
	}

	if (n == 5)
	{
		for ( i = 0; i < count; i++)
		{
			sf = fscanf(fp, "%lf%lf%lf%lf%lf", &tmp1, &tmp2, &tmp3, &tmp4, &tmp5);
			ary[i][0]=tmp1;
			ary[i][1]=tmp2;
			ary[i][2]=tmp3;
			ary[i][3]=tmp4;
			ary[i][4]=tmp5;
		}
	}

	rewind(fp);
	fclose(fp);
}

void Get_kernel(int Ns, int Nr, int Nt, double dt)
{
	int nm, nm1, nm2, magnfi;
	long i, j, fs, count, count1, count2;
	char flnm1[100], flnm2[100], FF[100];
	double dx, dy, flgd1, flgd2, flgd3, flgd4, flgd5, **arrd2;
	FILE *fp;

	magnfi = 5;

	dx = 1.0;  dy = 1.0;   nm = Nt / magnfi;
//	dt = 10.0*dt;
	nm1 = 0;   nm2 = nm + 1;

	sprintf(flnm1,"../Output/Snapshot_syn/source1/Ux1.dat");

	count = count_file(flnm1);
	arrd2 = Arraydouble2(count, 5);
	rdfile(flnm1, count, 5, arrd2);

	for (i = 0; i < count; i++)
	{
		arrd2[i][2] = 0.0;
		arrd2[i][3] = 0.0;
		arrd2[i][4] = 0.0;
	}

	for (i = 0; i < nm; i++)
	{
		nm1++, nm2--;

		sprintf(flnm1, "../Output/Snapshot_bianry_adj/source%d/receiver%d/Uu%d.bin", Ns, Nr, nm1);
		sprintf(flnm2, "../Output/Snapshot_bianry_syn/source%d/Uu%d.bin", Ns, nm2);

		multi_bin(flnm1, flnm2, count, arrd2);
	}

	sprintf(flnm1, "../Output/Kernel/Kernel_vel_s%d_r%d.dat", Ns, Nr);
	sprintf(FF, "cd ../Output/Snapshot_bianry_adj/source%d/receiver%d; rm *", Ns, Nr);

	if ((fp = fopen(flnm1, "w")) == NULL)
	{
		printf("Open filed! -- %s\n", flnm1), exit(1);
	}

	for (i = 0; i < count; i++)
	{
		arrd2[i][2] = (dx*dy)*dt*magnfi*arrd2[i][2];
		arrd2[i][3] = (dx*dy)*dt*magnfi*arrd2[i][3];

		fprintf(fp, "%.8e%20.8e%20.8e\n", arrd2[i][0], arrd2[i][1], arrd2[i][2]), fflush(fp);
	}

	fs = system(FF);
	fclose(fp);
	Freedouble2(count, arrd2);
}

void multi_bin(char *fl1, char *fl2, long count, double **arrd2)
{
	int nm1, nm2, sf, fs;
	long i, j, count1, count2;
	float a[4], b[4];
	double cc1, cc11, cc12, cc2, cc21, cc22;
	FILE *fp1, *fp2;

	sf = sizeof(float);

	fp1 = fopen(fl1, "rb"); fp2 = fopen(fl2, "rb");
	if (fp1 == NULL || fp2 == NULL)
	{
		printf("Can't open Uu.bin file --\n %s\n%s\n", fl1, fl2), exit(1);
	}

	for (i = 0; i < count; i++)
	{
		fs = fread(a, sf, 4, fp1);
		fs = fread(b, sf, 4, fp2);

		cc11 = a[0] + a[3];
		cc12 = b[0] + b[3];
		cc21 = a[2] + a[3];
		cc22 = b[2] + b[3];

		cc1 = cc11*cc12;
		cc2 = cc21*cc22;

		arrd2[i][2] += cc1;
		arrd2[i][3] += cc2;
	}

	fclose(fp1), fclose(fp2);
}

void Average_rank(int rank, int size, int Ns, int *t1, int *t2)
{
	int i, es, rs, ss, s1, s2;

	es = Ns / size;
	rs = Ns - es*size;
	ss = 0; s1 = 0;

	for(i = 0; i <= rank; i++)
	{
		ss += es;
		if (i < rs)
		{
			ss++;
		}
		if (i == rank - 1)
		{
			s1 = ss;
		}
	}

	s2 = ss;
	*t1 = s1;
	*t2 = s2;
}