#include "SEM.h"
#include <stdio.h>
void Build_element_centers(short N2, long elem, long **Elem_array, double **node_array, double *elem_center_x, double *elem_center_z);
long Find_nearest_id_medium_fast( long elem, double *elem_center_x, double *elem_center_z, double x0, double z0);


void getKernel(int Ns, int Nr, int TimeIntv, int numDt, long numNode, long numUu, float dt,
	float *arrayUsyn, float *arrayUadj, double **node_array)
{
	int i, j, k, l, count, sizei, sizel, sizef, nm, nm1, nm2, fs;
	int numPdata, numSoomtk, numKenl, numVp2d;
	long id1, id2;
	float a[4], b[4];
	double cc1, cc2, cc11, cc12, cc21, cc22, souDep;
	double **arrayPdata, **arrayMesh;
	char filePdata[100], fileKernlf[100], fileKernli[100], fileKernl2d[100];
	char FF[210], moment[100];
	FILE *fp1, *fp2, *fp3;

	/// Load prameters ///
	//////
	sprintf(fileKernl2d, "../Output/Kernel/Kernel_vel_s%d_r%d.dat", Ns, Nr);

	count = numNode;
	arrayMesh = Arraydouble2(count, 5);
	for (i = 0; i < count; i++) {
		arrayMesh[i][0] = node_array[i][0];
		arrayMesh[i][1] = -1.0*node_array[i][1];
		arrayMesh[i][2] = 0.0;
		arrayMesh[i][3] = 0.0;
		arrayMesh[i][4] = 0.0;
	}
	///
	/// Get kernel ///
	nm1 = numDt;
	nm = nm1;
	nm1 = 0;
	nm2 = nm - 1;

	for (j = 0; j < nm; j++) {
		id1 = nm1*4*count;
		id2 = nm2*4*count;

		for (k = 0; k < count; k++) {
			for (l = 0; l < 4; l++) {
				a[l] = arrayUsyn[id1];
				b[l] = arrayUadj[id2];

				id1++;
				id2++;
			}

			cc11 = a[0] + a[3];
			cc12 = b[0] + b[3];
			cc21 = a[2] + a[3];
			cc22 = b[2] + b[3];

			cc1 = cc11*cc12;
			cc2 = cc21*cc22;
			arrayMesh[k][2] += cc1;
			arrayMesh[k][3] += cc2;
		}

		nm1++;
		nm2--;
	}

	fp3 = Open_file(fileKernl2d, "w");
	for (j = 0; j < count; j++) {
		arrayMesh[j][2] *= dt*TimeIntv;
		arrayMesh[j][3] *= dt*TimeIntv;
		
		fprintf(fp3, "%.3f %.3f %.8e\n", arrayMesh[j][0], arrayMesh[j][1], arrayMesh[j][2]);
	}
	fclose(fp3);

	Freedouble2(count, arrayMesh);
}

void Get_kernel_pho(int Ns, int Nr, int TimeIntv, int numDt, int N2, long numNode, long numUu, long elem, long **Elem_array,
	float dt, float *arrayUsyn, float *arrayUadj, float *rho, float *lamda, float *mu, double **node_array)
{
   int i, j, k, l, count, sizei, sizel, sizef, nm, nm1, nm2, fs;
   int numPdata, numSoomtk, numKenl, numVp2d;
   long id1, id2, id_nearest;
   float a[4], b[4];
   double cc1, cc2, cc11, cc12, cc21, cc22, souDep;
   double *elem_center_x, *elem_center_z, **arrayPdata, **arrayMesh;
   char filePdata[100], fileKernlf[100], fileKernli[100], fileKernl2d[100];
   char FF[210], moment[100];
   FILE *fp1, *fp2, *fp3;

   /// allocate memory for element centers ///
   elem_center_x = Arraydouble1(elem);
   elem_center_z = Arraydouble1(elem);

   Build_element_centers(N2, elem, Elem_array, node_array, elem_center_x, elem_center_z);
   if (elem_center_x == NULL || elem_center_z == NULL)	{
	   printf("Cannot allocate elem center arrays.\n"), exit(1);
   }
   //////
   /// Load prameters ///
   //////
   sprintf(fileKernl2d, "../Output/Kernel/Kernel_vel_s%d_r%d.dat", Ns, Nr);

   count = numNode;
   arrayMesh = Arraydouble2(count, 7);
   for (i = 0; i < count; i++) {
	   arrayMesh[i][0] = node_array[i][0];
	   arrayMesh[i][1] = -1.0*node_array[i][1];
	   arrayMesh[i][2] = 0.0;
	   arrayMesh[i][3] = 0.0;

	   id_nearest = Find_nearest_id_medium_fast(elem, elem_center_x, elem_center_z, node_array[i][0], node_array[i][1]);

	   arrayMesh[i][4] = rho[id_nearest]; // rho
	   arrayMesh[i][5] = sqrt((lamda[id_nearest] + 2.0*mu[id_nearest]) / rho[id_nearest]); // Vp
	   arrayMesh[i][6] = sqrt(mu[id_nearest]/rho[id_nearest]); // Vs
   }
   ///
   /// Get kernel ///
   nm = numDt;
   nm1 = 0;
   nm2 = nm - 1;

   for (j = 0; j < nm; j++) {
	   id1 = nm1*4*count;
	   id2 = nm2*4*count;

	   for (k = 0; k < count; k++) {
		   for (l = 0; l < 4; l++) {
			   a[l] = arrayUsyn[id1];
			   b[l] = arrayUadj[id2];

			   id1++;
			   id2++;
		   }

		   cc11 = a[0] + a[3];
		   cc12 = b[0] + b[3];
		   cc21 = a[2] + a[3];
		   cc22 = b[2] + b[3];

		   cc1 = cc11*cc12;
		   cc2 = cc21*cc22;
		   arrayMesh[k][2] += cc1;
		   arrayMesh[k][3] += cc2;
	   }

	   nm1++;
	   nm2--;
   }

   fp3 = Open_file(fileKernl2d, "w");
   for (j = 0; j < count; j++) {
	   arrayMesh[j][2] = arrayMesh[j][2] * dt*TimeIntv * (2.0*arrayMesh[j][4]*arrayMesh[j][5]*arrayMesh[j][5]);
	   arrayMesh[j][3] = arrayMesh[j][3] * dt*TimeIntv * (2.0*arrayMesh[j][4]*arrayMesh[j][5]*arrayMesh[j][5]);
	   
	   fprintf(fp3, "%.3f %.3f %.8e\n", arrayMesh[j][0], arrayMesh[j][1], arrayMesh[j][2]);
   }
   fclose(fp3);

   free(elem_center_x);
   free(elem_center_z);
   Freedouble2(count, arrayMesh);
}

void Multi_bin(const char *flnm1, const char *flnm2, long count, double **arrd2)
{
	int i, j, count1, count2, nm1, nm2, sf, fs;
	float a[4], b[4];
	double cc1, cc11, cc12, cc2, cc21, cc22;
	FILE *fp1, *fp2;

	sf = sizeof(float);

	fp1 = fopen(flnm1, "rb"); fp2 = fopen(flnm2, "rb");
	if (fp1 == NULL || fp2 == NULL) {
		printf("Can't open Uu.bin file --\n %s\n%s\n", flnm1, flnm2); exit(1);
	}

	for (i = 0; i < count; i++) {
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

void Build_element_centers(short N2, long elem, long **Elem_array, double **node_array, double *elem_center_x, double *elem_center_z)
{
    int num, n1, n2, n3, n4;
    long i, id1, id2, id3, id4;

    num = floor(sqrt((double)N2) + 0.5);

    if (num * num != N2) {
        printf("Error: N2 = %d is not a square number.\n", N2), exit(1);
    }

    n1 = 0;
    n2 = num * (num - 1);
    n3 = num * num - 1;
    n4 = num - 1;

    for (i = 0; i < elem; i++) {
        id1 = Elem_array[i][n1] - 1;
        id2 = Elem_array[i][n2] - 1;
        id3 = Elem_array[i][n3] - 1;
        id4 = Elem_array[i][n4] - 1;

        elem_center_x[i] = 0.25*(
            node_array[id1][0] +
            node_array[id2][0] +
            node_array[id3][0] +
            node_array[id4][0]
        );

        elem_center_z[i] = 0.25*(
            node_array[id1][1] +
            node_array[id2][1] +
            node_array[id3][1] +
            node_array[id4][1]
        );
    }
}

long Find_nearest_id_medium_fast( long elem, double *elem_center_x, double *elem_center_z, double x0, double z0)
{
    long i, id_min;
    double dx, dz, dist2, dist2_min;

    id_min = 0;
	dist2_min = INFINITY;

    for (i = 0; i < elem; i++)  {
        dx = elem_center_x[i] - x0;
        dz = elem_center_z[i] - z0;

        dist2 = dx*dx + dz*dz;

        if (dist2 < dist2_min)  {
            dist2_min = dist2;
            id_min = i;
        }
    }

    return id_min;
}