#include "SEM.h"

int getKernel(int Ns, int Nr, int TimeIntv, int numDt, long numNode, long numUu, float dt, float *arrayUusyn, float *arrayUuadj)
{
	int i, j, k, l, count, sizei, sizel, sizef, nm, nm1, nm2, fs;
	int numPdata, numSoomtk, numKenl, numVp2d;
	long id1, id2;
	// short N, Pflag, ndiv, nsnap, nstep, sflag, stype, sf;
	float a[4], b[4];
	// double dt, F, Vp, Layx, Layz, R0, Lscale;
	double cc1, cc2, cc11, cc12, cc21, cc22, souDep;
	double **arrayPdata, **arrayMesh;
	char filePdata[100], fileKernlf[100], fileKernli[100], fileMesh[100], fileKernl2d[100];
	char FF[210], moment[100];
	FILE *fp1, *fp2, *fp3;

	/// Load prameters ///
	//////
	sprintf(fileMesh, "../Output/Snapshot_syn/source1/Ux1.dat");
	sprintf(fileKernl2d, "../Output/Kernel/Kernel_vel_s%d_r%d.dat", Ns, Nr);
	souDep = 18000.0;

	count = count_file(fileMesh);
	arrayMesh = Arraydouble2(count, 5);
	rdfile(fileMesh, count, 5, arrayMesh);

	//////
	/// Load syn data ///
	//////
	if (numNode != count) {
		printf("The syn and adj data is not marching!\n"); return 1;
	}

	for (i = 0; i < count; i++) {
		arrayMesh[i][2] = 0.0;
		arrayMesh[i][3] = 0.0;
		arrayMesh[i][4] = 0.0;
	}
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
				a[l] = arrayUusyn[id1];
				b[l] = arrayUuadj[id2];

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
		// printf("cc1\n"), sleep(1000);
		nm1++;
		nm2--;
	}

	fp3 = Open_file(fileKernl2d, "w");
	for (j = 0; j < count; j++)
	{
		arrayMesh[j][2] *= dt*TimeIntv;
		arrayMesh[j][3] *= dt*TimeIntv;

		if (arrayMesh[j][1] < -souDep) {
			arrayMesh[j][2] = -arrayMesh[j][2];
		}
		
		fprintf(fp3, "%.3f %.3f %.8e\n", arrayMesh[j][0], arrayMesh[j][1], arrayMesh[j][2]);
	}
	fclose(fp3);

	Freedouble2(count, arrayMesh);
	return 0;
}

void Multi_bin(const char *flnm1, const char *flnm2, long count, double **arrd2)
{
	int i, j, count1, count2, nm1, nm2, sf, fs;
	float a[4], b[4];
	double cc1, cc11, cc12, cc2, cc21, cc22;
	FILE *fp1, *fp2;

	sf = sizeof(float);

	fp1 = fopen(flnm1, "rb"); fp2 = fopen(flnm2, "rb");
	if (fp1 == NULL || fp2 == NULL)
	{
		printf("Can't open Uu.bin file --\n %s\n%s\n", flnm1, flnm2), exit(1);
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