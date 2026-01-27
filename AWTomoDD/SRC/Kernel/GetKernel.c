#include "Common.h"
#include "Rfile.h"
void multi_bin(char *fl1, char *fl2, long count, double **arrd2);

int main(int argc,char *argv[])
{
	int i, j, k, nm, nm1, nm2, count, count1, count2, magnfi, numTime, fs;
	int nt, ns, nr, Nt, Ns, Nr, numX, numZ, numStep, numCore;
	double x1, x2, z1, z2, valPer, valDam, valSmo;
	double T, dx, dy, dt, souDep, **arrd2, **arrayData;
	char flnm1[100], flnm2[100], FF[100];
	FILE *fp1;

	dt = 0.01;
	dx = 1.0;
	dy = 1.0;
	magnfi = 5;       

	sprintf(flnm1, "../Output/Snapshot_syn/source1/Ux1.dat");
	sprintf(flnm2, "../Input/arrival_p.txt");

	count1 = Count_file(flnm1);
	numTime = Count_file(flnm2);

	arrd2 = Arraydouble2(count1, 5);
	arrayData = Arraydouble2(numTime, 3);

	Read_file(flnm1, count1, 5, arrd2);
	Read_file(flnm2, numTime, 3, arrayData);


	for (i = 0; i < numTime; i++) {
		ns = arrayData[i][0];
		nr = arrayData[i][1];
		T = arrayData[i][2];
		Nt = T / dt;
		nm = Nt / magnfi;
		souDep = 50000.0;

		// printf("%d %d %d %d\n", ns, nr, nm, count1), sleep(1000);

		nm1 = 0;
		nm2 = nm + 1;

		for (j = 0; j < count1; j++) {
			arrd2[j][2] = 0.0;
			arrd2[j][3] = 0.0;
			arrd2[j][4] = 0.0;
		}
		for (j = 0; j < nm; j++) {
			nm1++, nm2--;
			sprintf(flnm1, "../Output/Snapshot_bianry_adj/source%d/Uu%d_%d.bin", ns, nr, nm1);
			sprintf(flnm2, "../Output/Snapshot_bianry_syn/source%d/Uu%d_%d.bin", ns, nr, nm2);
			multi_bin(flnm1, flnm2, count1, arrd2);
		}

		sprintf(flnm1, "../Output/Kernel/Kernel_vel_s%d_r%d.dat", ns, nr);
		
		fp1 = Open_file(flnm1, "w");
		for (j = 0; j < count1; j++) {
			arrd2[j][2] = (dx*dy)*dt*magnfi*arrd2[j][2];
			arrd2[j][3] = (dx*dy)*dt*magnfi*arrd2[j][3];

			if (arrd2[j][1] < -souDep) {
				arrd2[j][2] = -arrd2[j][2];
			}

			fprintf(fp1, "%.8e%20.8e%20.8e\n", arrd2[j][0], arrd2[j][1], arrd2[j][2]);
		}
		fclose(fp1);
	}

	Freedouble2(count1, arrd2);
	Freedouble2(numTime, arrayData);
	return 0;
}

void multi_bin(char *fl1, char *fl2, long count, double **arrd2)
{
	int nm1, nm2, sf, fs;
	long i, j, count1, count2;
	float a[4], b[4];
	double cc1, cc11, cc12, cc2, cc21, cc22;
	FILE *fp1, *fp2;

	sf = sizeof(float);

	fp1 = Open_file(fl1, "rb");
	fp2 = Open_file(fl2, "rb");

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
