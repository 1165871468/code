#include "Common.h"
#include "Rfile.h"

int main(int argc,char *argv[])
{
	int i, j, k, l, c, Ns, Nr, cent, numData;
	int fs, flgi1, flgi2, flgi3;
	int count1, count2, count3;
	double flgd1, flgd2, flgd3, flgd4, flgd5, flgd6, flgd7, flgd8, flgd9, max1, max2, max3, max4, dt, error;
	double **arrayResults, **timePdatd, **waveSyn, **waveObs, **arrayError, **arrd1, **arrd2, **arrd3, **arrayData;
	double T, a1[3] = {0}, a2[3] = {0};
	char fileData[100], fileError[100], fileSyn[100], fileObs[100], fileSou[100], moment[100];
	FILE *fp1, *fp2, *fp3, *fp4, *fp5, *fp6;

	dt = 0.01;
	error = 0.0;

	sprintf(fileData, "../Input/arrival_p.txt");
	sprintf(fileError, "../Output/Error/Error.txt");
	numData = Count_file(fileData);
	arrayData = Arraydouble2(numData, 3);
	Read_file(fileData, numData, 3, arrayData);

	fp1 = Open_file(fileError, "w");

	for (i = 0; i < numData; i++) {
		Ns = arrayData[i][0];
		Nr = arrayData[i][1];
		T = arrayData[i][2];
		count1 = T / dt;

		sprintf(fileSyn, "../Output/Receivers_syn/source%d/D%d", Ns, Nr);
		sprintf(fileObs, "../Obs/Receivers/source%d/D%d", Ns, Nr);
		sprintf(fileSou, "../Output/Source_adj/Source_times%d_r%d.txt", Ns, Nr);

		count2 = Count_file(fileSyn);
		count3 = Count_file(fileObs);

		if (count2 != count1 || count2 != count3) {
			printf("Error, data is not marching! (s%d_r%d): %d %d %d\n", Ns, Nr, count1, count2, count3), exit(1);
		}

		waveSyn = Arraydouble2(count2, 3);
		waveObs = Arraydouble2(count2, 3);
		arrayResults = Arraydouble2(count2, 3);
		Read_file(fileSyn, count2, 3, waveSyn);
		Read_file(fileObs, count2, 3, waveObs);

		max1 = fabs(waveSyn[0][1]);
		max2 = fabs(waveObs[0][1]);
		max3 = fabs(waveSyn[0][2]);
		max4 = fabs(waveObs[0][2]);

		for (k = 0; k < count2; k++) {
			if (max1 < fabs(waveSyn[k][1])) {
				max1 = fabs(waveSyn[k][1]);
			}
			if (max2 < fabs(waveObs[k][1])) {
				max2 = fabs(waveObs[k][1]);
			}
			if (max3 < fabs(waveSyn[k][2])) {
				max3 = fabs(waveSyn[k][2]);
			}
			if (max4 < fabs(waveObs[k][2])) {
				max4 = fabs(waveObs[k][2]);
			}
		}

		error = 0.0;
		for (k = 0; k < count2; k++) {
			arrayResults[k][0] = waveSyn[k][0];
			//arrayResults[k][1] = (waveSyn[k][1]/max1 - waveObs[k][1]/max2)*max1;
			// arrayResults[k][1] = (waveSyn[k][1]/max1 - waveObs[k][1]/max2)*max1;
			arrayResults[k][1] = waveSyn[k][1] - waveObs[k][1];
			// arrayResults[k][2] = (waveSyn[k][2]/max3 - waveObs[k][2]/max4)*max3;
			arrayResults[k][2] = waveSyn[k][2] - waveObs[k][2];
			// arrayResults[k][1] = waveSyn[k][1] - waveObs[k][1];
			// arrayResults[k][2] = waveSyn[k][2] - waveObs[k][2];
			error += dt*(arrayResults[k][2]*arrayResults[k][2]) / 2.0;
		}

		fp2 = Open_file(fileSou, "w");
		for (k = 0; k < count2; k++) {
			l = count2 - k - 1;
			fprintf(fp2,"%.4f %.12e %.12e\n", arrayResults[k][0], arrayResults[l][2], arrayResults[l][2]);
		}
		fclose(fp2);

		fprintf(fp1,"%d %d %.12e\n", Ns, Nr, error);
	}

	fclose(fp1);
	

	Freedouble2(numData, arrayData);
	return 0;
}
