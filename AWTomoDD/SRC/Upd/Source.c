
#include "SEM.h"

int getAdjSource(int Ns, int Nr, double dt, double T)
{
	int i, j, k, l, c, cent, numData;
	int fs, flgi1, flgi2, flgi3;
	int count1, count2, count3;
	double flgd1, flgd2, flgd3, flgd4, flgd5, flgd6, flgd7, flgd8, flgd9, max1, max2, max3, max4, error;
	double **arrayResults, **timePdatd, **waveSyn, **waveObs, **arrayError, **arrd1, **arrd2, **arrd3, **arrayData;
	double a1[3] = {0}, a2[3] = {0};
	char fileData[100], fileError[100], fileSyn[100], fileObs[100], fileSor[100], moment[100];
	FILE *fp1, *fp2, *fp3, *fp4, *fp5, *fp6;
    //////
	sprintf(fileError, "../Output/Error/Error.txt");
    sprintf(fileSyn, "../Output/Receivers_syn/source%d/D%d", Ns, Nr);
    sprintf(fileObs, "../Obs/Receivers/source%d/D%d", Ns, Nr);
    sprintf(fileSor, "../Output/Source_adj/Source_times%d_r%d.txt", Ns, Nr);
	fp1 = Open_file(fileError, "a");
    //////
    count1 = T / dt;
    count2 = count_file(fileSyn);
    count3 = count_file(fileObs);

    if (count2 != count1 || count2 != count3) {
        printf("Error, data is not marching! (s%d_r%d): %d %d %d\n", Ns, Nr, count1, count2, count3), exit(1);
    }

    waveSyn = Arraydouble2(count2, 3);
    waveObs = Arraydouble2(count2, 3);
    arrayResults = Arraydouble2(count2, 3);

    rdfile(fileSyn, count2, 3, waveSyn);
    rdfile(fileObs, count2, 3, waveObs);
    //////
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
    //////
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
    //////
    fp2 = Open_file(fileSor, "w");
    for (k = 0; k < count2; k++) {
        l = count2 - k - 1;
        fprintf(fp2,"%.4f %.12e %.12e\n", arrayResults[k][0], arrayResults[l][2], arrayResults[l][2]);
    }
    fclose(fp2);
    //////
    fprintf(fp1,"%d %d %.12e\n", Ns, Nr, error);
    fclose(fp1);
    //////
    Freedouble2(count2, waveSyn);
    Freedouble2(count2, waveObs);
    Freedouble2(count2, arrayResults);

	return 0;
}
