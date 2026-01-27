#include "Common.h"
#include "Rfile.h"

int main(int argc,char *argv[])
{
    int i, j, k, l, flg1, flg2, id1, id2, numNr, numNs, numData;
    double tt, error1, error2, error3, **arrayError;
    char fileError[100], fileErrorDD[100];
    FILE *fp1, *fp2;

    sprintf(fileError, "../Output/Error/Error.txt");
    sprintf(fileErrorDD, "../Output/Error/Error_dd.txt");

    numData = Count_file(fileError);
    arrayError = Arraydouble2(numData, 4);
    Read_file(fileError, numData, 4, arrayError);

    numNr = 1;
    numNs = 1;

    for (i = 0; i < numData; i++) {
        if (arrayError[i][0] > numNs) {
           numNs = arrayError[i][0];
        }
        if (arrayError[i][1] > numNr) {
           numNr = arrayError[i][1];
        }
    }

    fp2 = Open_file(fileErrorDD, "w");
    for (i = 1; i <= numNs; i++) {
        for (j = 1; j <= numNr; j++) {
            id1 = j;

            for (k = 0; k <= numNr; k++)  {
                if (k > j)  {
                    id2 = k;
                    flg1 = 0;
                    flg2 = 0;

                    for (l = 0; l < numData; l++) {
                        if (flg1 == 0 && (fabs(arrayError[l][0] - i) < 0.01) && (fabs(arrayError[l][1] - id1) < 0.01) ) {
                            error1 = arrayError[l][2];
                            flg1 = 1;
                        }
                        if (flg2 == 0 && (fabs(arrayError[l][0] - i) < 0.01) && (fabs(arrayError[l][1] - id2) < 0.01) )
                        {
                            error2 = arrayError[l][2];
                            flg2 = 1;
                        }

                        if (flg1 != 0 && flg2 != 0)  {
                            error3 = error1 - error2;
                            fprintf(fp2, "%d %d %d %.16e\n", i, id1, id2, error3);
                            break;
                        }
                    }
                }
            }
        }
    }
    fclose(fp2);
    Freedouble2(numData, arrayError);
    return 0;
}