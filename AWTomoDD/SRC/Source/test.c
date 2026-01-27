#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#define pi acos(-1.0)
#define earth_r 6371393.0

int main(int argc,char *argv[]) {
    int i, j, k, numData;
    double w, t0, t1, t2, tt, dt, f1, f2, f3, f4;
    char file1[100], file2[100];
    FILE *fp1, *fp2;

    sprintf(file1, "D-o.txt");
    sprintf(file2, "D-s.txt");

    fp1 = fopen(file1, "w");
    fp2 = fopen(file2, "w");

    t0 = 5.0;
    t1 = 1.3;
    t2 = 2.0;
    dt = 0.01;
    numData = t0/dt + 1;

    for (i = 0; i < numData; i++)  {
        tt = i*dt;

        f1 = 0.0;
        f2 = 0.0;
        if (tt >= t1 && tt <= t1+1)  {
            f1 = sin((tt - t1)*2.0*pi);
        }
        if (tt >= t2 && tt <= t2+1)  {
            f2 = sin((tt - t2)*2.0*pi);
        }

        fprintf(fp1, "%.4f %.4f %.4f\n", tt, f1, f1);
        fprintf(fp2, "%.4f %.4f %.4f\n", tt, f2, f2);
    }

    fclose(fp1);
    fclose(fp2);

    return 0;
}