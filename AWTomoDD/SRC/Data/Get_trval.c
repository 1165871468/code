#include "Common.h"
#include "Rfile.h"

int main()
{
	int i, j;
	int numX, numY, Ns, Nr, numStep, numCore;
	double time;
	double x1, x2, y1, y2, maxPer, dampFct, soomFct, extFct;
	char fileTimes[100];
	
	FILE *fp1;

	Read_parInver2d(&numX, &numY, &Ns, &Nr, &numStep, &numCore, &x1, &x2, &y1, &y2, &maxPer, &dampFct, &soomFct, &extFct);

	time = 20.0;

	sprintf(fileTimes, "../Input/arrival_p.txt");

	fp1 = Open_file(fileTimes, "w");
	
	for (i = 1; i <= Ns; i++) {
		for (j = 1; j <= Nr; j++) {
			fprintf(fp1, "%d %5d %10.2f\n", i, j, time);
		}
	}

	fclose(fp1);
	return 0;
}
