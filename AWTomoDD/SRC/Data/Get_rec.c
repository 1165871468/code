#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

int main()
{
	int i, numNs, numNr;
	double x1, x2, z1, z2, xx, zz, dep, delt;
	FILE *fp;

	numNs = 3;
	numNr = 3;

	x1 = 0.0;      x2 = 100000.0;
	z1 = 0.0;      z2 = 40000.0;
	dep = 0.0;

	delt = (x2 - x1) / (numNr + 1) ;
	fp = fopen("receiver", "w");
	for(i = 1; i < numNr + 1; i++)
	{
		xx = i*delt + x1;
		zz = dep ;
		fprintf(fp, "%.1f%10.1f\n", xx, zz), fflush(fp);
	}
	fclose(fp);

	x1 = 0.0;      x2 = 100000.0;
	z1 = 0.0;      z2 = 40000.0;
	dep = 32000.0;

	delt = (x2 - x1) / (numNs + 1) ;
	fp = fopen("source", "w");
	for(i = 1; i < numNs + 1; i++)
	{
		xx = i*delt + x1;
		zz = dep;
		fprintf(fp, "%.1f%10.1f\n", xx, zz), fflush(fp);
	}
	fclose(fp);

	return 0;
}
