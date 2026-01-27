#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#define PI acos(-1.0)

int Count_file(const char *flnm);
double *Arraydouble1(long Nx);
double **Arraydouble2(long Nx,long Nz);
void CrossCrelationX(char *F1, char *F2, double x1, double y1, double x2, double y2,
	double *time, double *ratio);
void CrossCrelationZ(char *F1, char *F2, double x1, double y1, double x2, double y2,
	double *time, double *ratio);
void SourecAdj(int numData, double timeResx, double timeResz, double **dataOri,
	double **dataDer1, double **dataDer2, double **sourceAdj);
void FirstOrderDerivative(int numData, double **dataOri,  double **dataDer1);
void Freedouble2(long nx, double **Array);
void Read_file(char *flnm, int row, int col, double **arrayData);
void SecondOrderDerivative(int numData, double **dataOri,  double **dataDer2);
FILE *Open_file(const char *flnm, const char *mode);

int main(int argc,char *argv[])
{
	int i, j, numData;
	char FF[500], F1[100], F2[100];
	double yy, zz, x1, y1, z1, x2, y2, z2, result1, result2, timeX, timeZ, ratioX, ratioZ;
	double sx1, sy1, sx2, sy2, **dataOri, **dataDer1, **dataDer2, **sourceAdj;

	FILE *fp1, *fp2, *fp3;

	sprintf(F1, "./D-o.txt");
	sprintf(F2, "./D-s.txt");
	x1 = 0.0;    y1 = 5.0;
	x2 = 0.0;    y2 = 5.0;

	sprintf(F1, "./D5-o");
	sprintf(F2, "./D5-s");
	x1 = 1.0;    y1 = 11.0;
	x2 = 1.0;    y2 = 11.0;

	numData = Count_file(F1);
	dataOri = Arraydouble2(numData, 3);
	dataDer1 = Arraydouble2(numData, 3);
	dataDer2 = Arraydouble2(numData, 3);
	sourceAdj = Arraydouble2(numData, 3);

	Read_file(F1, numData, 3, dataOri);

	sprintf(FF, "%s-t", F1);

	CrossCrelationX(F1, F2, x1, y1, x2, y2, &timeX, &ratioX);

	CrossCrelationZ(F1, F2, x1, y1, x2, y2, &timeZ, &ratioZ);

	FirstOrderDerivative(numData, dataOri,  dataDer1);

	SecondOrderDerivative(numData, dataOri,  dataDer2);

	SourecAdj(numData, timeX, timeZ, dataOri, dataDer1,  dataDer2,  sourceAdj);

	fp1 = Open_file(FF, "w");
	for (i = 0; i < numData; i++)  {
		fprintf(fp1, "%.4f %.4f %.4f\n", sourceAdj[i][0], sourceAdj[i][1], sourceAdj[i][2]);
	}
	fclose(fp1);

	Freedouble2(numData, dataOri);
	Freedouble2(numData, dataDer1);
	Freedouble2(numData, dataDer2);
	Freedouble2(numData, sourceAdj);
	return 0;
}

void SourecAdj(int numData, double timeResx, double timeResz, double **dataOri,
	double **dataDer1, double **dataDer2, double **sourceAdj)
{
	int i, j;
	double dt, Nrx, Nrz;

	dt = dataOri[1][0] - dataOri[0][0];

	Nrx = 0.0;
	Nrz = 0.0;
	for (i = 0; i < numData; i++)  {
		Nrx += dataOri[i][1]*dataDer2[i][1]*dt;
		Nrz += dataOri[i][2]*dataDer2[i][2]*dt;
	}
	Nrx = timeResx/Nrx;
	Nrz = timeResz/Nrz;

	printf("%.4e  %.4e\n", Nrx, Nrz);

	for (i = 0; i < numData; i++)  {
		sourceAdj[i][0] = dataDer1[numData-i-1][0];
		sourceAdj[i][1] = dataDer1[numData-i-1][1]*Nrx;
		sourceAdj[i][2] = dataDer1[numData-i-1][2]*Nrz;
	}
}

void FirstOrderDerivative(int numData, double **dataOri,  double **dataDer1)
{
	int i, j;
	double dx, fy1, fy2, fz1, fz2, fy, fz;
	FILE *fp1;

	dx = dataOri[1][0] - dataOri[0][0];

	for (i = 0; i < numData; i++)  {
		if (i == 0)  {
			fy1 = dataOri[i][1];
			fy2 = dataOri[i+1][1];
			fz1 = dataOri[i][2];
			fz2 = dataOri[i+1][2];

			fy = (fy2 - fy1)/dx;
			fz = (fz2 - fz1)/dx;
		}
		else if (i == numData-1)  {
			fy1 = dataOri[i-1][1];
			fy2 = dataOri[i][1];
			fz1 = dataOri[i-1][2];
			fz2 = dataOri[i][2];

			fy = (fy2 - fy1)/dx;
			fz = (fz2 - fz1)/dx;
		}
		else  {
			fy1 = dataOri[i-1][1];
			fy2 = dataOri[i+1][1];
			fz1 = dataOri[i-1][2];
			fz2 = dataOri[i+1][2];

			fy = (fy2 - fy1)/(2.0*dx);
			fz = (fz2 - fz1)/(2.0*dx);
		}

		dataDer1[i][0] = dataOri[i][0];
		dataDer1[i][1] = fy;
		dataDer1[i][2] = fz;
	}
}

void SecondOrderDerivative(int numData, double **dataOri,  double **dataDer2)
{
	int i, j;
	double dx, fy1, fy2, fy3, fz1, fz2, fz3, fy, fz;

	dx = dataOri[1][0] - dataOri[0][0];

	for (i = 0; i < numData; i++)  {
		if (i == 0)  {
			fy1 = dataOri[i][1];
			fy2 = dataOri[i+1][1];
			fy3 = dataOri[i+2][1];
			fz1 = dataOri[i][2];
			fz2 = dataOri[i+1][2];
			fz3 = dataOri[i+2][2];

			fy = (fy3 - 2.0*fy2 + fy1)/(dx*dx);
			fz = (fz3 - 2.0*fz2 + fz1)/(dx*dx);
		}
		else if (i == numData-1)  {
			fy1 = dataOri[i-2][1];
			fy2 = dataOri[i-1][1];
			fy3 = dataOri[i][1];
			fz1 = dataOri[i-2][2];
			fz2 = dataOri[i-1][2];
			fz3 = dataOri[i][2];

			fy = (fy3 - 2.0*fy2 + fy1)/(dx*dx);
			fz = (fz3 - 2.0*fz2 + fz1)/(dx*dx);
		}
		else  {
			fy1 = dataOri[i-1][1];
			fy2 = dataOri[i][1];
			fy3 = dataOri[i+1][1];
			fz1 = dataOri[i-1][2];
			fz2 = dataOri[i][2];
			fz3 = dataOri[i+1][2];

			fy = (fy3 - 2.0*fy2 + fy1)/(dx*dx);
			fz = (fz3 - 2.0*fz2 + fz1)/(dx*dx);
		}

		dataDer2[i][0] = dataOri[i][0];
		dataDer2[i][1] = fy;
		dataDer2[i][2] = fz;
	}
}


void CrossCrelationX(char *F1, char *F2, double x1, double y1, double x2, double y2, double *time, double *ratio)
{
	int i, j, ij, id, ti, d1, d2, count1, count2, count3, count4, count5, count6, NN1, NN2, sf;
	double xx, yy1, cy1, cy2, M1, M2, tm1, tm2, tot, tt, delta1, delta2;
	double *datax1, *datay1, *datax2, *datay2, *Ny1, *Ny2;
	char *ff, FF[100];
	FILE *fp1, *fp2;

	if(y1-x1 < 0.1 || y2-x2 < 0.1)  {
		printf("The waveform data error!");  exit(0);
	}

	fp1 = fopen(F1,"r");
	fp2 = fopen(F2,"r");

	FF[0] = '\0';
	count1 = 0;
	while((ff = fgets(FF, 90, fp1)) != NULL)  {
		count1++;
	}
	rewind(fp1);

	count2 = 0;
	while((ff = fgets(FF, 90, fp2)) != NULL)  {
		count2++;
	}
	rewind(fp2);

	datax1 = Arraydouble1(count1);
	datay1 = Arraydouble1(count1);
	datax2 = Arraydouble1(count2);
	datay2 = Arraydouble1(count2);

	for(i = 0; i < count1; i++)  {
		// sf=fscanf(fp1, "%lf %lf", &datax1[i], &datay1[i]);
		sf=fscanf(fp1, "%lf %lf %lf", &datax1[i], &datay1[i], &tt);
	}
	fclose(fp1);

	for(i = 0; i < count2; i++)  {
		// sf=fscanf(fp2, "%lf %lf", &datax2[i], &datay2[i]);
		sf=fscanf(fp2, "%lf %lf %lf", &datax2[i], &datay2[i], &tt);
	}
	fclose(fp2);

	for(i = 0; i < count1; i++)  {
		if(datax1[i] <= x1 && datax1[i+1] >= x1)  {
			count3 = i;
		}

		if(datax1[i] <= y1 && datax1[i+1] >= y1)  {
			count4 = i;
		}
	}

	for(i = 0; i < count2; i++)  {
		if(datax2[i] <= x2 && datax2[i+1] >= x2)  {
			count5 = i;
		}

		if(datax2[i] <= y2 && datax2[i+1] >= y2)  {
			count6 = i;
		}
	}

	delta1 = 100.0;
	for(i = count3; i < count4-1; i++)  {
		tm1 = datax1[i+1] - datax1[i];

		if(tm1 < delta1)  {
			delta1=tm1;
		}
	}
	delta2 = 100.0;
	for(i = count5; i < count6-1; i++)  {
		tm1 = datax2[i+1] - datax2[i];

		if(tm1 < delta2)  {
			delta2 = tm1;
		}
	}

	if(delta1 > delta2)  {
		delta1 = delta2;
	}

	NN1 = (int)((y1 - x1)/delta1 + 0.0000001) + 1;
	NN2 = (int)((y2 - x2)/delta1 + 0.0000001) + 1;

	Ny1 = Arraydouble1(NN1);
	Ny2 = Arraydouble1(NN2);

	id = count3;
	for(i = 0; i < NN1; i++)  {
		xx = x1 + i*delta1;

		if(xx >= datax1[id+1])  {
			id++;
		}

		tm1 = datax1[id];
		tm2 = datax1[id+1];

		Ny1[i] = (tm2 - xx)/(tm2 - tm1)*datay1[id] + (xx - tm1)/(tm2 - tm1)*datay1[id+1];
	}
	free(datax1);
	free(datay1);

	d1 = (int)(0.05*NN1);
	if(d1 == 0)  {
		d1 = 1;
	}

	for(i = 0; i < d1; i++)  {
		Ny1[i] = 1.0*i/d1*Ny1[i];
		Ny1[i+NN1-d1] = 1.0*(d1 - i - 1)/d1*Ny1[i + NN1 - d1];
	}

	M1 = 0.0;
	for(i = 0; i < NN1; i++)  {
		if(fabs(Ny1[i]) > M1)  {
			M1 = fabs(Ny1[i]);
		}
	}

	for(i = 0; i < NN1; i++)  {
		Ny1[i] = Ny1[i] / M1;
	}
	
	id = count5;
	for(i = 0; i < NN2; i++)  {
		xx = x2 + i*delta1;

		if(xx >= datax2[id+1])  {
			id++;
		}
		
		tm1 = datax2[id];
		tm2 = datax2[id+1];
		
		Ny2[i] = (tm2 - xx)/(tm2 - tm1)*datay2[id] + (xx - tm1)/(tm2 - tm1)*datay2[id+1];
	}
	free(datax2);
	free(datay2);

	d1 = (int)(0.05*NN2);
	if(d1 == 0)  {
		d1 = 1;
	}

	for(i = 0; i < d1; i++)  {
		Ny2[i] = 1.0*i/d1*Ny2[i];
		Ny2[i+NN2-d1] = 1.0*(d1 - i - 1)/d1*Ny2[i+NN2-d1];
	}

	M1 = 0.0;
	for(i = 0; i < NN2; i++)  {
		if(fabs(Ny2[i]) > M1)  {
			M1 = fabs(Ny2[i]);
		}
	}

	for(i = 0; i < NN2; i++)  {
		Ny2[i] = Ny2[i]/M1;
	}

	ti = 0;
	tt = 0.0;	
	d1 = NN1;
	d2 = NN2;

	cy1 = 0.0;
	for(i = 0; i < NN1; i++)  {
		cy1 += Ny1[i]*Ny1[i];
	}
	cy2 = 0.0;
	for(i = 0; i < NN2; i++)  {
		cy2 += Ny2[i]*Ny2[i];
	}

	cy1 = sqrt(cy1);
	cy2 = sqrt(cy2);

	for(i = -d1; i < d2; i++)  {		
		tot = 0.0;

		for(j = 0; j < NN1; j++)  {
			yy1 = Ny1[j];
			ij = i + j;	

			if(ij >= 0 && ij <= NN2-1)  {
				tot += yy1*Ny2[ij];
			}
		}

		if(tot > tt)  {
			tt = tot;
			ti = i;
		}
	}

	*time = x1 - x2 - ti*delta1;
	*ratio = tt/(cy1*cy2);

	free(Ny1);
	free(Ny2);
}

void CrossCrelationZ(char *F1, char *F2, double x1, double y1, double x2, double y2, double *time, double *ratio)
{
	int i, j, ij, id, ti, d1, d2, count1, count2, count3, count4, count5, count6, NN1, NN2, sf;
	double xx, yy1, cy1, cy2, M1, M2, tm1, tm2, tot, tt, delta1, delta2;
	double *datax1, *datay1, *datax2, *datay2, *Ny1, *Ny2;
	char *ff, FF[100];
	FILE *fp1, *fp2;

	if(y1-x1 < 0.1 || y2-x2 < 0.1)  {
		printf("The waveform data error!");  exit(0);
	}

	fp1 = fopen(F1,"r");
	fp2 = fopen(F2,"r");

	FF[0] = '\0';
	count1 = 0;
	while((ff = fgets(FF, 90, fp1)) != NULL)  {
		count1++;
	}
	rewind(fp1);

	count2 = 0;
	while((ff = fgets(FF, 90, fp2)) != NULL)  {
		count2++;
	}
	rewind(fp2);

	datax1 = Arraydouble1(count1);
	datay1 = Arraydouble1(count1);
	datax2 = Arraydouble1(count2);
	datay2 = Arraydouble1(count2);

	for(i = 0; i < count1; i++)  {
		//sf=fscanf(fp1, "%lf %lf", &datax1[i], &datay1[i]);
		sf=fscanf(fp1, "%lf %lf %lf", &datax1[i], &tt, &datay1[i]);
	}
	fclose(fp1);

	for(i = 0; i < count2; i++)  {
		//sf=fscanf(fp2, "%lf %lf", &datax2[i], &datay2[i]);
		sf=fscanf(fp2, "%lf %lf %lf", &datax2[i], &tt, &datay2[i]);
	}
	fclose(fp2);

	for(i = 0; i < count1; i++)  {
		if(datax1[i] <= x1 && datax1[i+1] >= x1)  {
			count3 = i;
		}

		if(datax1[i] <= y1 && datax1[i+1] >= y1)  {
			count4 = i;
		}
	}

	for(i = 0; i < count2; i++)  {
		if(datax2[i] <= x2 && datax2[i+1] >= x2)  {
			count5 = i;
		}

		if(datax2[i] <= y2 && datax2[i+1] >= y2)  {
			count6 = i;
		}
	}

	delta1 = 100.0;
	for(i = count3; i < count4-1; i++)  {
		tm1 = datax1[i+1] - datax1[i];

		if(tm1 < delta1)  {
			delta1=tm1;
		}
	}
	delta2 = 100.0;
	for(i = count5; i < count6-1; i++)  {
		tm1 = datax2[i+1] - datax2[i];

		if(tm1 < delta2)  {
			delta2 = tm1;
		}
	}

	if(delta1 > delta2)  {
		delta1 = delta2;
	}

	NN1 = (int)((y1 - x1)/delta1 + 0.0000001) + 1;
	NN2 = (int)((y2 - x2)/delta1 + 0.0000001) + 1;

	Ny1 = Arraydouble1(NN1);
	Ny2 = Arraydouble1(NN2);

	id = count3;
	for(i = 0; i < NN1; i++)  {
		xx = x1 + i*delta1;

		if(xx >= datax1[id+1])  {
			id++;
		}

		tm1 = datax1[id];
		tm2 = datax1[id+1];

		Ny1[i] = (tm2 - xx)/(tm2 - tm1)*datay1[id] + (xx - tm1)/(tm2 - tm1)*datay1[id+1];
	}
	free(datax1);
	free(datay1);

	d1 = (int)(0.05*NN1);
	if(d1 == 0)  {
		d1 = 1;
	}

	for(i = 0; i < d1; i++)  {
		Ny1[i] = 1.0*i/d1*Ny1[i];
		Ny1[i+NN1-d1] = 1.0*(d1 - i - 1)/d1*Ny1[i + NN1 - d1];
	}

	M1 = 0.0;
	for(i = 0; i < NN1; i++)  {
		if(fabs(Ny1[i]) > M1)  {
			M1 = fabs(Ny1[i]);
		}
	}

	for(i = 0; i < NN1; i++)  {
		Ny1[i] = Ny1[i] / M1;
	}
	
	id = count5;
	for(i = 0; i < NN2; i++)  {
		xx = x2 + i*delta1;

		if(xx >= datax2[id+1])  {
			id++;
		}
		
		tm1 = datax2[id];
		tm2 = datax2[id+1];
		
		Ny2[i] = (tm2 - xx)/(tm2 - tm1)*datay2[id] + (xx - tm1)/(tm2 - tm1)*datay2[id+1];
	}
	free(datax2);
	free(datay2);

	d1 = (int)(0.05*NN2);
	if(d1 == 0)  {
		d1 = 1;
	}

	for(i = 0; i < d1; i++)  {
		Ny2[i] = 1.0*i/d1*Ny2[i];
		Ny2[i+NN2-d1] = 1.0*(d1 - i - 1)/d1*Ny2[i+NN2-d1];
	}

	M1 = 0.0;
	for(i = 0; i < NN2; i++)  {
		if(fabs(Ny2[i]) > M1)  {
			M1 = fabs(Ny2[i]);
		}
	}

	for(i = 0; i < NN2; i++)  {
		Ny2[i] = Ny2[i]/M1;
	}

	ti = 0;
	tt = 0.0;	
	d1 = NN1;
	d2 = NN2;

	cy1 = 0.0;
	for(i = 0; i < NN1; i++)  {
		cy1 += Ny1[i]*Ny1[i];
	}
	cy2 = 0.0;
	for(i = 0; i < NN2; i++)  {
		cy2 += Ny2[i]*Ny2[i];
	}

	cy1 = sqrt(cy1);
	cy2 = sqrt(cy2);

	for(i = -d1; i < d2; i++)  {		
		tot = 0.0;

		for(j = 0; j < NN1; j++)  {
			yy1 = Ny1[j];
			ij = i + j;	

			if(ij >= 0 && ij <= NN2-1)  {
				tot += yy1*Ny2[ij];
			}
		}

		if(tot > tt)  {
			tt = tot;
			ti = i;
		}
	}

	*time = x1 - x2 - ti*delta1;
	*ratio = tt/(cy1*cy2);

	free(Ny1);
	free(Ny2);
}

FILE *Open_file(const char *flnm, const char *mode)
{
    FILE *fp1;

    if ((fp1 = fopen(flnm, mode)) == NULL)  {
        perror(flnm); exit(1);
    }
    return fp1;
}


int Count_file(const char *flnm)
{
    int count;
    char ch;
    FILE *fp;

    fp = Open_file(flnm, "r");
    count = 0;
    while ((ch = fgetc(fp)) != EOF)  {
        if(ch == '\n')  {
            count++;
        }
    }
    fclose(fp);
    return count;
}

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

void Freedouble2(long nx, double **Array)
{
    long i;

    for(i = 0; i < nx; i++)
    {
        free(Array[i]);
    }
    free(Array);
}

void Read_file(char *flnm, int row, int col, double **arrayData)
{
    int i, j, sf;
    double tt[100];
    char ch, FF[100], *ff;
    FILE *fp;

    fp = Open_file(flnm, "r");
    for (i = 0; i < row; i++)  {
        for (j = 0; j < col; j++)  {
            sf = fscanf(fp, "%lf", &tt[j]);
            arrayData[i][j] = tt[j];
        }
    }
    fclose(fp);
}