#include "Common.h"

void CrossCrelationX(char *F1, char *F2, double x1, double y1, double x2, double y2,
	double *time, double *ratio);
void CrossCrelationZ(char *F1, char *F2, double x1, double y1, double x2, double y2,
	double *time, double *ratio);
void SourecAdj(int numData, double timeResx, double timeResz, double **dataOri,
	double **dataDer1, double **dataDer2, double **sourceAdj);
void FirstOrderDerivative(int numData, double **dataOri,  double **dataDer1);
void SecondOrderDerivative(int numData, double **dataOri,  double **dataDer2);

int main(int argc,char *argv[])
{
	int i, j, flg, numDt, numData, Ns, Nr;
	double thodCC, x1, y1, z1, x2, y2, z2, timeX, timeZ, ratioX, ratioZ, error, dt, T, simData;
	double  **dataOri, **dataDer1, **dataDer2, **sourceAdj, **arrayData;
	char fileData[100], fileError[100], fileSyn[100], fileObs[100], fileSor[100];
	FILE *fp1, *fp2, *fp3;

	dt = 0.01;
	thodCC = 0.90;

	sprintf(fileData, "../Input/arrival_p.txt");
	sprintf(fileError, "../Output/Error/Error.txt");

	numData = Count_file(fileData);
	arrayData = Arraydouble2(numData, 3);
	Read_file(fileData, numData, 3, arrayData);

	flg = 0;
	fp1 = Open_file(fileError, "w");
	for (i = 0; i < numData; i++) {
		Ns = arrayData[i][0];
		Nr = arrayData[i][1];
		T = arrayData[i][2];
		numDt = T / dt;

		sprintf(fileObs, "../Obs/Receivers/source%d/D%d", Ns, Nr);
		sprintf(fileSyn, "../Output/Receivers_syn/source%d/D%d", Ns, Nr);

		x1 = 0.1;  y1 = T;
		x2 = x1;   y2 = y1;

		//CrossCrelationX(fileObs, fileSyn, x1, y1, x2, y2, &timeX, &ratioX);
		CrossCrelationZ(fileObs, fileSyn, x1, y1, x2, y2, &timeZ, &ratioZ);
		error = timeZ*timeZ/2.0;
		//simData = fabs(timeZ/timeX);

		if (ratioZ >= thodCC)  {
			fprintf(fp1,"%d %5d %20.12e %15.4e\n", Ns, Nr, error, ratioZ);
			flg++;
		}
	}
	printf("%10s %d %d\n", " ", numData, flg);
	fclose(fp1);
	Freedouble2(numData, arrayData);
	return 0;
	////////////////////////////
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
		Nrx += dataOri[i][1]*dataDer2[i][1];
		Nrz += dataOri[i][2]*dataDer2[i][2];
	}
	Nrx = timeResx/(Nrx*dt);
	Nrz = timeResz/(Nrz*dt);
	// Nrx = 1.0/(Nrx*dt);
	// Nrz = 1.0/(Nrz*dt);


	for (i = 0; i < numData; i++)  {
		sourceAdj[i][0] = dataDer1[i][0];
		sourceAdj[i][1] = -dataDer1[numData-i-1][1]*Nrx;
		sourceAdj[i][2] = -dataDer1[numData-i-1][2]*Nrz;
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

	fp1 = fopen(F1, "r");
	fp2 = fopen(F2, "r");

	if (fp1 == NULL || fp2 == NULL) {
		printf("Can not open waveform data!! --CrossCrelation\n"); exit(-1);
	}

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

	if (fp1 == NULL || fp2 == NULL) {
		printf("Can not open waveform data!! --CrossCrelation\n"); exit(-1);
	}
	

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
