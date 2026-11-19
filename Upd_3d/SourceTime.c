#include "SEM.h"
#include <stdio.h>

#define PI acos(-1.0)
#define EPS 1.0e-12
#define CROSS_CORRELATION_WINDOW_LEN    4.0     // 互相关窗口长度（秒）

void CrossCrelationX(char *F1, char *F2, double x1, double y1, double x2, double y2,
	double *time, double *ratio);
void CrossCrelationZ(char *F1, char *F2, double x1, double y1, double x2, double y2,
	double *time, double *ratio);
void CrossCorrelation(const char *file_obs, const char *file_syn, double t_start, 
	double t_end, int component, double *time_shift, double *coeff);
void SourecAdj(int numData, double timeResx, double timeResz, double **dataOri,
	double **dataDer1, double **dataDer2, double **sourceAdj);
void FirstOrderDerivative(int numData, double **dataOri,  double **dataDer1);
void SecondOrderDerivative(int numData, double **dataOri,  double **dataDer2);

void GetSourceTimeAdj(int Ns, int Nr, double dt, double T)
{
	int i, j, numData;
	double yy, zz, x1, y1, z1, x2, y2, z2, timeX, timeZ, time_win, ratioX, ratioZ, error;
	double  **dataOri, **dataDer1, **dataDer2, **sourceAdj;
	char fileError[200], fileSyn[200], fileObs[200], fileSor[200];

	FILE *fp1, *fp2, *fp3;

	sprintf(fileError, "../Output/Error/Error.txt");
    sprintf(fileObs, "../Obs/Receivers/source%d/D%d", Ns, Nr);
    sprintf(fileSyn, "../Output/Receivers_syn/source%d/D%d", Ns, Nr);
    sprintf(fileSor, "../Output/Source_adj/Source_times%d_r%d.txt", Ns, Nr);


	time_win = CROSS_CORRELATION_WINDOW_LEN;
	double win_start = T - time_win;
	double win_end   = T;

	if (win_start < 0.0) win_start = 0.0;
	y1 = win_end;  x1 = win_start;
	y2 = win_end;  x2 = win_start;
	///////// Load waveform data //////
	numData = count_file(fileSyn);
	dataOri = Arraydouble2(numData, 3);
	dataDer1 = Arraydouble2(numData, 3);
	dataDer2 = Arraydouble2(numData, 3);
	sourceAdj = Arraydouble2(numData, 3);
	rdfile(fileSyn, numData, 3, dataOri);
	///
	///////// Calculate adjoint source //////
	// CrossCrelationX(fileObs, fileSyn, x1, y1, x2, y2, &timeX, &ratioX);
	// CrossCrelationZ(fileObs, fileSyn, x1, y1, x2, y2, &timeZ, &ratioZ);

	// printf("%.6e %.6e\n", timeX, ratioX);
	// printf("%.6e %.6e\n", timeZ, ratioZ);

    CrossCorrelation(fileObs, fileSyn, win_start, win_end, 1, &timeX, &ratioX);
    CrossCorrelation(fileObs, fileSyn, win_start, win_end, 2, &timeZ, &ratioZ);

	// printf("%.6e %.6e\n", timeX, ratioX);
	// printf("%.6e %.6e\n", timeZ, ratioZ);
	// sleep(1000);

	FirstOrderDerivative(numData, dataOri,  dataDer1);
	SecondOrderDerivative(numData, dataOri,  dataDer2);

	SourecAdj(numData, timeX, timeZ, dataOri, dataDer1,  dataDer2,  sourceAdj);
	///
	///////// Output adjoint source //////
	fp1 = Open_file(fileSor, "w");
	for (i = 0; i < numData; i++)  {
		fprintf(fp1, "%.4f %.4f %.4f\n", sourceAdj[i][0], sourceAdj[i][1], sourceAdj[i][2]);
	}
	fclose(fp1);
	///
	///////// Output error //////
	// fp2 = Open_file(fileError, "a");
	// error = timeZ*timeZ/2.0;
    // fprintf(fp2,"%d %d %.12e\n", Ns, Nr, error);
    // fclose(fp2);
	///
	///////// Free arrays //////
	Freedouble2(numData, dataOri);
	Freedouble2(numData, dataDer1);
	Freedouble2(numData, dataDer2);
	Freedouble2(numData, sourceAdj);
	////////////////////////////
}

void SourecAdj(int numData, double timeResx, double timeResz, double **dataOri,
	double **dataDer1, double **dataDer2, double **sourceAdj)
{
	int i, j;
	double dt, Nrx, Nrz;

	dt = dataOri[1][0] - dataOri[0][0];

	if (numData < 2 ||fabs(dt) < EPS) {
		printf("The waveform data is to short or the time interval is too small !\n");
		exit(-1);
	}

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

// 一阶导数（内部点四阶中心差分，边界点二阶精度）
void FirstOrderDerivative(int numData, double **dataOri, double **dataDer1) {
    double dt = dataOri[1][0] - dataOri[0][0];
    int i;
    
    for (i = 0; i < numData; i++) {
        dataDer1[i][0] = dataOri[i][0];   // 时间列
    }
    
    // 处理 X 分量（索引 1）和 Z 分量（索引 2）
    for (int comp = 1; comp <= 2; comp++) {
        double *y = (double*)malloc(numData * sizeof(double));
        for (i = 0; i < numData; i++) y[i] = dataOri[i][comp];
        
        double *dy = (double*)malloc(numData * sizeof(double));
        
        // 边界点：使用二阶中心差分（单侧边界用前向/后向差分）
        // i = 0: 前向二阶差分
        dy[0] = (-3*y[0] + 4*y[1] - y[2]) / (2*dt);
        // i = 1: 中心差分（二阶精度）
        dy[1] = (y[2] - y[0]) / (2*dt);
        // i = numData-2: 中心差分（二阶精度）
        dy[numData-2] = (y[numData-1] - y[numData-3]) / (2*dt);
        // i = numData-1: 后向二阶差分
        dy[numData-1] = (3*y[numData-1] - 4*y[numData-2] + y[numData-3]) / (2*dt);
        
        // 内部点：四阶中心差分公式
        // f'(i) = (-f(i+2) + 8f(i+1) - 8f(i-1) + f(i-2)) / (12*dt)
        for (i = 2; i < numData-2; i++) {
            dy[i] = (-y[i+2] + 8*y[i+1] - 8*y[i-1] + y[i-2]) / (12.0 * dt);
        }
        
        // 写回 dataDer1
        for (i = 0; i < numData; i++) dataDer1[i][comp] = dy[i];
        
        free(y); free(dy);
    }
}

// 二阶导数（内部点四阶中心差分，边界点二阶精度）
void SecondOrderDerivative(int numData, double **dataOri, double **dataDer2) {
    double dt = dataOri[1][0] - dataOri[0][0];
    double dt2 = dt * dt;
    int i;
    
    for (i = 0; i < numData; i++) {
        dataDer2[i][0] = dataOri[i][0];
    }
    
    for (int comp = 1; comp <= 2; comp++) {
        double *y = (double*)malloc(numData * sizeof(double));
        for (i = 0; i < numData; i++) y[i] = dataOri[i][comp];
        
        double *d2y = (double*)malloc(numData * sizeof(double));
        
        // 边界点：使用二阶精度差分
        // i = 0: 前向三点差分
        d2y[0] = (2*y[0] - 5*y[1] + 4*y[2] - y[3]) / dt2;
        // i = 1: 三点中心差分（二阶精度）
        d2y[1] = (y[0] - 2*y[1] + y[2]) / dt2;
        // i = numData-2: 三点中心差分
        d2y[numData-2] = (y[numData-3] - 2*y[numData-2] + y[numData-1]) / dt2;
        // i = numData-1: 后向三点差分
        d2y[numData-1] = (2*y[numData-1] - 5*y[numData-2] + 4*y[numData-3] - y[numData-4]) / dt2;
        
        // 内部点：四阶中心差分公式
        // f''(i) = (-f(i+2) + 16f(i+1) - 30f(i) + 16f(i-1) - f(i-2)) / (12*dt^2)
        for (i = 2; i < numData-2; i++) {
            d2y[i] = (-y[i+2] + 16*y[i+1] - 30*y[i] + 16*y[i-1] - y[i-2]) / (12.0 * dt2);
        }
        
        for (i = 0; i < numData; i++) dataDer2[i][comp] = d2y[i];
        
        free(y); free(d2y);
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
	if (fp1 == NULL || fp2 == NULL) {
		printf("Open waveform file failed!\n");
		exit(1);
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
	if (count1 < 2 || count2 < 2) {
		printf("Waveform data points are too few!\n");
		exit(1);
	}

	datax1 = Arraydouble1(count1);
	datay1 = Arraydouble1(count1);
	datax2 = Arraydouble1(count2);
	datay2 = Arraydouble1(count2);

	for(i = 0; i < count1; i++)  {
		// sf=fscanf(fp1, "%lf %lf", &datax1[i], &datay1[i]);
		sf=fscanf(fp1, "%lf %lf %lf", &datax1[i], &datay1[i], &tt);
		if (sf != 3) {
			printf("Read waveform data failed: %s\n", F1);
			exit(1);
		}
	}
	fclose(fp1);

	for(i = 0; i < count2; i++)  {
		// sf=fscanf(fp2, "%lf %lf", &datax2[i], &datay2[i]);
		sf=fscanf(fp2, "%lf %lf %lf", &datax2[i], &datay2[i], &tt);
		if (sf != 3) {
			printf("Read waveform data failed: %s\n", F2);
			exit(1);
		}
	}
	fclose(fp2);

	count3 = -1;
	count4 = -1;
	count5 = -1;
	count6 = -1;
	for(i = 0; i < count1 - 1; i++)  {
		if(datax1[i] <= x1 && datax1[i+1] >= x1)  {
			count3 = i;
		}

		if(datax1[i] <= y1 && datax1[i+1] >= y1)  {
			count4 = i;
		}
	}

	for(i = 0; i < count2 - 1; i++)  {
		if(datax2[i] <= x2 && datax2[i+1] >= x2)  {
			count5 = i;
		}

		if(datax2[i] <= y2 && datax2[i+1] >= y2)  {
			count6 = i;
		}
	}
	if (count3 < 0 || count4 < 0 || count5 < 0 || count6 < 0 || count4 <= count3 || count6 <= count5) {
		printf("Window [x,y] is out of waveform range!\n");
		exit(1);
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
	if (delta1 <= EPS) {
		printf("Invalid time interval in waveform data!\n");
		exit(1);
	}

	NN1 = (int)((y1 - x1)/delta1 + 0.0000001) + 1;
	NN2 = (int)((y2 - x2)/delta1 + 0.0000001) + 1;

	Ny1 = Arraydouble1(NN1);
	Ny2 = Arraydouble1(NN2);

	id = count3;
	for(i = 0; i < NN1; i++)  {
		xx = x1 + i*delta1;

		if(id < count1 - 2 && xx >= datax1[id+1])  {
			id++;
		}

		tm1 = datax1[id];
		tm2 = datax1[id+1];

		if (fabs(tm2 - tm1) < EPS) {
			Ny1[i] = datay1[id];
		} else {
			Ny1[i] = (tm2 - xx)/(tm2 - tm1)*datay1[id] + (xx - tm1)/(tm2 - tm1)*datay1[id+1];
		}
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
	if (M1 < EPS) {
		M1 = 1.0;
	}

	for(i = 0; i < NN1; i++)  {
		Ny1[i] = Ny1[i] / M1;
	}
	
	id = count5;
	for(i = 0; i < NN2; i++)  {
		xx = x2 + i*delta1;

		if(id < count2 - 2 && xx >= datax2[id+1])  {
			id++;
		}
		
		tm1 = datax2[id];
		tm2 = datax2[id+1];
		
		if (fabs(tm2 - tm1) < EPS) {
			Ny2[i] = datay2[id];
		} else {
			Ny2[i] = (tm2 - xx)/(tm2 - tm1)*datay2[id] + (xx - tm1)/(tm2 - tm1)*datay2[id+1];
		}
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
	if (M1 < EPS) {
		M1 = 1.0;
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
	if (cy1*cy2 < EPS) {
		*ratio = 0.0;
	} else {
		*ratio = tt/(cy1*cy2);
	}

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
		printf("Open waveform file failed!\n");
		exit(1);
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
	if (count1 < 2 || count2 < 2) {
		printf("Waveform data points are too few!\n");
		exit(1);
	}

	datax1 = Arraydouble1(count1);
	datay1 = Arraydouble1(count1);
	datax2 = Arraydouble1(count2);
	datay2 = Arraydouble1(count2);

	for(i = 0; i < count1; i++)  {
		//sf=fscanf(fp1, "%lf %lf", &datax1[i], &datay1[i]);
		sf=fscanf(fp1, "%lf %lf %lf", &datax1[i], &tt, &datay1[i]);
		if (sf != 3) {
			printf("Read waveform data failed: %s\n", F1);
			exit(1);
		}
	}
	fclose(fp1);

	for(i = 0; i < count2; i++)  {
		//sf=fscanf(fp2, "%lf %lf", &datax2[i], &datay2[i]);
		sf=fscanf(fp2, "%lf %lf %lf", &datax2[i], &tt, &datay2[i]);
		if (sf != 3) {
			printf("Read waveform data failed: %s\n", F2);
			exit(1);
		}
	}
	fclose(fp2);

	count3 = -1;
	count4 = -1;
	count5 = -1;
	count6 = -1;
	for(i = 0; i < count1 - 1; i++)  {
		if(datax1[i] <= x1 && datax1[i+1] >= x1)  {
			count3 = i;
		}

		if(datax1[i] <= y1 && datax1[i+1] >= y1)  {
			count4 = i;
		}
	}

	for(i = 0; i < count2 - 1; i++)  {
		if(datax2[i] <= x2 && datax2[i+1] >= x2)  {
			count5 = i;
		}

		if(datax2[i] <= y2 && datax2[i+1] >= y2)  {
			count6 = i;
		}
	}
	if (count3 < 0 || count4 < 0 || count5 < 0 || count6 < 0 || count4 <= count3 || count6 <= count5) {
		printf("Window [x,y] is out of waveform range!\n");
		exit(1);
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
	if (delta1 <= EPS) {
		printf("Invalid time interval in waveform data!\n");
		exit(1);
	}

	NN1 = (int)((y1 - x1)/delta1 + 0.0000001) + 1;
	NN2 = (int)((y2 - x2)/delta1 + 0.0000001) + 1;

	Ny1 = Arraydouble1(NN1);
	Ny2 = Arraydouble1(NN2);

	id = count3;
	for(i = 0; i < NN1; i++)  {
		xx = x1 + i*delta1;

		if(id < count1 - 2 && xx >= datax1[id+1])  {
			id++;
		}

		tm1 = datax1[id];
		tm2 = datax1[id+1];

		if (fabs(tm2 - tm1) < EPS) {
			Ny1[i] = datay1[id];
		} else {
			Ny1[i] = (tm2 - xx)/(tm2 - tm1)*datay1[id] + (xx - tm1)/(tm2 - tm1)*datay1[id+1];
		}
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
	if (M1 < EPS) {
		M1 = 1.0;
	}

	for(i = 0; i < NN1; i++)  {
		Ny1[i] = Ny1[i] / M1;
	}
	
	id = count5;
	for(i = 0; i < NN2; i++)  {
		xx = x2 + i*delta1;

		if(id < count2 - 2 && xx >= datax2[id+1])  {
			id++;
		}
		
		tm1 = datax2[id];
		tm2 = datax2[id+1];
		
		if (fabs(tm2 - tm1) < EPS) {
			Ny2[i] = datay2[id];
		} else {
			Ny2[i] = (tm2 - xx)/(tm2 - tm1)*datay2[id] + (xx - tm1)/(tm2 - tm1)*datay2[id+1];
		}
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
	if (M1 < EPS) {
		M1 = 1.0;
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
	if (cy1*cy2 < EPS) {
		*ratio = 0.0;
	} else {
		*ratio = tt/(cy1*cy2);
	}

	free(Ny1);
	free(Ny2);
}

void CrossCorrelation(const char *file_obs, const char *file_syn, double t_start, 
	double t_end, int component, double *time_shift, double *coeff)
{
	char *ff, FF[100];

	if (t_end - t_start < 0.1) {
	fprintf(stderr, "Cross-correlation window too short.\n");
	exit(1);
	}
	FILE *fp1 = fopen(file_obs, "r");
	FILE *fp2 = fopen(file_syn, "r");
	if (!fp1 || !fp2) {
	fprintf(stderr, "Cannot open waveform files.\n");  exit(-1);
	}

	// 读取观测文件
	FF[0] = '\0';
	int count1 = 0;
	while((ff = fgets(FF, 90, fp1)) != NULL)  {
		count1++;
	}
	rewind(fp1);
	int n1 = count1;
	if (n1 < 3) {
	fprintf(stderr, "Error: Insufficient data in %s -- %d\n", file_obs, count1); exit(-1);
	}

	double *t1 = Arraydouble1(n1);
	double *a1 = Arraydouble1(n1);
	double tmp;

	for (int i = 0; i < n1; i++) {
		int nread;
		if (component == 1)
			nread = fscanf(fp1, "%lf %lf %lf", &t1[i], &a1[i], &tmp);
		else
			nread = fscanf(fp1, "%lf %lf %lf", &t1[i], &tmp, &a1[i]);
		if (nread != 3) {
			fprintf(stderr, "Error: Failed to read 3 columns from %s at row %d\n", file_obs, i + 1);
			fclose(fp1);
			exit(1);
		}
	}
	fclose(fp1);

	// 读取合成文件
	FF[0] = '\0';
	int count2 = 0;
	while((ff = fgets(FF, 90, fp2)) != NULL)  {
		count2++;
	}
	rewind(fp2);
	if (count2 < 3) {
	fprintf(stderr, "Error: Insufficient data in %s -- %d\n", file_syn, count2); exit(-1);
	}
	rewind(fp2);
	double *t2 = Arraydouble1(count2);
	double *a2 = Arraydouble1(count2);
	for (int i = 0; i < count2; i++) {
		int nread;
		if (component == 1)
			nread = fscanf(fp2, "%lf %lf %lf", &t2[i], &a2[i], &tmp);
		else
			nread = fscanf(fp2, "%lf %lf %lf", &t2[i], &tmp, &a2[i]);
		if (nread != 3) {
			fprintf(stderr, "Error: Failed to read 3 columns from %s at row %d\n", file_syn, i + 1);
			fclose(fp2);
			free(t1);
			free(a1);
			free(t2);
			free(a2);
			exit(1);
		}
	}
	fclose(fp2);

	// 找到窗口内的起始索引（注意 i < n-1 防止越界）
	int idx1_start = -1, idx1_end = -1;
	for (int i = 0; i < count1-1; i++) {
	if (t1[i] <= t_start && t1[i+1] >= t_start) idx1_start = i;
	if (t1[i] <= t_end   && t1[i+1] >= t_end)   idx1_end = i;
	}
	int idx2_start = -1, idx2_end = -1;
	for (int i = 0; i < count2-1; i++) {
	if (t2[i] <= t_start && t2[i+1] >= t_start) idx2_start = i;
	if (t2[i] <= t_end   && t2[i+1] >= t_end)   idx2_end = i;
	}
	if (idx1_start < 0 || idx1_end < 0 || idx2_start < 0 || idx2_end < 0) {
	fprintf(stderr, "Window out of data range.\n");
	free(t1); free(a1); free(t2); free(a2);
	exit(1);
	}

	// 计算最小采样间隔，用于重采样
	double dt1 = 1e9, dt2 = 1e9;
	for (int i = idx1_start; i < idx1_end; i++) {
	double d = t1[i+1] - t1[i];
	if (d < dt1) dt1 = d;
	}
	for (int i = idx2_start; i < idx2_end; i++) {
	double d = t2[i+1] - t2[i];
	if (d < dt2) dt2 = d;
	}
	double delta = (dt1 < dt2) ? dt1 : dt2;

	// 重采样点数
	int n_new1 = (int)((t_end - t_start) / delta) + 1;
	int n_new2 = (int)((t_end - t_start) / delta) + 1;
	double *y1 = (double*)calloc(n_new1, sizeof(double));
	double *y2 = (double*)calloc(n_new2, sizeof(double));

	// 线性插值观测波形
	int id = idx1_start;
	for (int k = 0; k < n_new1; k++) {
	double tt = t_start + k * delta;
	while (id < idx1_end && tt >= t1[id+1]) id++;
	double w = (tt - t1[id]) / (t1[id+1] - t1[id]);
	y1[k] = (1-w) * a1[id] + w * a1[id+1];
	}
	// 线性插值合成波形
	id = idx2_start;
	for (int k = 0; k < n_new2; k++) {
	double tt = t_start + k * delta;
	while (id < idx2_end && tt >= t2[id+1]) id++;
	double w = (tt - t2[id]) / (t2[id+1] - t2[id]);
	y2[k] = (1-w) * a2[id] + w * a2[id+1];
	}

	free(t1); free(a1); free(t2); free(a2);

	// 两端 Taper 5%
	int ntaper = (int)(0.05 * n_new1);
	if (ntaper < 1) ntaper = 1;
	for (int k = 0; k < ntaper; k++) {
	double w1 = (double)k / ntaper;
	double w2 = (double)(ntaper - k - 1) / ntaper;
	y1[k] *= w1;
	y1[n_new1-1-k] *= w2;
	y2[k] *= w1;
	y2[n_new2-1-k] *= w2;
	}

	// 振幅归一化
	double max1 = 0.0, max2 = 0.0;
	for (int k = 0; k < n_new1; k++) {
	if (fabs(y1[k]) > max1) max1 = fabs(y1[k]);
	}
	for (int k = 0; k < n_new2; k++) {
	if (fabs(y2[k]) > max2) max2 = fabs(y2[k]);
	}
	if (max1 > 0) for (int k = 0; k < n_new1; k++) y1[k] /= max1;
	if (max2 > 0) for (int k = 0; k < n_new2; k++) y2[k] /= max2;

	// 互相关
	double max_corr = -1e9;
	int best_shift = 0;
	for (int shift = -n_new1; shift < n_new2; shift++) {
	double corr = 0.0;
	for (int j = 0; j < n_new1; j++) {
	int idx = j + shift;
	if (idx >= 0 && idx < n_new2) corr += y1[j] * y2[idx];
	}
	if (corr > max_corr) {
	max_corr = corr;
	best_shift = shift;
	}
	}

	// 计算能量归一化系数
	double norm1 = 0.0, norm2 = 0.0;
	for (int k = 0; k < n_new1; k++) norm1 += y1[k]*y1[k];
	for (int k = 0; k < n_new2; k++) norm2 += y2[k]*y2[k];
	norm1 = sqrt(norm1);
	norm2 = sqrt(norm2);
	*coeff = max_corr / (norm1 * norm2);
	*time_shift = t_start - t_start - best_shift * delta;  // 观测减合成

	free(y1); free(y2);
}