#include "Common.h"
#include "Interpolation.h"
#include "Rfile.h"
#include <mpi.h>
#include <time.h>
// int Count_file(char *FL);
// double **Arraydouble2(long Nx,long Nz);
// FILE *Open_file(const char *flnm, const char *mode);
// void Average_rank(int rank, int size, int Ns, int *t1, int *t2);
// void Freedouble2(long Nx, double **Array);
void Matrix_infor(int *Ns, int *Nr, int *numx, int *numz, double *x0, double *x1, double *z0, double *z1, double *dx, double *dz);
// void Read_file(char *flnm, int row, int col, double **array_data);

int main(int argc, char *argv[])
{
    clock_t start, finish;
    int i, j, k, l, n, s1, s2, rank, size, thd, nmx, nmz, count1, count2, count3, count4, flgi, fs;
    int id1, Ns, Nr, numArri, numData, numVel, numModIn, numKerF, numId[4];
    int numX, numY, numNs, numNr, numStep, numCore;
    int numX1, numY1, numKerT;
    double dd, dx, dx1, dy, dy1, deg, rr1, rr2, xx1, xx2, yy1, yy2, ff1, ff2, ff3, ff4;
    double x1, x2, y1, y2, maxPer, dampFct, soomFct, extFct, waveLg;
    double *arraySenSum, **arrayArri, **arrayData, **arrayVel, **arrayVelIn, **arrayKerF, **arrayKerT;
    double **arrayVp, **arrayVs, **arrayPho;
    double fVp, fVs, fPho, fKf, fa, fb, fc, fd, locNode[2];
    char fileArri[100], fileVel[100], fileKerF[100], fileKerSI[100], fileSen[100], fileSenIn[100], fileSenSum[100];
    Vec2 vecP, vecA, vecB, vecC, vecD;
    FILE *fp1, *fp2, *fp3;
    /// Time start ///
    start=clock();
    ///////////////// Load inversion parameters /////////////////
    Read_parInver2d(&numX, &numY, &numNs, &numNr, &numStep, &numCore, &x1, &x2, &y1, &y2, &maxPer, &dampFct, &soomFct, &extFct, &waveLg);
    /////////////////// Load travel time data ///////////////////
    sprintf(fileArri, "../Input/arrival_p.txt");
    numArri = Count_file(fileArri);
    arrayArri = Arraydouble2(numArri, 3);
    Read_file(fileArri, numArri, 3, arrayArri);
    ///
    ///////////////////////// Mpi start! ////////////////////////
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Average_rank(rank, size, numArri, &s1, &s2);
    // Average_rank(rank, size, 1, &s1, &s2);
    ///
    //////////// The array of the model for inversion ///////////
    dx = (x2 - x1)/(numX - 1);
    dy = (y2 - y1)/(numY - 1);

    dx1 = dx/2.0;
    dy1 = dy/2.0;
    numX1 = (x2 - x1)/dx1 + 1.0;
    numY1 = (y2 - y1)/dy1 + 1.0;

    if (numX1 < 100) {
        numX1 = 100;
        dx1 = (x2 - x1)/(numX1 - 1);
    }
    if (numY1 < 100) {
        numY1 = 100;
        dy1 = (y2 - y1)/(numY1 - 1);
    }

    numModIn = numX*numY;
    numKerT = numX1*numY1;

    arrayVelIn = Arraydouble2(numModIn, 3);
    arrayKerT = Arraydouble2(numKerT, 3);

    id1 = 0;
    for (i = 0; i < numX; i++) {
        for (j = 0; j < numY; j++) {
            arrayVelIn[id1][0] = x1 + i*dx;
            arrayVelIn[id1][1] = y1 + j*dy;
            arrayVelIn[id1][2] = 0.0;
            id1++;
        }
    }
    id1 = 0;
    for (i = 0; i < numX1; i++) {
        for (j = 0; j < numY1; j++) {
            arrayKerT[id1][0] = x1 + i*dx1;
            arrayKerT[id1][1] = y1 + j*dy1;
            arrayKerT[id1][2] = 0.0;
            id1++;
        }
    }
    ///
    /////////// Output the soomth kernel for inversion //////////
    for (i = s1; i < s2; i++) {
        Ns = arrayArri[i][0];
        Nr = arrayArri[i][1];

        //////////// Get the vauel of kernel (forward) /////////////
        sprintf(fileVel, "../Output/Model/Forver/vel.dat");
        sprintf(fileKerF, "../Output/Kernel/Kernel_vel_s%d_r%d.dat", Ns, Nr);
        sprintf(fileKerSI, "../Output/Kernel/Smooth_s%d_r%d.dat", Ns, Nr);
        sprintf(fileSen, "../Output/Kernel/Sensitivity_s%d_r%d.dat", Ns, Nr);

        numVel = Count_file(fileVel);
        numKerF = Count_file(fileKerF);

        arrayVel = Arraydouble2(numVel, 5);
        arrayVp = Arraydouble2(numVel, 3);
        arrayVs = Arraydouble2(numVel, 3);
        arrayPho = Arraydouble2(numVel, 3);
        arrayKerF = Arraydouble2(numKerF, 3);

        Read_file(fileVel, numVel, 5, arrayVel);
        Read_file(fileKerF, numKerF, 3, arrayKerF);

        for (j = 0; j < numVel; j++) {
            arrayVp[j][0] = arrayVel[j][0];
            arrayVp[j][1] = arrayVel[j][1];
            arrayVp[j][2] = arrayVel[j][2];

            arrayVs[j][0] = arrayVel[j][0];
            arrayVs[j][1] = arrayVel[j][1];
            arrayVs[j][2] = arrayVel[j][3];

            arrayPho[j][0] = arrayVel[j][0];
            arrayPho[j][1] = arrayVel[j][1];
            arrayPho[j][2] = arrayVel[j][4];
        }

        for (j = 0; j < numKerF; j++) {
            arrayKerF[j][1] = -1.0*arrayKerF[j][1];
        }
        //////
        for (j = 0; j < numKerT; j++) {
            locNode[0] = arrayKerT[j][0];
            locNode[1] = arrayKerT[j][1];
            /// Interpolation ///
            flgi = Fine_node_square(numVel, numId, arrayVp, locNode);
            /// Vp ///
            if (flgi != -1) {
                vecP.x = locNode[0];            vecP.y = locNode[1];
                vecA.x = arrayVp[numId[0]][0];  vecA.y = arrayVp[numId[0]][1];  fa = arrayVp[numId[0]][2];
                vecB.x = arrayVp[numId[1]][0];  vecB.y = arrayVp[numId[1]][1];  fb = arrayVp[numId[1]][2];
                vecC.x = arrayVp[numId[2]][0];  vecC.y = arrayVp[numId[2]][1];  fc = arrayVp[numId[2]][2];
                vecD.x = arrayVp[numId[3]][0];  vecD.y = arrayVp[numId[3]][1];  fd = arrayVp[numId[3]][2];
                fVp = Interplot_2din(vecP, vecA, vecB, vecC, vecD, fa, fb, fc, fd);
            }
            if (flgi == -1 || isnan(fVp) == 1) {
                id1 = Nearest_point2d(locNode, arrayVp, numVel);
                fVp = arrayVp[id1][2];
            }
            /// Vs ///
            if (flgi != -1) {
                vecP.x = locNode[0];             vecP.y = locNode[1];
                vecA.x = arrayVs[numId[0]][0];  vecA.y = arrayVs[numId[0]][1];  fa = arrayVs[numId[0]][2];
                vecB.x = arrayVs[numId[1]][0];  vecB.y = arrayVs[numId[1]][1];  fb = arrayVs[numId[1]][2];
                vecC.x = arrayVs[numId[2]][0];  vecC.y = arrayVs[numId[2]][1];  fc = arrayVs[numId[2]][2];
                vecD.x = arrayVs[numId[3]][0];  vecD.y = arrayVs[numId[3]][1];  fd = arrayVs[numId[3]][2];
                fVs = Interplot_2din(vecP, vecA, vecB, vecC, vecD, fa, fb, fc, fd);
            }
            if (flgi == -1 || isnan(fVs) == 1) {
                id1 = Nearest_point2d(locNode, arrayVs, numVel);
                fVs = arrayVs[id1][2];
            }
            /// Pho ///
            if (flgi != -1) {
                vecP.x = locNode[0];             vecP.y = locNode[1];
                vecA.x = arrayPho[numId[0]][0];  vecA.y = arrayPho[numId[0]][1];  fa = arrayPho[numId[0]][2];
                vecB.x = arrayPho[numId[1]][0];  vecB.y = arrayPho[numId[1]][1];  fb = arrayPho[numId[1]][2];
                vecC.x = arrayPho[numId[2]][0];  vecC.y = arrayPho[numId[2]][1];  fc = arrayPho[numId[2]][2];
                vecD.x = arrayPho[numId[3]][0];  vecD.y = arrayPho[numId[3]][1];  fd = arrayPho[numId[3]][2];
                fPho = Interplot_2din(vecP, vecA, vecB, vecC, vecD, fa, fb, fc, fd);
            }
            if (flgi == -1 || isnan(fPho) == 1) {
                id1 = Nearest_point2d(locNode, arrayPho, numVel);
                fPho = arrayPho[id1][2];
            }
            /// Uxx + Uzz ///
            flgi = Fine_node_square(numKerF, numId, arrayKerF, locNode);
            if (flgi != -1) {
                vecP.x = locNode[0];              vecP.y = locNode[1];
                vecA.x = arrayKerF[numId[0]][0];  vecA.y = arrayKerF[numId[0]][1];  fa = arrayKerF[numId[0]][2];
                vecB.x = arrayKerF[numId[1]][0];  vecB.y = arrayKerF[numId[1]][1];  fb = arrayKerF[numId[1]][2];
                vecC.x = arrayKerF[numId[2]][0];  vecC.y = arrayKerF[numId[2]][1];  fc = arrayKerF[numId[2]][2];
                vecD.x = arrayKerF[numId[3]][0];  vecD.y = arrayKerF[numId[3]][1];  fd = arrayKerF[numId[3]][2];
                fKf = Interplot_2din(vecP, vecA, vecB, vecC, vecD, fa, fb, fc, fd);
            }
            if (flgi == -1 || isnan(fKf) == 1) {
                id1 = Nearest_point2d(locNode, arrayKerF, numVel);
                fKf = arrayKerF[id1][2];
            }
            arrayKerT[j][2] = 2.0*fVp*fVp*fPho*fKf*dx*dy;
            //////
        }

        /////////// Get the soomth kernel (inversion)  ///////////
        fp1 = Open_file(fileSen, "w");
        // fp1 = Open_file("./Sens/vel.dat", "w");
        for (j = 0; j < numModIn; j++) {
            dd = waveLg/1.0;
            // dd = dx + dy;
            thd = 20.0;
            xx1 = arrayVelIn[j][0];
            yy1 = arrayVelIn[j][1];
            fa = 4.0/(pi*dd*dd);

            ff1 = 0.0;
            ff2 = 0.0;
            for (k = 0; k < numKerT; k++) {
                xx2 = arrayKerT[k][0];
                yy2 = arrayKerT[k][1];

                rr1 = fabs(xx2 - xx1) / dd;
                rr2 = fabs(yy2 - yy1) / dd;

                if (rr1 > thd && rr2 > thd) {
                    continue;
                }

                fb = -4.0*((xx2 - xx1)*(xx2 - xx1) + (yy2 - yy1)*(yy2 - yy1))/(dd*dd);
                fc = fa*exp(fb);

                ff1 += arrayKerT[k][2]*fc;
                ff2 += fc;
            }

            arrayVelIn[j][2] = (ff1/ff2);

            fprintf(fp1, "%.4f %10.4f %20.8e\n", arrayVelIn[j][0], -1.0*arrayVelIn[j][1], arrayVelIn[j][2]);
        }
        fclose(fp1);
        /////////// Output the progress of the program  ////////////
        if (rank == 0) {
            deg = 100.0*(i+1)/s2;
            printf("%10.1f%%\n", deg);
        }
        //////
        Freedouble2(numVel, arrayVel);
        Freedouble2(numVel, arrayVp);
        Freedouble2(numVel, arrayVs);
        Freedouble2(numVel, arrayPho);
        Freedouble2(numKerF, arrayKerF);
    }
    //////////////////////// Free arrays  //////////////////////////
    Freedouble2(numArri, arrayArri);
    Freedouble2(numModIn, arrayVelIn);
    Freedouble2(numKerT, arrayKerT);
    //////
    ///////////////////////// Data collect /////////////////////////
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        int ns, nr1, nr2, numDD;
        double Sens1[3], Sens2[3], **arrayDD;
        char fileTimeDD[100], fileSenDD[100], fileSens1[100], fileSens2[100];
        FILE *fp4, *fp5;

		sprintf(fileTimeDD, "../Output/Error/Error_dd.txt");
		sprintf(fileSenDD, "../Output/Model/Inver/Sensitivity_in_dd.dat");
        sprintf(fileSenIn, "../Output/Model/Inver/Sensitivity_in.dat");
        sprintf(fileSenSum, "../Output/Model/Inver/Sensitivity.dat");
        sprintf(fileArri, "../Input/arrival_p.txt");

        numArri = Count_file(fileArri);
        numDD = Count_file(fileTimeDD);

        arraySenSum = Arraydouble1(numModIn);
        arrayVelIn = Arraydouble2(numModIn, 3);
        arrayDD = Arraydouble2(numDD, 4);

        Read_file(fileTimeDD, numDD, 4, arrayDD);
        Read_file("../Output/Kernel/Sensitivity_s1_r1.dat", numModIn, 3, arrayVelIn);

        fp1 = Open_file(fileSenDD, "w");
        fp2 = Open_file(fileSenSum, "w");
        for (i = 0; i < numDD; i++) {
			ns = arrayDD[i][0];
			nr1 = arrayDD[i][1];
			nr2 = arrayDD[i][2];

			sprintf(fileSens1, "../Output/Kernel/Sensitivity_s%d_r%d.dat", ns, nr1);
			sprintf(fileSens2, "../Output/Kernel/Sensitivity_s%d_r%d.dat", ns, nr2);

			fp4 = Open_file(fileSens1, "r");
			fp5 = Open_file(fileSens2, "r");
			for (j = 0; j < numModIn; j++)  {
				fs = fscanf(fp4, "%lf %lf %lf", &Sens1[0], &Sens1[1], &Sens1[2]);
				fs = fscanf(fp5, "%lf %lf %lf", &Sens2[0], &Sens2[1], &Sens2[2]);

                arrayVelIn[j][2] = Sens1[2] - Sens2[2];
                arraySenSum[j] += arrayVelIn[j][2];

                fprintf(fp1, "%.8e\n", arrayVelIn[j][2]);
			}
			fclose(fp4);
			fclose(fp5);
        }

        for (j = 0; j < numModIn; j++) {
            fprintf(fp2, "%.4f %10.4f %.16e\n", arrayVelIn[j][0], -1.0*arrayVelIn[j][1], arraySenSum[j]);
        }

        fclose(fp1);
        fclose(fp2);
        free(arraySenSum);
        Freedouble2(numModIn, arrayVelIn);
        Freedouble2(numDD, arrayDD);
    }
    /////////// Mpi End  ///////////
    MPI_Finalize();
    return 0;
    /////////
}
