#include "Inversion_tt.h"
#include <stdio.h>
#include <unistd.h>

int main(int argc, char *argv[])
{
    int i, j, k, l, m, n, nz, nmx, nmz, num_cor, row, col, count1, count2, step, step_max, flgi1, Ns, Nr;
    double xx1, xx2, x0, x1, zz1, zz2, z0, z1, dx, dz, dis, rr;
    double error, error1, ratio, beta, vt, vt1, g1, gama, tmp1, tmp2, tmp3, maxx, eta, eta1, eta2;
    double *velc, *deltaX, *deltaX1;
    double **arrd_A, *arrd_b, *arrd_x, *data, *Unc, ***arrin;
    Sarray *Fd, *Ft;
    FILE *fp1, *fp2;

    Matrix_infor(&num_cor, &Ns, &Nr, &row, &col, &nz, &step_max, &eta);

    arrd_A = Arraydouble2(row, col);   arrd_b = Arraydouble1(row);    arrd_x = Arraydouble1(col);
    deltaX = Arraydouble1(col);           deltaX1 = Arraydouble1(col);
    data = Arraydouble1(row);             Unc = Arraydouble1(row);        velc = Arraydouble1(col);
    Fd=(Sarray *)malloc(nz*sizeof(Sarray));
    Ft=(Sarray *)malloc(nz*sizeof(Sarray));

    Initial_sarray(nz, Fd);
    Initial_sarray(nz, Ft);
    Initial_velocity0(col, velc);
    Initial_velocity1(col, velc);

    for (i = 0; i < row; i++) {
        Unc[i] = 1.0;
    }
    for (i = 0; i < col; i++) {
        deltaX[i] = 0.0;
    }

    for(step = 0; step < step_max; step++) {
        printf("\n\nInversion step = %d\n", step+1);

        if (step > 0) {
            Get_gradient(num_cor);
            Initial_velocity1(col, velc);
        }

        Update_matrix(&row, &col, &nz);

        Collect_matrix(step, Ns, Nr, row, col, data, arrd_A);

        Sensitivity_matrix(row, col, nz, arrd_A, Fd);

        Transpose(row, col, nz, Fd, Ft);
 
        CG(nz, row, col, eta, Fd, Ft, Unc, data, deltaX);
        // printf("ccc1\n"); sleep(1000);
        Model_update(step+1, col, deltaX, velc);

        Output_delta(step+1, nz, row, col, nz, deltaX);
        // printf("ccc2\n"); sleep(1000);
    }

    free(arrd_b);  free(arrd_x);  free(deltaX);  free(deltaX1);
    free(data);    free(Unc);     free(velc);
    free(Fd);      free(Ft);
    Freedouble2(row, arrd_A);
}

void CG(long nz, long row, long col, double eta, Sarray *Fd, Sarray *Ft, double *Unc, double *data, double *deltaX)
{
    int i, step;
    double error, error1, fixe, ratio, beta, vt, vt1, g1, gama, tmp1, tmp2, tmp3, maxx, eta1, eta2;
    double *Ta, *Pcg, *Gr, *Gr1, *deltaX1;

    Ta = Arraydouble1(row);      Pcg = Arraydouble1(col);
    Gr = Arraydouble1(col);      Gr1 = Arraydouble1(col);
    deltaX1 = Arraydouble1(col);

    for (i = 0; i < col; i++) {
        deltaX[i] = 0.0;
        deltaX1[i] = 0.0;
        Pcg[i] = 0.0;
        Gr[i] = 0.0;
        Gr1[i] = 0.0;
    }
    for (i = 0; i < row; i++) {
        Ta[i] = 0.0;
    }

    Error0(row, data, &error);
    fprintf(stdout, "RMS=%e(s)\n", sqrt(error/row)); fflush(stdout);

    if (error < 1.0e-100) {
        printf("The error is too small!!! \n"); exit(-1);
    }

    eta1 = 1.0e-15; eta2 = 0.0; maxx = 10.0;

    while ((maxx > eta) || (maxx < 0.2*eta && eta1 >= 1.0)) {
        step=0; ratio=1.0;

        for (i = 0; i < col; i++) {
            deltaX[i] = 0.0;
        }

        while (step < 100 && ratio > 0.01) {
            Total_error(nz, row, col, eta1, Fd, data, deltaX, Unc, &error);
            Gradient(row, col, nz, eta1, Ta, Fd, Ft, data, Unc, deltaX, Gr);

            // eta1 = 0.0;

            if (step == 0) {
                fixe = error;

                for (i = 0; i < col; i++) {
                    Pcg[i] = -1.0*Gr[i];
                }
            }
            else {
                tmp1 = 0.0;
                tmp2 = 0.0;
                tmp3 = 0.0;

                for (i = 0; i < col; i++) {
                    tmp1 += Gr[i]*Gr[i];
                    tmp2 += Gr[i]*Gr1[i];
                    tmp3 += Gr1[i]*Gr1[i];
                }
            
                beta = (tmp1 - tmp2) / (tmp3 + 1.0e-100);

                if (beta < 0.0) {
                    beta = 0.0;
                }

                for (i = 0; i < col; i++) {
                    Pcg[i] = -1.0*Gr[i] + beta*Pcg[i];
                }
            }

            g1 = Vector_product(col, Gr, Pcg);

            if(fabs(g1) < 1.0e-100) {
                vt = 0.0;
            }
            else {
                vt = -2.0*error / g1;
            }

            Delta_update(col, vt, Pcg, deltaX, deltaX1);
            Total_error(nz, row, col, eta1, Fd, data, deltaX1, Unc, &error1);

            while (error1 >= error*(1.0 + 1.0e-6)) {
                if (error1 + error < 1.0e-100) {
                    vt1 = 0.0;
                }
                else {
                    vt1 = error / (error1 + error)*vt;
                }

                Delta_update(col, vt1, Pcg, deltaX, deltaX1);
                Total_error(nz, row, col, eta1, Fd, data, deltaX1, Unc, &error1);
                vt = vt1;
            }

            gama = 2.0*((error1 - error) - g1*vt) / (vt*vt + 1.0e-100);

            if (vt*vt < 1.0e-100) {
                vt = 0.0;
            }
            else if (gama >= 1.0e-100) {
                vt = (-1.0*g1)/gama;
            }

            Delta_update(col, vt, Pcg, deltaX, deltaX);
            Total_error(nz, row, col, eta1, Fd, data, deltaX, Unc, &error);

            ratio = error / ((1.0 + 1.0e-3)*fixe);
            step++;

            for (i = 0; i < col; i++) {
                Gr1[i] = Gr[i];
            }
        }
        maxx = Max_abs_vec(col, deltaX);

        if (maxx > eta) {
            eta2 = eta1;
            eta1 = 2.0*eta1;
        }
        else if (maxx < 0.1*eta) {    
            eta1 = 0.5*(eta1 + eta2);
        }
        else {
            break;
        }
    }
    printf("before: %.4e\nafter: %.4e\n", fixe, error);
    fprintf(stdout, "damping=%e; maximal perturbation=%e\n\n", eta1, maxx); fflush(stdout);
    
    free(Ta);   free(Pcg);
    free(Gr);   free(Gr1);
    free(deltaX1);
}

void Model_update(int step, int col, double *deltaX, double *velc)
{
    int nx, ny, nz, i, j, k, count2, offx, offy, offz, expected, fs, base[3];
    double lo[3], hi[3], spacing[3], fraction[3], part[3];
    double **forward, **initial, point[3], vp, vp_per, perturbation, weight;
    char input_inversion_name[128], input_forward_name[128], output_inversion_name[128];
    char output_forward_name[128], output_result_name[128], output_perturb_name[128];
    char ff_command[300];
    FILE *fp_result, *fp_forward, *fp_update, *fp_perturb;

    snprintf(input_inversion_name, sizeof input_inversion_name,
        "../Output/Model/Inver/vel_ini3d.txt");
    snprintf(input_forward_name, sizeof input_forward_name,
            "../Output/Model/Forver/vel3d.txt");
    snprintf(output_inversion_name, sizeof output_inversion_name,
            "../Output/Model/Inver/vel3d.txt");
    snprintf(output_forward_name, sizeof output_forward_name,
            "../Output/Model/Forver/vel3d.txt");
    snprintf(output_result_name, sizeof output_result_name,
        "../Output/Result/vel3d%d.txt", step);
    snprintf(output_perturb_name, sizeof output_perturb_name,
        "../Output/Result/vel_Per%d.txt", step);

    Model_grid_infor(&nx, &ny, &nz, lo, hi);
    expected = nx * ny * nz;
    if (expected != col) {
        fprintf(stderr, "Grid parameters %d differ from matrix columns %d\n", expected, col);
        exit(1);
    }

    spacing[0] = (hi[0]-lo[0])/(nx-1);
    spacing[1] = (hi[1]-lo[1])/(ny-1);
    spacing[2] = (hi[2]-lo[2])/(nz-1);

    initial = Arraydouble2(col, 6);
    Rdfile(input_inversion_name, col, 6, initial);

    fp_result = Open_file(output_result_name, "w");
    fp_update = Open_file(output_inversion_name, "w");
    for (i = 0; i < col; i++) {
        int ix = i/(ny*nz), iy = (i/nz)%ny, iz = i%nz;
        double expected_point[3] = {lo[0]+ix*spacing[0], lo[1]+iy*spacing[1], lo[2]+iz*spacing[2]};
        
        for (j = 0; j < 3; j++) {
            if (fabs(initial[i][j]-expected_point[j]) > 0.051 + 1e-10*fabs(spacing[j])) {
                fprintf(stderr, "%s: grid order mismatch at row %d\n", input_inversion_name, i+1); exit(1);
            }
        }

        velc[i] *= 1.0 - deltaX[i];

        fprintf(fp_result, "%.16e\n", velc[i]);
        fprintf(fp_update, "%.15e %.15e %.15e %.15e %.15e %.15e\n",
                initial[i][0], initial[i][1], initial[i][2],
                velc[i], initial[i][4], initial[i][5]);
    }
    fclose(fp_result);
    fclose(fp_update);

    count2 = Count_file(input_forward_name);
    if (count2 < 1) { fprintf(stderr, "Missing forward velocity model: %s %d\n", input_forward_name, count2); exit(1); }

    forward = Arraydouble2(count2, 6);
    Rdfile(input_forward_name, count2, 6, forward);

    fp_forward = Open_file(output_forward_name, "w");
    fp_perturb = Open_file(output_perturb_name, "w");
    for (i = 0; i < count2; i++) {
        int inside = 1;

        for (j = 0; j < 3; j++) {
            point[j] = forward[i][j];
            if (point[j] < lo[j] || point[j] > hi[j]) inside = 0;
            fraction[j] = (point[j]-lo[j])/spacing[j];
            base[j] = (int)floor(fraction[j]);
        }
        
        vp = forward[i][3];

        if (inside) {
            int size[3] = {nx, ny, nz};
            for (j = 0; j < 3; j++) {
                if (base[j] == size[j]-1) base[j]--;
                part[j] = fraction[j]-base[j];
            }

            vp_per = 0.0;

            for (offx = 0; offx < 2; offx++)
                for (offy = 0; offy < 2; offy++)
                    for (offz = 0; offz < 2; offz++) {
                        k = ((base[0] + offx)*ny + base[1] + offy)*nz + base[2] + offz;

                        weight = (offx ? part[0] : 1.0-part[0]) *
                                 (offy ? part[1] : 1.0-part[1]) *
                                 (offz ? part[2] : 1.0-part[2]);

                        vp_per += weight*deltaX[k];
                    }

            vp = (1.0 - vp_per)*vp;
        }

        if (forward[i][3] == 0.0) {
            fprintf(stderr, "Zero initial Vp at forward model point %d\n", i); exit(1);
        }

        perturbation = (vp - forward[i][3])/forward[i][3];

        fprintf(fp_forward, "%.15e %.15e %.15e %.15e %.15e %.15e\n",
                point[0], point[1], point[2], vp, forward[i][4], forward[i][5]);
        fprintf(fp_perturb, "%.15e %.15e %.15e %.15e %.15e %.15e\n",
                point[0], point[1], point[2], perturbation, 0.0, 0.0);
    }
    fclose(fp_forward);
    fclose(fp_perturb);

    snprintf(output_result_name, sizeof output_result_name,
        "../Output/Result/vel3d_%d.txt", step);
    snprintf(ff_command, sizeof ff_command,
        "cp %s %s", output_forward_name, output_result_name);
    fs = system(ff_command);

    Freedouble2(col, initial);
    Freedouble2(count2, forward);
}

void Delta_update(long col, double vt, double *Pcg, double *deltaX, double *deltaX1)
{
    long i;
    for(i = 0; i < col; i++) {
        deltaX1[i] = deltaX[i] + vt*Pcg[i];
    }    
}

double Vector_product(long len, double *V1, double *V2)
{
    long i;
    double tmp;

    tmp=0.0;
    for(i=0;i<len;i++)
    {
        tmp+=V1[i]*V2[i];
    }

    return tmp;
}

double Max_abs_vec(long len, double *vector)
{
    long i;
    double tmp;

    tmp = 0.0;
    for (i = 0; i < len; i++)
    {
        if (tmp < fabs(vector[i]))
        {
            tmp = fabs(vector[i]);
        }
    }
    return tmp;
}

void Error0(long row, double *data, double *error)
{
    long i;
    double err;

    err = 0.0;
    for(i = 0; i < row; i++) {
        err += data[i]*data[i];
    }
    *error = err / 2.0;
}

void Gradient(long row, long col, long nz, double eta, double *Ta, Sarray *Fd, Sarray *Ft, double *data, double *Unc, double *deltaX, double *Gr)
{
    long i;

    for(i = 0; i < row; i++) {
        Ta[i] = 0.0;
    }
    for(i = 0; i < col; i++) {
        Gr[i] = 0.0;
    }

    for(i = 0; i < nz; i++) {
        Ta[Fd[i].row] += Fd[i].val*deltaX[Fd[i].col];
    }

    for(i = 0; i < row; i++) {
        Ta[i] -= data[i];
        Ta[i] = Ta[i] / Unc[i];
    }

    for(i = 0; i < nz; i++) {
        Gr[Ft[i].row] += Ft[i].val*Ta[Ft[i].col];
    }

    for(i = 0; i < col; i++) {
        Gr[i] = Gr[i] + eta*deltaX[i];
    }
}

void Total_error(long nz, long row, long col, double eta, Sarray *Fd, double *data, double *deltaX, double *Unc, double *error)
{
    long i;
    double err, err1, *tmp;

    tmp = Arraydouble1(row);

    for(i = 0; i < nz; i++) {
        tmp[Fd[i].row] += Fd[i].val*deltaX[Fd[i].col];
    }

    err = 0.0; err1 = 0.0;

    for(i = 0; i < row; i++) {
        err += (tmp[i] - data[i])*(tmp[i] - data[i]) / Unc[i];
    }

    for(i = 0; i < col; i++) {
//        err1 += 0.5*eta*deltaX[i]*deltaX[i];
        err1 += eta*deltaX[i]*deltaX[i];
//        printf("%10.4f %10.4f\n", eta, deltaX[i]);
    }
    *error = (err + err1) / 2.0;
    // *error = err / 2.0;

//    printf("error+++++++++++++++++: %20.4e%20.4e%20.4e\n", err, err1, *error);

    free(tmp);
}

void Sensitivity_matrix(int row, int col, int nz, double **array_a, Sarray *Fd)
{
    int i, j, count1;

    count1 = 0;
    for(i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            Fd[count1].row = i;
            Fd[count1].col = j;
            Fd[count1].val = array_a[i][j];
            count1++;
        }
    }
}

void Transpose(int row, int col, int nz, Sarray *Fd, Sarray *Ft)
{
    int i, trans, *tmp1, *tmp2;

    tmp1=Arrayint1(col);
    tmp2=Arrayint1(col);

    for(i = 0; i < nz; i++) {
        tmp1[Fd[i].col]++;
    }

    for(i = 1; i < col; i++) {
        tmp2[i] = tmp2[i-1] + tmp1[i-1];
    }

    for(i = 0; i < nz; i++) {
        trans = tmp2[Fd[i].col]++;
        Ft[trans].col = Fd[i].row;
        Ft[trans].row = Fd[i].col;
        Ft[trans].val = Fd[i].val;
    }

    free(tmp1); free(tmp2);
}

double interpolate(int n, double x, double y, double *arrfx, double *arrfz, double **arrfy)
{
    int i, j;
    double shapx, shapz, result;
    result = 0.0;

    for (i = 0; i < n - 1; i++)
    {
        if (arrfx[i] == arrfx[i + 1])
            arrfx[i + 1] = arrfx[i] + 0.000001;

        if (arrfz[i] == arrfz[i + 1])
            arrfz[i + 1] = arrfz[i] + 0.000001;
    }
    
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            shapx = lagrange_value(n, i, x, arrfx);
            shapz = lagrange_value(n, j, y, arrfz);

            result += shapx * shapz * arrfy[i][j];
        }
    }

    return result;
}

double interpolate2D(int ix, int jy, double x, double y, double ***arrayData)
{
    int i, j;
    double cof1, cof2, ff;

    ff = 0.0;
    for (i = 0; i < 2; i++) {
        for (j = 0; j < 2; j++) {
            cof1 = (x - arrayData[ix+i][jy+j][0]) / (arrayData[ix+1][jy+1][0] - arrayData[ix][jy][0]);
            cof2 = (y - arrayData[ix+i][jy+j][1]) / (arrayData[ix+1][jy+1][1] - arrayData[ix][jy][1]);
            cof1 = 1.0 - fabs(cof1);
            cof2 = 1.0 - fabs(cof2);

            ff += arrayData[ix+i][jy+j][2]*cof1*cof2;
        }
    }

    return ff;
}

double lagrange_value(int n, int i, double x, double *data_x)
{
    int j;
    double term;

    term = 1.0;

    for (j = 0; j < n; j++)
    {
        if (i != j)
        {
            term *= (x - data_x[j]) / (data_x[i] - data_x[j]);
        }
    }

    return term;
}

int Count_file(const char *filename)
{
    int count = 0;
    int ch;
    int last = '\n';
    FILE *fp = fopen(filename, "r");

    if (fp == NULL) {
        perror(filename);
        return -1;
    }

    while ((ch = fgetc(fp)) != EOF) {
        if (ch == '\n') count++;
        last = ch;
    }

    if (ferror(fp)) {
        perror(filename);
        fclose(fp);
        return -1;
    }

    if (last != '\n') count++;  // 计入没有换行符的最后一行

    fclose(fp);
    return count;
}

FILE *Open_file(const char *flnm, const char *mode)
{
    FILE *fp1;

    if ((fp1 = fopen(flnm, mode)) == NULL)
    {
        perror(flnm); exit(1);
    }
    return fp1;
}

void Rdfile(char *flnm, int row, int col, double **array_data)
{
    int i, j, sf;
    double tt[100];
    char ch, FF[100], *ff;
    FILE *fp;

    fp = Open_file(flnm, "r");
    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            sf = fscanf(fp, "%lf", &tt[j]);
            array_data[i][j] = tt[j];
        }
    }
    fclose(fp);
}

void Initial_sarray(long nz, Sarray *Fd)
{
    long i;
    for(i = 0; i < nz; i++)
    {
        Fd[i].row = 0;
        Fd[i].col = 0;
        Fd[i].val = 0.0;
    }
}

int *Arrayint1(long Nx)
{
    long i;
    int *Array;
    Array=(int *)malloc(Nx*sizeof(int));
    for(i=0;i<Nx;i++)
        Array[i]=0;
    return Array;
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

float *Arrayfloat1(long Nx)
{
    long i;
    float *Array;
    Array=(float *)malloc(Nx*sizeof(float));
    for(i=0;i<Nx;i++)
        Array[i]=0.0;
    return Array;
}

long *Arraylong1(long Nx)
{
    long i;
    long *Array;
    Array = (long *)malloc(Nx*sizeof(long));
    for(i=0;i<Nx;i++)
        Array[i]=0;
    return Array;
}

double **Arraydouble2(long Nx, long Nz)
{
    long i;
    double **Array;

    if (Nx <= 0 || Nz <= 0) {
        fprintf(stderr, "Error in Arraydouble2: Nx=%ld Nz=%ld\n", Nx, Nz);
        exit(EXIT_FAILURE);
    }

    Array = malloc(Nx * sizeof(double *));

    for (i = 0; i < Nx; i++)
        Array[i] = calloc(Nz, sizeof(double));

    return Array;
}

double ***Arraydouble3(long N1,long N2,long N3)
{
    long i,j,k;
    double ***Array;

    Array = (double ***)malloc(N1*sizeof(double **));

    for(i=0;i<N1;i++) {
        Array[i] = (double **)malloc(N2*sizeof(double *));
        for(j=0;j<N2;j++) {
            Array[i][j] = (double *)malloc(N3*sizeof(double));
            for(k=0;k<N3;k++) {
                Array[i][j][k]=0.0;
            }
        }
    }
    return Array;
}

void Freedouble2(long Nx, double **Array)
{
    long i;
    for(i=0;i<Nx;i++)
        free(Array[i]);
    free(Array);
}

void Freedouble3(long N1, long N2, double ***Array)
{
    long i, j;
    for(i = 0; i < N1; i++)
    {
        for (j = 0; j < N2; j++)
        {
            free(Array[i][j]);
        }
        
        free(Array[i]);
    }

    free(Array);
}

void Output_delta(int step, int nz, int row, int col, int Nz, double *deltaX)
{
    int i, j, fs;
    char flnm1[100];
    double xx, yy, zz, ff;
    FILE *fp1;

    sprintf(flnm1,"../Output/Result/vel_per_in%d.txt", step);

    if ((fp1 = fopen(flnm1, "w")) == NULL) {
        printf("Cannot open %s\n", flnm1); exit(1);
    }
    for (i = 0; i < col; i++) {
        fprintf(fp1, "%.16e\n", deltaX[i]);
    }
    fclose(fp1);
}

void Collect_matrix(int step, int Ns, int Nr, int row, int col, double *data, double **A)
{
    int i, j, count1, count2;
    int ii1, ii2, fs;
    double tt1, tt2;
    char  file_sen[100],  file_error[100], cf[250];
    FILE *fp1, *fp2;

    sprintf( file_sen, "../Output/Model/Inver/Sensitivity_in_tt.dat");
    sprintf( file_error, "../Output/Error/Error.txt");

    count1 = Count_file( file_sen);
    count2 = Count_file( file_error);

    if (row <= 0 || col <= 0 || count1 != row*col) {
        fprintf(stderr, "%s: expected %d values; found %d\n",  file_sen, row*col, count1); exit(1);
    }
    if (count2 != row) {
        fprintf(stderr, "%s: expected %d rows; found %d\n",  file_error, row, count2); exit(1);
    }

    fp1 = Open_file( file_error, "r");
    fp2 = Open_file( file_sen, "r");
    for (i = 0; i < row; i++) {
        if (fscanf(fp1, "%d %d %lf %lf", &ii1, &ii2, &tt1, &tt2) != 4) {
            fprintf(stderr, "%s: invalid row %d\n",  file_error, i+1); exit(1);
        }

        data[i] = -1.0*tt1;

        for (j = 0; j < col; j++) {
            if (fscanf(fp2, "%lf", &A[i][j]) != 1) {
                fprintf(stderr, "%s: invalid value at row %d, col %d\n",  file_sen, i+1, j+1); exit(1);
            }
        }
    }
    fclose(fp1);
    fclose(fp2);

    sprintf( file_sen, "../Output/Error/Error_1.txt");
    sprintf( file_error, "../Output/Error/Error%d.txt", step);
    sprintf(cf, "cp %s %s", file_sen, file_error);
    fs = system(cf);
    if (fs != 0) { fprintf(stderr, "Failed to copy Error_1.txt\n"); exit(1); }
}

/*
void Collect_matrix(int step, int Ns, int Nr, long row, long col, double *data, double **A)
{
    int i, j, ii1, ii2, ii3, id1, id2, fs, count1, count2, count3;
    double tt1, tt2, ff1, **arrd1;
    char flnm1[100], flnm2[100], cf[250];
    FILE *fp1;

    sprintf(flnm1, "../Output/Model/Inver/Sensitivity_in.dat");
    count1 = Count_file(flnm1);
    arrd1 = Arraydouble2(count1, 1);
    Rdfile(flnm1, count1, 1, arrd1);

    if (count1 != row*col) {
        printf("The data is not marching! -- Collect_matrix\n"); exit(1);
    }

    sprintf(flnm2, "../Output/Error/Error.txt");
    count2 = Count_file(flnm2);

    if (count2 != row) {
        printf("The data is not marching! -- Collect_matrix\n"); exit(1);
    }
    else {
        sprintf(flnm2, "../Output/Error/Error.txt");

        fp1 = fopen(flnm2, "r");
        if (fp1 == NULL) {
            printf("Can not open tisn file -- %s\n", flnm2); exit(1);
        }
        for (i = 0; i < row; i++) {
            fs = fscanf(fp1, "%d %d %lf %lf", &ii1, &ii2, &tt1, &tt2);
            data[i] = -1.0*tt1;
            // data[i] = -1.0;
        }
        fclose(fp1);
    }

    id2 = 0;
    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            A[i][j] = arrd1[id2][0];
            id2++;
        }
    }

    Freedouble2(count1, arrd1);
    sprintf(flnm1, "../Output/Error/Error_1.txt");
    sprintf(flnm2, "../Output/Error/Error%d.txt", step);
    sprintf(cf, "cp %s %s", flnm1, flnm2);
    fs = system(cf);
}
*/

void Matrix_infor(int *num_cor, int *Ns, int *Nr, int *row, int *col, int *nz, int *step, double *eta)
{
    char rows[18][500];
    int count = Read_config_rows(rows), offset = count - 17;
    int nx, ny, nz3;
    double lo[3], hi[3];
    long long columns, entries;

    Model_grid_infor(&nx, &ny, &nz3, lo, hi);
    if (sscanf(rows[4+offset], "%d", Ns) != 1 ||
        sscanf(rows[5+offset], "%d", Nr) != 1 ||
        sscanf(rows[10+offset], "%d", step) != 1 ||
        sscanf(rows[11+offset], "%lf", eta) != 1 ||
        sscanf(rows[16+offset], "%d", num_cor) != 1 ||
        *Ns < 1 || *Nr < 1 || *step < 1 || *num_cor < 1 || *eta < 0.0) {
        fprintf(stderr, "Invalid source, receiver, inversion, or runtime parameters\n"); exit(1);
    }
    columns = (long long)nx*ny*nz3;
    entries = (long long)*Ns * *Nr * columns;
    if (columns > INT_MAX || entries > INT_MAX) {
        fprintf(stderr, "Inversion matrix exceeds int index limit\n"); exit(1);
    }
    *col = (int)columns;
    *row = *Ns * *Nr;  /* 最大容量；各轮有效行数由 Update_matrix 更新 */
    *nz = (int)entries;
}

void Update_matrix(int *row, int *col, int *nz)
{
    double d1, d2;

    d1 = Count_file("../Output/Error/Error.txt");
    d2 = *col;

    *row = d1;
    *col = d2;
    *nz = d1*d2;
}

void Get_gradient(int num)
{
    int i, fs;
    char ff[512], python_command[128]; ;

    Read_python_command(python_command);

    /// Extracting 2D slices from 3D ///
    sprintf(ff, "cd ./Model/ && %s extract_2d_slices_from_3d.py && cd ../", python_command);
    fs = system(ff);
    if (fs != 0) { fprintf(stderr, "extract_2d_slices_from_3d.py failed\n"); exit(1); }
 
    /// Calculating forward and adjoint wavefields ///
    sprintf(ff, "mpirun -np %d ./Update", num);
    fs = system(ff);
    if (fs != 0) { fprintf(stderr, "Update failed\n"); exit(1); }

   /// Sorting error ///
    snprintf(ff, sizeof(ff), "cd ./Source/ && %s time_diff.py && cd ../", python_command);
    fs = system(ff);
    if (fs != 0) { fprintf(stderr, "time_diff.py failed\n"); exit(1); }

    /// Calculating sensitivity matrix ///
    snprintf(ff, sizeof(ff),
    "cd ./Sens && "
    "%s ./soomth_optimized.py && "
    "%s ./integral2d_to_inversion3d.py && "
    "%s ./kernel_tt.py && "
    "cd ../",
    python_command, python_command, python_command);
    fs = system(ff);
    if (fs != 0) { fprintf(stderr, "3D sensitivity pipeline failed\n"); exit(1); }
}

void Initial_velocity0(int col, double *velc)
{
    int i, fs;
    char flnm1[100];
    double tt1, tt2, tt3, tt4, tt5, tt6;
    FILE *fp1;

    sprintf(flnm1, "../Output/Model/Inver/vel_ini3d.txt");
    
    if ((fp1 = fopen(flnm1, "r")) == NULL) {
        printf("Can not open this file -- %s -- Initial_velocity\n", flnm1); exit(1);
    }

    for (i = 0; i < col; i++) {
        fs = fscanf(fp1, "%lf%lf%lf%lf%lf%lf", &tt1, &tt2, &tt3, &tt4, &tt5, &tt6);
        if (fs != 6) { fprintf(stderr, "%s: invalid row %d\n", flnm1, i+1); exit(1); }
        velc[i] = tt4;
    }
    
    fclose(fp1);
}

void Initial_velocity1(int col, double *velc)
{
    int i, count, fs;
    char flnm1[100];
    double tt1, tt2, tt3, tt4, tt5, tt6;
    FILE *fp1;

    sprintf(flnm1, "../Output/Model/Inver/vel3d.txt");
    if ((fp1 = fopen(flnm1, "r")) == NULL) {
        printf("Can not open this file -- %s -- Initial_velocity1\n", flnm1); exit(1);
    }

    count = Count_file(flnm1);
    if (count != col) {
        printf("The data is not marching! -- Initial_velocity1\n"); exit(1);
    }

    for (i = 0; i < col; i++) {
        fs = fscanf(fp1, "%lf%lf%lf%lf%lf%lf", &tt1, &tt2, &tt3, &tt4, &tt5, &tt6);
        if (fs != 6) { fprintf(stderr, "%s: invalid row %d\n", flnm1, i+1); exit(1); }
        velc[i] = tt4;
    }
    
    fclose(fp1);
}

void modInforIn(int *numx, int *numz, double *x0, double *x1, double *z0, double *z1)
{
    int i, j, ii1, ii2, sf1;
    double dd1, dd2;
    char cc1[500], flnm1[100], flnm2[100], *sf2;
    FILE *fp1;

    sprintf(flnm1, "../Input/parameters_inversion.config");

    if ((fp1 = fopen(flnm1, "r")) == NULL)
    {
        printf("Can not read this file -- %s", flnm1); exit(1);
    }
    ////////////////////////////////////////////////////////////////////////////
    sf1 = fscanf(fp1, "%lf%lf", &dd1, &dd2);  sf2 = fgets(cc1, 500, fp1); // 1
    *x0 = dd1; *x1 = dd2;
    sf1 = fscanf(fp1, "%lf%lf", &dd1, &dd2);  sf2 = fgets(cc1, 500, fp1); // 2
    *z0 = dd1; *z1 = dd2;
    sf1 = fscanf(fp1, "%d%d", &ii1, &ii2);  sf2 = fgets(cc1, 500, fp1);   // 3
    *numx = ii1; *numz = ii2;
}

/* 配置行按 read_file.py 的 18 行格式解析；c: 后的内容是注释。 */
static int Read_config_rows(char rows[18][500])
{
    char line[500], *comment, *start;
    int count = 0;
    FILE *fp = fopen("../Input/parameters_inversion3d.config", "r");
    if (!fp) { perror("parameters_inversion3d.config"); exit(1); }
    while (fgets(line, sizeof(line), fp)) {
        comment = strstr(line, "c:");
        if (comment) *comment = '\0';
        start = line;
        while (*start == ' ' || *start == '\t') start++;
        if (*start == '\0' || *start == '\n') continue;
        if (count == 18) { fprintf(stderr, "Too many config rows\n"); exit(1); }
        strcpy(rows[count++], start);
    }
    fclose(fp);
    if (count != 18) {
        fprintf(stderr, "Expected 18 config rows; got %d\n", count); exit(1);
    }
    return count;
}

static void Model_grid_infor(int *nx, int *ny, int *nz, double lo[3], double hi[3])
{
    char rows[18][500];
    Read_config_rows(rows);
    if (sscanf(rows[0], "%lf %lf", &lo[0], &hi[0]) != 2 ||
        sscanf(rows[1], "%lf %lf", &lo[1], &hi[1]) != 2 ||
        sscanf(rows[2], "%lf %lf", &lo[2], &hi[2]) != 2 ||
        sscanf(rows[3], "%d %d %d", nx, ny, nz) != 3 ||
        *nx < 2 || *ny < 2 || *nz < 2 ||
        hi[0] <= lo[0] || hi[1] <= lo[1] || hi[2] <= lo[2]) {
        fprintf(stderr, "Invalid 3D inversion grid in config\n"); exit(1);
    }
}

static void Read_python_command(char command[128])
{
    char rows[18][500];
    Read_config_rows(rows);
    if (sscanf(rows[16], "%127s", command) != 1) {
        fprintf(stderr, "Missing Python command in parameters_inversion3d.config\n");
        exit(1);
    }
}