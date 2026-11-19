#include "Common.h"

int *Arrayint1(long nx)
{
    long i;
    int *array;

    array = (int *)malloc(nx * sizeof(int));
    for ( i = 0; i < nx; i++)
    {
        array[i] = 0;
    }
    return array;
}

int **Arrayint2(long nx, long ny)
{
    long i,j;
    int **array;

    array = (int **)malloc(nx * sizeof(int *));
    for ( i = 0; i < nx; i++)
    {
        array[i] = (int *)malloc(ny * sizeof(int));
    }

    for ( i = 0; i < nx; i++)
    {
        for ( j = 0; j < ny; j++)
        {
            array[i][j] = 0;
        }
        
    }
    return array;
}

int Count_file(const char *flnm)
{
    int count;
    char ch;
    FILE *fp;

    fp = Open_file(flnm, "r");
    count = 0;
    while ((ch = fgetc(fp)) != EOF)
    {
        if(ch == '\n')
        {
            count++;
        }
    }
    fclose(fp);
    return count;
}

int File_exist(const char* flnm)
{
    if (access(flnm, F_OK) != -1)
    {
        return 1;
    }
    return 0;
}

double *Arraydouble1(long nx)
{
    long i;
    double *array;

    array = (double *)malloc(nx * sizeof(double));
    for ( i = 0; i < nx; i++)
    {
        array[i] = 0.0;
    }
    return array;
}

double **Arraydouble2(long nx, long ny)
{
    long i,j;
    double **array;

    array = (double **)malloc(nx * sizeof(double *));
    for ( i = 0; i < nx; i++)
    {
        *(array + i) = (double *)malloc(ny * sizeof(double));
    }

    for ( i = 0; i < nx; i++)
    {
        for ( j = 0; j < ny; j++)
        {
            array[i][j] = 0.0;
        }
        
    }
    return array;
}


double ***Arraydouble3(long N1,long N2,long N3)
{
	long i,j,k;
	double ***Array;

	Array=(double ***)malloc(N1*sizeof(double **));

	for(i=0;i<N1;i++)
	{
		Array[i]=(double **)malloc(N2*sizeof(double *));
		for(j=0;j<N2;j++)
		{
			Array[i][j]=(double *)malloc(N3*sizeof(double *));
			for(k=0;k<N3;k++)
			{
				Array[i][j][k]=0.0;
			}
		}
	}

	return Array;
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

void Average_rank(int rank, int size, int Ns, int *t1, int *t2)
{
	int i, es, rs, ss, s1, s2;

	es = Ns / size;
	rs = Ns - es*size;
	ss = 0; s1 = 0;
	for(i = 0; i <= rank; i++)
	{
		ss += es;
		if (i < rs)
		{
			ss++;
		}
		if (i == rank - 1)
		{
			s1 = ss;
		}
	}

	s2 = ss;
	*t1 = s1;
	*t2 = s2;
}

void Freeint2(long nx, int **Array)
{
    long i;

    for(i = 0; i < nx; i++)
    {
        free(Array[i]);
    }
    free(Array);
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

void Get_model_infor(FILE *fp1, int num_s, int num_r, double loc_s[], double loc_r[], double loc_new[],
    double rot_theta[], double y_o, double model_dis, double model_dep)
{
    fprintf(fp1, "%d %10d ", num_s, num_r);
    fprintf(fp1, "%15.4f %10.4f %10.4f %.4f %10.4f %10.4f ", earth_r-loc_s[0], loc_s[1], loc_s[2], earth_r-loc_r[0], loc_r[1], loc_r[2]);
    fprintf(fp1, "%10.4f %10.4f %10.4f %10.4f ", loc_new[1], loc_new[0]-earth_r, loc_new[4], loc_new[3]-earth_r);
    fprintf(fp1, "%10.8f %10.8f %10.8f %10.8f ", rot_theta[0], rot_theta[1], rot_theta[2], rot_theta[3]);
    fprintf(fp1, "%10.4f %10.1f %10.1f\n", y_o, model_dis, model_dep);
}

void Model_parameter(int *Ns, int *Nr, double *elong_factor)
{
    int i, j, ii1, ii2, ii3, num_core, numr, nums, numx, numy, numz, num_stp, sf1;
    double extd, dd1, dd2, deltax, deltay, damping, smooth, per_max, xx1, xx2, yy1, yy2, zz1, zz2;
    char cc1[500], flnm1[100], flnm2[100], *sf2;
    FILE *fp1;

    sprintf(flnm1, "../Input/parameters_inversion.config");
    fp1 = Open_file(flnm1, "r");
    ////////////////////////////////////////////////////////////////////////////
    sf1 = fscanf(fp1, "%lf%lf", &xx1, &xx2);           sf2 = fgets(cc1, 500, fp1);   // 1
    sf1 = fscanf(fp1, "%lf%lf", &yy1, &yy2);           sf2 = fgets(cc1, 500, fp1);   // 2
    sf1 = fscanf(fp1, "%lf%lf", &zz1, &zz2);           sf2 = fgets(cc1, 500, fp1);   // 2
    sf1 = fscanf(fp1, "%d%d%d", &numx, &numy, &numz);  sf2 = fgets(cc1, 500, fp1);   // 3
    sf1 = fscanf(fp1, "%d", &nums);                    sf2 = fgets(cc1, 500, fp1);   // 4
    sf1 = fscanf(fp1, "%d", &numr);                    sf2 = fgets(cc1, 500, fp1);   // 5
    sf1 = fscanf(fp1, "%d", &num_stp);                 sf2 = fgets(cc1, 500, fp1);   // 6
    sf1 = fscanf(fp1, "%lf", &per_max);                sf2 = fgets(cc1, 500, fp1);   // 7
    sf1 = fscanf(fp1, "%lf", &damping);                sf2 = fgets(cc1, 500, fp1);   // 8
    sf1 = fscanf(fp1, "%lf", &smooth);                 sf2 = fgets(cc1, 500, fp1);   // 9
    sf1 = fscanf(fp1, "%d", &num_core);                sf2 = fgets(cc1, 500, fp1);   // 10
    sf1 = fscanf(fp1, "%lf", &extd);                   sf2 = fgets(cc1, 500, fp1);   // 11
    ////////////////////////////////////////////////////////////////////////////
    *Nr = numr;
    *Ns = nums;
    *elong_factor = extd;
    fclose(fp1);
}

void Read_file(char *flnm, int row, int col, double **array_data)
{
    int i, j, sf;
    double tt[100];
    char ch, FF[100], *ff;
    FILE *fp;

    fp = Open_file(flnm, "r");
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < col; j++)
        {
            sf = fscanf(fp, "%lf", &tt[j]);
            array_data[i][j] = tt[j];
        }
    }
    fclose(fp);
}
