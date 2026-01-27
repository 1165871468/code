#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "SEM.h"

double Read_par2d(int Nr, int Ns, double *X, double *Z, double *Rx, double *Rz, double *Sx, double *Sz)
{
	int i, id, count, Nr1, Ns1;
	double Rx1, Rz1, Sx1, Sz1, X1, Z1, **modelInfors;
	char flnm[100];

	sprintf(flnm, "../Output/Mesh/Model_infor.txt");
	count = count_file(flnm);
	modelInfors = Arraydouble2(count, 19);
	rdfile(flnm, count, 19, modelInfors);

	id = -1;
	for (i = 0; i < count; i++)
	{
		Ns1 = modelInfors[i][0];    Nr1 = modelInfors[i][1];
		Sx1 = modelInfors[i][8];    Sz1 = modelInfors[i][9];
		Rx1 = modelInfors[i][10];   Rz1 = modelInfors[i][11];
		X1 = modelInfors[i][17];    Z1 = modelInfors[i][18];
		if (Ns == Ns1 && Nr == Nr1)
		{
			id = 1;  break;
		}
	}
	*X = X1;    *Z = Z1; 
	*Rx = Rx1;  *Rz = -Rz1;
	*Sx = Sx1;  *Sz = -Sz1;

	Freedouble2(count, modelInfors);
	return id;
}

void VS_StrTrim(char*pStr)
{  
    char *pTmp = pStr;  
     
    while (*pStr != '\0')   
    {  
        if (*pStr != ' ')  
        {  
            *pTmp++ = *pStr;
        }  
        ++pStr;  
    }  
    *pTmp = '\0';  
}

void rdfile(char *flnm, int row, int col, double **arrays)
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
            arrays[i][j] = tt[j];
        }
    }
    fclose(fp);
}

void Read_par(short *N,short *Pflag,short *ndiv,short *nstep,short *nsnap,short *sflag,short *stype,long *Nr_s,long *Nr,double *X,double *Z,double *T,double *dt,double *F,
	double *Vp,double *Layx,double *Layz,double *R0,double *Lscale,char *source,char *receiver,char *mesh,char *moment)
{
	short ss;
	long nn;
	int i1,i2,i3,ff1;
	double t1,t2,t3;
	char ct[500],file[200],*ff;
	FILE *fp;

	file[0]='\0';
	strcat(file,"../Input/parameters");				//打开parameters

	if((fp=fopen(file,"r"))==NULL)
	{
		printf("Cannot open the files parameters\n"); exit(-1);
	}
	
	ff1=fscanf(fp,"%hd",&ss);
	*N=ss;
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%lf %lf",&t1,&t2);
	*X=t1; *Z=t2;
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%lf",&t1);
	*T=t1;
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%lf",&t1);
	*dt=t1;
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%lf",&t1);
	*F=t1;
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%hd",&ss);
	*nstep=ss;
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%hd",&ss);
	*nsnap=ss;
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%hd",&ss);
	*sflag=ss;
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%lf",&t1);
	*Vp=t1;
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%lf",&t1);
	*Layx=t1;
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%lf",&t1);
	*Layz=t1;
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%lf",&t1);
	*R0=t1;
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%hd",&ss);
	*Pflag=ss;
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%lf",&t1);
	*Lscale=t1;
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%ld",&nn);
	*Nr_s=nn;
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%ld",&nn);
	*Nr=nn;
	ff=fgets(ct,500,fp);

	ff=fgets(source,30,fp);
	ff=fgets(ct,500,fp);
	VS_StrTrim(source);

	ff=fgets(receiver,30,fp);
	ff=fgets(ct,500,fp);
	VS_StrTrim(receiver);
	
	ff=fgets(mesh,30,fp);
	ff=fgets(ct,500,fp);
	VS_StrTrim(mesh);

	ff1=fscanf(fp,"%hd",&ss);
	*ndiv=ss;
	ff=fgets(ct,500,fp);

	ff=fgets(moment,30,fp);
	ff=fgets(ct,500,fp);
	VS_StrTrim(moment);

	ff1=fscanf(fp,"%hd",&ss);
	*stype=ss;
	ff=fgets(ct,500,fp);
		
	fclose(fp);
}

void Read_parForwar(short *N, short *Pflag, short *ndiv, short *nstep, short *nsnap, short *sflag, short *stype, double *dt,
	double *F, double *Vp, double *Layx, double *Layz, double *R0, double *Lscale, double *X, double *Z, char *moment)
{
    short N1, Pflag1, ndiv1, nsnap1, nstep1, sflag1, stype1, sf1;
	char cc1[500], file[200], *sf2;
	FILE *fp1;

	sprintf(file, "../Input/parameters_forward.config");
	fp1 = Open_file(file, "r");
    ////////////////////////////////////////////////////////////////////////////
    sf1 = fscanf(fp1, "%hd", N);                      sf2 = fgets(cc1, 500, fp1);   // 1
    sf1 = fscanf(fp1, "%lf", dt);                     sf2 = fgets(cc1, 500, fp1);   // 2
    sf1 = fscanf(fp1, "%lf", F);                      sf2 = fgets(cc1, 500, fp1);   // 3
    sf1 = fscanf(fp1, "%hd", nstep);                  sf2 = fgets(cc1, 500, fp1);   // 4
    sf1 = fscanf(fp1, "%hd", nsnap);                  sf2 = fgets(cc1, 500, fp1);   // 5
    sf1 = fscanf(fp1, "%hd", sflag);                  sf2 = fgets(cc1, 500, fp1);   // 6
    sf1 = fscanf(fp1, "%lf", Vp);                     sf2 = fgets(cc1, 500, fp1);   // 7
    sf1 = fscanf(fp1, "%lf", Layx);                   sf2 = fgets(cc1, 500, fp1);   // 8
    sf1 = fscanf(fp1, "%lf", Layz);                   sf2 = fgets(cc1, 500, fp1);   // 9
    sf1 = fscanf(fp1, "%lf", R0);                     sf2 = fgets(cc1, 500, fp1);   // 10
    sf1 = fscanf(fp1, "%hd", Pflag);                  sf2 = fgets(cc1, 500, fp1);   // 11
    sf1 = fscanf(fp1, "%lf", Lscale);                 sf2 = fgets(cc1, 500, fp1);   // 12
    sf1 = fscanf(fp1, "%hd", ndiv);                   sf2 = fgets(cc1, 500, fp1);   // 13
    sf1 = fscanf(fp1, "%s", moment);                  sf2 = fgets(cc1, 500, fp1);   // 14
    sf1 = fscanf(fp1, "%hd", stype);                  sf2 = fgets(cc1, 500, fp1);   // 15
	sf1 = fscanf(fp1, "%lf%lf", X, Z);                sf2 = fgets(cc1, 500, fp1);   // 16
    ////////////////////////////////////////////////////////////////////////////
	fclose(fp1);
}

void Read_location1(double *Sx,double *Sz,char *source)
{
	int ff1;
	char file[200];
	double ftmp1,ftmp2;
	FILE *fp;

	file[0]='\0';
	strcat(file,"../Input/");
	strcat(file,source);
	
	if((fp=fopen(file,"r"))==NULL)
	{
		fprintf(stdout,"Cannot open the file Location\n");
		fflush(stdout);
		exit(0);
	}
	ff1=fscanf(fp,"%lf%lf",&ftmp1,&ftmp2);
	*Sx=ftmp1;
	*Sz=ftmp2;
	
	fclose(fp);
}

void Read_location2(long Nr,double **Rloc,char *source)
{
	int i,ff1;
	char file[200];
	double ftmp1,ftmp2;
	FILE *fp;

	file[0]='\0';
	strcat(file,"../Input/");
	strcat(file,source);
	
	if((fp=fopen(file,"r"))==NULL)
	{
		fprintf(stdout,"Cannot open the file Location\n");
		fflush(stdout);
		exit(0);
	}

	for(i=0;i<Nr;i++)
	{
		ff1=fscanf(fp,"%lf%lf",&ftmp1,&ftmp2);
		Rloc[i][0]=ftmp1;
		Rloc[i][1]=ftmp2;
	}
	
	fclose(fp);
}

void Read_tensor(double **Tensor)
{
	short ss;
	int i1,i2,i3,ff1;
	double t1,t2,t3;
	char ct[100],file[200],*ff;
	FILE *fp;

	file[0]='\0';
	strcat(file,"../Input/");
	strcat(file,"Moment-tensor");

	if((fp=fopen(file,"r"))==NULL)
	{
		fprintf(stdout,"Cannot open the files Moment-tensor\n");
		fflush(stdout);
		exit(0);
	}
	ff=fgets(ct,100,fp);
	ff=fgets(ct,100,fp);
	ff1=fscanf(fp,"%lf%lf",&Tensor[0][0],&Tensor[0][1]);
	ff1=fscanf(fp,"%lf%lf",&Tensor[1][0],&Tensor[1][1]);
	fclose(fp);
}

void loadRecAndSor(int Nr, int Ns, double *Rx, double *Rz, double *Sx, double *Sz)
{
	int i, id, numRec, numSor, Nr1, Ns1;
	double Rx1, Rz1, Sx1, Sz1, X1, Z1, **locRec, **locSor;
	char fileRec[100], fileSor[100];

	sprintf(fileRec, "../Input/Receiver.txt");
	sprintf(fileSor, "../Input/Source.txt");

	numRec = count_file(fileRec);
	numSor = count_file(fileSor);

	locRec = Arraydouble2(numRec, 2);
	locSor = Arraydouble2(numSor, 2);

	rdfile(fileRec, numRec, 2, locRec);
	rdfile(fileSor, numSor, 2, locSor);

	for (i = 0; i < Nr; i++) {
		*Rx = locRec[i][0];
		*Rz = locRec[i][1];
	}

	for (i = 0; i < Ns; i++) {
		*Sx = locSor[i][0];
		*Sz = locSor[i][1];
	}

	Freedouble2(numRec, locRec);
	Freedouble2(numSor, locSor);
}

