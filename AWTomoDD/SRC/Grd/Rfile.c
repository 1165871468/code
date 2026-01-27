#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "SEM.h"

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
		fprintf(stdout,"Cannot open the files parameters\n");
		fflush(stdout);
		exit(0);
	}								//读文件错误
	
	ff1=fscanf(fp,"%hd",&ss);					//将文件的各个量写入数组
	*N=ss;								//插值阶数
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%lf %lf",&t1,&t2);
	*X=t1; *Z=t2;							//计算域的水平和垂直比例
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%lf",&t1);
	*T=t1;								//模拟的时间
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%lf",&t1);
	*dt=t1;								//时间间隔
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%lf",&t1);
	*F=t1;								//瑞雷波的主频率
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%hd",&ss);
	*nstep=ss;							//每Nt*dt将时间显示在屏幕上
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%hd",&ss);
	*nsnap=ss;							//每Nt*dt将快照输出
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%hd",&ss);
	*sflag=ss;							//1=计算域中的波形文件 0=计算域加PML层中的波形文件
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%lf",&t1);
	*Vp=t1;								//PML层的Vp
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%lf",&t1);
	*Layx=t1;							//PML层的厚度（x）
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%lf",&t1);
	*Layz=t1;							//PML层的厚度（z）
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%lf",&t1);
	*R0=t1;								//阻尼系数
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%hd",&ss);
	*Pflag=ss;							//PML层的标志
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%lf",&t1);
	*Lscale=t1;							//元素大小的平均尺度
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%ld",&nn);
	*Nr_s=nn;							//震源的数量
	ff=fgets(ct,500,fp);

	ff1=fscanf(fp,"%ld",&nn);
	*Nr=nn;								//接收器的数量
	ff=fgets(ct,500,fp);

	ff=fgets(source,30,fp);						//源文件的文件名字
	ff=fgets(ct,500,fp);
	VS_StrTrim(source);

	ff=fgets(receiver,30,fp);
	ff=fgets(ct,500,fp);						//接收器位置信息的文件名字
	VS_StrTrim(receiver);
	
	ff=fgets(mesh,30,fp);
	ff=fgets(ct,500,fp);						//模型网格信息的文件名字
	VS_StrTrim(mesh);

	ff1=fscanf(fp,"%hd",&ss);
	*ndiv=ss;							//将元素分解为多少小元素，得到源和接收器的位置信息
	ff=fgets(ct,500,fp);

	ff=fgets(moment,30,fp);
	ff=fgets(ct,500,fp);						//矩张量的文件名字
	VS_StrTrim(moment);

	ff1=fscanf(fp,"%hd",&ss);
	*stype=ss;							//1=点源 2=爆炸源
	ff=fgets(ct,500,fp);
		
	fclose(fp);
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

