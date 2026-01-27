#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include "Readmesh.h"

////////////////////////////////read nodes and triangles///////////////////读取节点数量
void Mesh(long *node,long *elem,char *mesh)
{
	char *ff1, A[100]={"kkkkkk"},FF[100];
	FILE *fp;				//linux命令显示文件类型
	long i=0,j=0;
	int ff2, flag=0;
	
	// sprintf(FF, "../Output/Mesh/%s", mesh);

	if((fp=fopen(mesh,  "r")) == NULL) {
		printf("cannot opent Mesh file %s\n", mesh); exit(0);
	}

	ff1 = fgets(A,100,fp);				//读入一行
	ff2 = fscanf(fp,"%ld%ld\n",&i,&j);			//指定参数的分割格式 用，分开
	fclose(fp);
	*node=j;
	*elem=i;
}
///////////////////////////get parameter of the mesh///////////////////////////////////////////////
////////////////////////////////获取网格参数////////////////////////////////////////////	
void Get_par(long node,long Elem,double **node_array,long **Elem_array,short *Elem_par,char *mesh)
{
	char *ff1, A[60]={"kkkkkk"},FF[100];
	FILE *fp;
	int ff2, itmp6;
	long itmp1,itmp2,itmp3,itmp4,itmp5,i;
	double ftmp1,ftmp2;

	sprintf(FF, "%s", mesh);

	if((fp = fopen(FF,"r")) == NULL)
	{
		printf("Cannot opent this file: %s\n", FF); exit(0);
	}
	ff1 = fgets(A,59,fp);
	ff1 = fgets(A,59,fp);
	ff1 = fgets(A,59,fp);
	for(i=0;i<node;i++)
	{
		ff2 = fscanf(fp,"%ld %lf %lf",&itmp1,&ftmp1,&ftmp2);
		node_array[i][0]=(float)ftmp1;
		node_array[i][1]=(float)ftmp2;
	}
	ff1 = fgets(A,59,fp);
	ff1 = fgets(A,59,fp);
	ff1 = fgets(A,59,fp);
	for(i=0;i<Elem;i++)
	{
		ff2 = fscanf(fp,"%ld %ld %ld %ld %ld %d",&itmp1,&itmp2,&itmp3,&itmp4,&itmp5,&itmp6);
		Elem_array[i][0]=itmp2;
		Elem_array[i][1]=itmp3;
		Elem_array[i][2]=itmp4;
		Elem_array[i][3]=itmp5;

		Elem_par[i]=itmp6;
	}
	fclose(fp);
}
/////////////////////////////////////////////////////////////////////////////////////
//////////////////////////Count how many nonzeros in each row////////////////////////
void Count(long node,long Elem,long **Elem_array,short int *Count_array,short *n1_point)
{
	long i,i1,i2,i3,i4;
	int n1=0;
	for(i=0;i<Elem;i++)
	{
		i1=Elem_array[i][0]-1;
		i2=Elem_array[i][1]-1;
		i3=Elem_array[i][2]-1;
		i4=Elem_array[i][3]-1;

		Count_array[i1]++;
		Count_array[i2]++;
		Count_array[i3]++;
		Count_array[i4]++;
	}
	for(i=0;i<node;i++)
	{
		if(Count_array[i]>n1)
			n1=Count_array[i];
	}
	*n1_point=n1;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
void Count_edge(long node,long elem,short int max,long **Elem_array,long *Edge,long *Iedge,long *Jedge,long *edge_num)							//max就是n1=4
{
	long i,j,i1,i2,i3,i4,count,*Jtmp;
	
	Jtmp=Arraylong1(max);
	count=0;

	for(i=0;i<Iedge[node];i++)
		Jedge[i]=-1;

	for(i=0;i<elem;i++)
	{
		i1=Elem_array[i][0]-1;
		i2=Elem_array[i][1]-1;
		i3=Elem_array[i][2]-1;
		i4=Elem_array[i][3]-1;
		
		if(i1<i2)
			Edge_assemble(i1,i2,max,Jtmp,Iedge,Jedge);
		else
			Edge_assemble(i2,i1,max,Jtmp,Iedge,Jedge);
		
		if(i1<i4)
			Edge_assemble(i1,i4,max,Jtmp,Iedge,Jedge);
		else
			Edge_assemble(i4,i1,max,Jtmp,Iedge,Jedge); 

		if(i2<i3)
			Edge_assemble(i2,i3,max,Jtmp,Iedge,Jedge);
		else
			Edge_assemble(i3,i2,max,Jtmp,Iedge,Jedge);
	
		if(i3<i4)
			Edge_assemble(i3,i4,max,Jtmp,Iedge,Jedge);
		else
			Edge_assemble(i4,i3,max,Jtmp,Iedge,Jedge);

	}
	free(Jtmp);
	count=0;
	for(i=0;i<node;i++)
	{
		for(j=Iedge[i];j<Iedge[i+1];j++)
			if(Jedge[j]<0)
				count++;
		Edge[i]=count;
	}
	*edge_num=Iedge[node]-count;
}
////////////////////////////////////////////////////////////////////////////////////////
void Edge_assemble(long n1,long n2,int max,long *Jtmp,long *Istiff,long *Jstiff)
{
	int flag,flag1,flag2;
	long i,k;
	flag=0;
	flag1=0;
	flag2=0;

	for(i=0;i<max;i++)
		Jtmp[i]=-1;

	for(i=Istiff[n1];i<Istiff[n1+1];i++)
		if(Jstiff[i]>-1)
		{
			flag=1;
			break;
		}
	
	if(flag==1)
	{
		for(i=Istiff[n1];i<Istiff[n1+1];i++)
		{
			if(Jstiff[i]==n2)
			{
				flag1=1;
				break;
			}
		}	
		if(flag1==0)
		{
			flag2=0;
			for(i=Istiff[n1];i<Istiff[n1+1];i++)
			{
				if(n2<Jstiff[i])
				{
					flag2=1;
					for(k=i;k<Istiff[n1+1]-1;k++)
						Jtmp[k-i]=Jstiff[k];
					
					Jstiff[i]=n2;

					for(k=i;k<Istiff[n1+1]-1;k++)
						Jstiff[k+1]=Jtmp[k-i];
					
					break;
				}
			}		
			if(flag2==0)
			{
				for(i=Istiff[n1];i<Istiff[n1+1];i++)
					if(Jstiff[i]==-1)
					{	
						Jstiff[i]=n2;
						break;
					}
			}
		}
	}
	else
		Jstiff[Istiff[n1]]=n2;
}
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////count node for csr//////////////////////////////
void Count_node_csr(short num,long elem,long **Elem_array,short *Count_array,short *Count_array_add)
{
	long i,j,i1,i2,i3,i4,T;
	int N2;

	N2=num*num;					//36
	T=4*num-5;					//19

	for(i=0;i<elem;i++)
	{
		i1=Elem_array[i][0]-1;
		i2=Elem_array[i][1]-1;
		i3=Elem_array[i][2]-1;
		i4=Elem_array[i][3]-1;

//		printf("%d,%d,%d,%d\n",Count_array[i1],Count_array[i2],Count_array[i3],Count_array[i4]);  //4 or 2

		Count_array_add[i1]=N2*Count_array[i1]-(Count_array[i1]-1)*num;
		Count_array_add[i2]=N2*Count_array[i2]-(Count_array[i2]-1)*num;
		Count_array_add[i3]=N2*Count_array[i3]-(Count_array[i3]-1)*num;
		Count_array_add[i4]=N2*Count_array[i4]-(Count_array[i4]-1)*num;
		
//		printf("%d,%d,%d,%d\n",Count_array_add[i1],Count_array_add[i2],Count_array_add[i3],Count_array_add[i4]);//66 or 126

		for(j=4;j<=T;j++)
		{
			i1=Elem_array[i][j]-1;
			Count_array_add[i1]=2*N2-num;			//上面插值的4到19列全赋值66
		}

		for(j=4*num-4;j<N2;j++)
		{
			i1=Elem_array[i][j]-1;
			Count_array_add[i1]=N2;				//20到36列全为36
		}
	}
}
//////////////////////////////////resort nodes in each elements/////////////////////////
////////////////////////////////////对元素中的节点重新记录///////////////////////////////
void Resort_node_element(short num,long elem,long **Elem_array)
{
	long i,j,j1,j2,n1,n2,n3,n4,*sort;

	n1=0;
	n2=num*(num-1);							//6*5	30
	n3=num*num-1;							//6*6-1	35
	n4=num-1;							//6-1	5

	sort=Arraylong1(num*num);
	
	for(i=0;i<elem;i++)
	{
		sort[n1]=Elem_array[i][0];
		sort[n2]=Elem_array[i][1];
		sort[n3]=Elem_array[i][2];
		sort[n4]=Elem_array[i][3];

		for(j=4;j<=num+1;j++)					//4--7
			{
			sort[num*(j-3)]=Elem_array[i][j];			//6 12 18 24
//			printf("%d\n",num*(j-3));
//			getchar();
			}
		for(j=num+2;j<=2*num-1;j++)				//8--11
			{			
			sort[num*num-2*num-1+j]=Elem_array[i][j];		//31 32 33 34
//			printf("%d\n",num*num-2*num-1+j);
//			getchar();
			}		
		for(j=2*num;j<=3*num-3;j++)				//12--15
			{
			sort[(3*num-1-j)*num-1]=Elem_array[i][j];		//29 23 17 11
//			printf("%d\n",(3*num-1-j)*num-1);
//			getchar();
			}
		for(j=3*num-2;j<=4*num-5;j++)				//16--19
			{
			sort[4*num-4-j]=Elem_array[i][j];			//4 3 2 1
//			printf("%d\n",4*num-4-j);
//			getchar();
			}
		
		for(j1=1;j1<=num-2;j1++)				//1--4
			for(j2=1;j2<=num-2;j2++)			//1--4
				sort[j1*num+j2]=Elem_array[i][4*num-5+(num-2)*(j1-1)+j2];//7--28 20-35
		
		for(j=0;j<num*num;j++)					//0--35
			Elem_array[i][j]=sort[j];


		////////////////Change///////////////
		
	}
	free(sort);
}
/////////////////////////////////////count nodes in computational domaian/////////////////////////
/////////////////////////////////////计算各个区域中有多少计算点////////////////////////////////
void Count_material_num(short num,short Pflag,long node_add,long elem,long **Elem_array,short *Elem_par,long *node1_total,long *node_num,long *elem_num)
{
	short N2;
	long i,j,elem_num1,node_num1,*id;
	
	N2=num*num;
	elem_num1=0;
	node_num1=0;
	id=Arraylong1(N2);

	for(i=0;i<elem;i++)
	{
		if(Elem_par[i]!=Pflag)				//判定区域 Pflag是pml的标志
		{
			for(j=0;j<N2;j++)			//0--36
			{
				id[j]=Elem_array[i][j]-1;
				node1_total[id[j]]=id[j]+1;	//记录节点 那个值的位置就是值的大小
			}
			elem_num1++;
		}
	}

	for(i=0;i<node_add;i++)
	{
//		printf("%d\n",node1_total[i]);
		if(node1_total[i]>0)
			node_num1++;				//非pml层有多少节点

	}
	*elem_num=elem_num1;
	*node_num=node_num1;
	free(id);
}

void Count_material_pml(short num,short Pflag,long node_add,long elem,long **Elem_array,short *Elem_par,long *node1_total,long *node_num,long *elem_num)
{
	short N2;
	long i,j,elem_num1,node_num1,*id;
	
	N2=num*num;
	elem_num1=0;
	node_num1=0;
	id=Arraylong1(N2);

	for(i=0;i<elem;i++)
	{
		if(Elem_par[i]==Pflag)
		{
			for(j=0;j<N2;j++)
			{
				id[j]=Elem_array[i][j]-1;
				node1_total[id[j]]=id[j]+1;
			}
			elem_num1++;
		}
	}

	for(i=0;i<node_add;i++)
	{
		if(node1_total[i]>0)
		{
			node_num1++;
		}
	}
	*elem_num=elem_num1;
	*node_num=node_num1;
	free(id);
}

/////////////////////////////////////Element in computational domaian/////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
void Elem_computation(short Pflag,long elem,short *Elem_par,long *Elem1_count)
{
	long i,count;

	count=0;
	for(i=0;i<elem;i++)
	{
		if(Elem_par[i]!=Pflag)
		{
			Elem1_count[count]=i+1;				//记录不在pml层中的元素编号
			count++;
		}
//		else
//		{
//			printf("pml %10d %10d\n", i, Elem1_count[i]);
//		}
	}
}

void Elem_pml(short Pflag,long elem,short *Elem_par,long *Elem2_count)
{
	long i,count;

	count=0;
	for(i=0;i<elem;i++)
	{
		if(Elem_par[i]==Pflag)
		{
			Elem2_count[count]=i+1;				//记录pml层中的元素编号
			count++;
		}
	}
}
///////////////////////////////////////////////////////////////////////////////
////////////////////////////////different domain sort///////////////////////////
void node_sort(long node_add,long node_num,long *node1_total,long *node1_array)
{
	long i,k;
	k=0;

	for(i=0;i<node_add;i++)
	{
		if(node1_total[i]>0)
		{
			node1_array[k]=node1_total[i];
			k++;
		}
	}

	for(i=0;i<node_add;i++)
	{
		node1_total[i]=0;
	}
	
	for(i=0;i<node_num;i++)						//将所有的计算点做好排序
	{
		k=node1_array[i]-1;
		node1_total[k]=i+1;
		// printf("%d%10d\n", k, node1_total[k]);
	}
//	getchar();
}

/////////////////////////////Two-Domain-Collection/////////////////////////////////
void Collection(long node_add,long *node1_total,long *node2_total,long *node3_total,long *cnode)
{
	long i,count;
	
	count=0;
	for(i=0;i<node_add;i++)
	{
		if(node1_total[i]>0&&node2_total[i]>0)
		{
			node3_total[i]=i+1;
			count++;
		}
	}
	*cnode=count;

}
////////////////////////////////Outside Boundary////////////////////////////////////
void Outside_boundary(long node_add,float X,float Z,double Layx,double Layz,double len,double **node_array,
	long *node4_num,long *node4_total)
{
	long i,count;
	double dx,dz,tmpx1,tmpx2,tmpz1;
	tmpx1=-Layx+len;					//pml层的厚度+单元尺寸
	tmpx2=X+Layx-len;					//计算域的厚度+pml层的厚度+单元尺寸
	tmpz1=Z+Layz-len;
	
	count=0;
	for(i=0;i<node_add;i++)
	{
		dx=node_array[i][0];
		dz=node_array[i][1];
		if(!(dx>tmpx1&&dx<tmpx2&&dz<tmpz1))		//边界
		{
			node4_total[i]=i+1;
			count++;
		}
	}
	*node4_num=count;
}

