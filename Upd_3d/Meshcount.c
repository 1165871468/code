#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include "Meshcount.h"
/////////////////////////////////re_node_trian///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void re_node_elem(short int num,long node,long elem,long edge_num,long **Elem_array,double **node_arrayt,double **node_array,long *Edge,long *Iedge,long *Jedge)	//num插值阶数；插值
{
	long i,i1,i2,i3,i4,j,j1,j2,count,num1,star,id,id1,tmp;
	double *P,x1,x2,z1,z2;

	P=Arraydouble1(num);
	Position(num-1,P);
	star=node+(num-2)*edge_num;

	for(i=0;i<node;i++)
	{
		node_array[i][0]=node_arrayt[i][0];
		node_array[i][1]=node_arrayt[i][1];
	}
	
	for(i=0;i<elem;i++)
	{
		i1=Elem_array[i][0]-1;
		i2=Elem_array[i][1]-1;
		i3=Elem_array[i][2]-1;
		i4=Elem_array[i][3]-1;

		/////////i1<i2//////////////
		if(i1<i2)
		{
			count=0;
			for(j=Iedge[i1];j<Iedge[i1+1];j++)
			{
				count++;
				if(Jedge[j]==i2)
					break;
			}
//			printf("%d\n",count);
			if(i1>0)
				tmp=Edge[i1-1];
			else
				tmp=0;

			num1=Iedge[i1]+count-tmp;			//插值后的边
			id=node+(num-2)*(num1-1);			//插值的点
//			printf("%d,%d,%d,%d,%d\n",num1,Iedge[i1],count,tmp,id);

			for(j=1;j<=num-2;j++)
			{
				Elem_array[i][3+j]=id+j;					//4--7
//				printf("%d,%d,%d\n",i,3+j,Elem_array[i][3+j]);
			}
			for(j=1;j<=num-2;j++)
			{
				node_array[id+j-1][0]=(float)((1.0-P[j])/2.0*node_array[i1][0]+(P[j]+1.0)/2.0*node_array[i2][0]);
				node_array[id+j-1][1]=(float)((1.0-P[j])/2.0*node_array[i1][1]+(P[j]+1.0)/2.0*node_array[i2][1]);						//插值点
			}
		}
		/////////i1>i2//////////////
		else
		{
			count=0;
			for(j=Iedge[i2];j<Iedge[i2+1];j++)
			{
				count++;
				if(Jedge[j]==i1)
					break;
			}
			if(i2>0)
				tmp=Edge[i2-1];
			else
				tmp=0;

			num1=Iedge[i2]+count-tmp;
			id=node+(num-2)*(num1-1);

			for(j=1;j<=num-2;j++)
				Elem_array[i][3+j]=id+num-1-j;

			for(j=1;j<=num-2;j++)
			{
				node_array[id+j-1][0]=(float)((1.0-P[j])/2.0*node_array[i2][0]+(P[j]+1.0)/2.0*node_array[i1][0]);
				node_array[id+j-1][1]=(float)((1.0-P[j])/2.0*node_array[i2][1]+(P[j]+1.0)/2.0*node_array[i1][1]);
			}		
		}
		/////////i2<i3//////////////
		if(i2<i3)
		{
			count=0;
			for(j=Iedge[i2];j<Iedge[i2+1];j++)
			{
				count++;
				if(Jedge[j]==i3)
					break;
			}
			if(i2>0)
				tmp=Edge[i2-1];
			else
				tmp=0;

			num1=Iedge[i2]+count-tmp;
			id=node+(num-2)*(num1-1);

			for(j=1;j<=num-2;j++)	
				Elem_array[i][num+1+j]=id+j;				//8--11
			
			for(j=1;j<=num-2;j++)
			{
				node_array[id+j-1][0]=(float)((1.0-P[j])/2.0*node_array[i2][0]+(P[j]+1.0)/2.0*node_array[i3][0]);
				node_array[id+j-1][1]=(float)((1.0-P[j])/2.0*node_array[i2][1]+(P[j]+1.0)/2.0*node_array[i3][1]);
			}
		}
		/////////i2>i3//////////////
		else
		{
			count=0;
			for(j=Iedge[i3];j<Iedge[i3+1];j++)
			{
				count++;
				if(Jedge[j]==i2)
					break;
			}
			if(i3>0)
				tmp=Edge[i3-1];
			else
				tmp=0;

			num1=Iedge[i3]+count-tmp;
			id=node+(num-2)*(num1-1);
			
			for(j=1;j<=num-2;j++)
				Elem_array[i][num+1+j]=id+num-1-j;

			for(j=1;j<=num-2;j++)
			{
				node_array[id+j-1][0]=(float)((1.0-P[j])/2.0*node_array[i3][0]+(P[j]+1.0)/2.0*node_array[i2][0]);
				node_array[id+j-1][1]=(float)((1.0-P[j])/2.0*node_array[i3][1]+(P[j]+1.0)/2.0*node_array[i2][1]);
			}
		}
		/////////i3<i4//////////////
		if(i3<i4)
		{
			count=0;
			for(j=Iedge[i3];j<Iedge[i3+1];j++)
			{
				count++;
				if(Jedge[j]==i4)
					break;
			}
			if(i3>0)
				tmp=Edge[i3-1];
			else
				tmp=0;

			num1=Iedge[i3]+count-tmp;
			id=node+(num-2)*(num1-1);

			for(j=1;j<=num-2;j++)
				Elem_array[i][2*num-1+j]=id+j;				//12-15
				
			for(j=1;j<=num-2;j++)
			{
				node_array[id+j-1][0]=(float)((1.0-P[j])/2.0*node_array[i3][0]+(P[j]+1.0)/2.0*node_array[i4][0]);
				node_array[id+j-1][1]=(float)((1.0-P[j])/2.0*node_array[i3][1]+(P[j]+1.0)/2.0*node_array[i4][1]);
			}
		}
		/////////i3>i4//////////////
		else
		{
			count=0;
			for(j=Iedge[i4];j<Iedge[i4+1];j++)
			{
				count++;
				if(Jedge[j]==i3)
					break;
			}
			if(i4>0)
				tmp=Edge[i4-1];
			else
				tmp=0;

			num1=Iedge[i4]+count-tmp;
			id=node+(num-2)*(num1-1);
			
			for(j=1;j<=num-2;j++)
				Elem_array[i][2*num-1+j]=id+num-1-j;

			for(j=1;j<=num-2;j++)
			{
				node_array[id+j-1][0]=(float)((1.0-P[j])/2.0*node_array[i4][0]+(P[j]+1.0)/2.0*node_array[i3][0]);
				node_array[id+j-1][1]=(float)((1.0-P[j])/2.0*node_array[i4][1]+(P[j]+1.0)/2.0*node_array[i3][1]);
			}
		}

		/////////i4<i1//////////////
		if(i4<i1)
		{
			count=0;
			for(j=Iedge[i4];j<Iedge[i4+1];j++)
			{
				count++;
				if(Jedge[j]==i1)
					break;
			}
			if(i4>0)
				tmp=Edge[i4-1];
			else
				tmp=0;

			num1=Iedge[i4]+count-tmp;
			id=node+(num-2)*(num1-1);

			for(j=1;j<=num-2;j++)
				Elem_array[i][3*num-3+j]=id+j;						//16--19
				
			for(j=1;j<=num-2;j++)
			{
				node_array[id+j-1][0]=(float)((1.0-P[j])/2.0*node_array[i4][0]+(P[j]+1.0)/2.0*node_array[i1][0]);
				node_array[id+j-1][1]=(float)((1.0-P[j])/2.0*node_array[i4][1]+(P[j]+1.0)/2.0*node_array[i1][1]);
			}
		}
		/////////i4>i1//////////////
		else
		{
			count=0;
			for(j=Iedge[i1];j<Iedge[i1+1];j++)
			{
				count++;
				if(Jedge[j]==i4)
					break;
			}
			if(i1>0)
				tmp=Edge[i1-1];
			else
				tmp=0;

			num1=Iedge[i1]+count-tmp;
			id=node+(num-2)*(num1-1);
			
			for(j=1;j<=num-2;j++)
				Elem_array[i][3*num-3+j]=id+num-1-j;

			for(j=1;j<=num-2;j++)
			{
				node_array[id+j-1][0]=(float)((1.0-P[j])/2.0*node_array[i1][0]+(P[j]+1.0)/2.0*node_array[i4][0]);
				node_array[id+j-1][1]=(float)((1.0-P[j])/2.0*node_array[i1][1]+(P[j]+1.0)/2.0*node_array[i4][1]);
			}
		}
////////////////////////////////center nodes of element////////////////////////////////////
/////////////////////////////////////元素的内部////////////////////////////////////////
		id=star+(num-2)*(num-2)*i;				//内部插值的点数

		for(j1=1;j1<=num-2;j1++)
		{
			id1=(j1-1)*(num-2);

			i1=Elem_array[i][3+j1]-1;
			i2=Elem_array[i][3*num-2-j1]-1;
			
			x1=node_array[i1][0];
			x2=node_array[i2][0];

			z1=node_array[i1][1];
			z2=node_array[i2][1];

			for(j2=1;j2<=num-2;j2++)			//20--35
			{
				Elem_array[i][4*num-5+id1+j2]=id+id1+j2;
				node_array[id+id1+j2-1][0]=(float)((1.0-P[j2])/2.0*x1+(P[j2]+1.0)/2.0*x2);
				node_array[id+id1+j2-1][1]=(float)((1.0-P[j2])/2.0*z1+(P[j2]+1.0)/2.0*z2);
//			printf("%d\n",id+id1+j2-1);
			}
		}
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////Count triangle in computational domain///////////////////////////

void Elem_in_count(short num,long elem,long elem_num,double dis,double X,double Z,long **Elem_array,double **node_array,long *elem_num1,long *Elem1_total)
{
	double x1,z1,x2,z2,x3,z3,x4,z4;
	long i,i1,i2,i3,i4,count;
	int n1,n2,n3,n4;
	
	for(i=0;i<elem_num;i++)
	{
		Elem1_total[i]=i+1;
	}

	n1=0;
	n2=num*(num-1);
	n3=num*num-1;
	n4=num-1;

	for(i=elem_num;i<elem;i++)
	{
		i1=Elem_array[i][n1]-1;
		i2=Elem_array[i][n2]-1;
		i3=Elem_array[i][n3]-1;
		i4=Elem_array[i][n4]-1;

		x1=node_array[i1][0];
		z1=node_array[i1][1];
		
		x2=node_array[i2][0];
		z2=node_array[i2][1];
		
		x3=node_array[i3][0];
		z3=node_array[i3][1];

		x4=node_array[i4][0];
		z4=node_array[i4][1];
		
		if((x1<=0&&x1>-dis)&&(z1>=0&&z1<=Z+dis))
		{
			Elem1_total[i]=i+1;
			continue;
		}
		else if((x1>0&&x1<X)&&(z1>=Z&&z1<=Z+dis))
		{
			Elem1_total[i]=i+1;
			continue;
		}
		else if((x1>=X&&x1<=X+dis)&&(z1>=0&&z1<=Z+dis))
		{
			Elem1_total[i]=i+1;
			continue;
		}

		if((x2<=0&&x2>-dis)&&(z2>=0&&z2<=Z+dis))
		{
			Elem1_total[i]=i+1;
			continue;
		}
		else if((x2>0&&x2<X)&&(z2>=Z&&z2<=Z+dis))
		{
			Elem1_total[i]=i+1;
			continue;
		}
		else if((x2>=X&&x2<=X+dis)&&(z2>=0&&z2<=Z+dis))
		{
			Elem1_total[i]=i+1;
			continue;
		}

		if((x3<=0&&x3>-dis)&&(z3>=0&&z3<=Z+dis))
		{
			Elem1_total[i]=i+1;
			continue;
		}
		else if((x3>0&&x3<X)&&(z3>=Z&&z3<=Z+dis))
		{
			Elem1_total[i]=i+1;
			continue;
		}
		else if((x3>=X&&x3<=X+dis)&&(z3>=0&&z3<=Z+dis))
		{
			Elem1_total[i]=i+1;
			continue;
		}

		if((x4<=0&&x4>-dis)&&(z4>=0&&z4<=Z+dis))
		{
			Elem1_total[i]=i+1;
			continue;
		}
		else if((x4>0&&x4<X)&&(z4>=Z&&z4<=Z+dis))
		{
			Elem1_total[i]=i+1;
			continue;
		}
		else if((x4>=X&&x4<=X+dis)&&(z4>=0&&z4<=Z+dis))
		{
			Elem1_total[i]=i+1;
			continue;
		}
	}
	count=0;
	for(i=elem_num;i<elem;i++)
	{
		if(Elem1_total[i]>0)
		count++;
	}
	*elem_num1=count+elem_num;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////count domain triangels to total domain/////////////////////
void Elem_to_total(long elem,long *Elem1_total,long *Elem1_count)
{
	long i,count;
	count=0;
	for(i=0;i<elem;i++)
	{
		if(Elem1_total[i]>0)
		{
			Elem1_count[count]=Elem1_total[i];
			count++;
		}
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
void node2_count(short num,long node_add,long elem_num1,long *Elem1_count,long **Elem_array,long *node2_total,long *node2_num)
{
	long i,j,k,*id,count,N2;

	N2=num*num;
	id=Arraylong1(N2);
	for(i=0;i<elem_num1;i++)
	{
		j=Elem1_count[i]-1;
		for(k=0;k<N2;k++)
		{
			id[k]=Elem_array[j][k]-1;
			node2_total[id[k]]=id[k]+1;
		}
	}
	count=0;
	for(i=0;i<node_add;i++)
	{
		if(node2_total[i]>0)
		{
			count++;
		}
	}
	*node2_num=count;
	free(id);
}
//////////////////////////////////////////////////////////////////////////////////////////////
void node_count_num(long node_add,long node_num,long *node2_total,long *node2_array)
{
	long i,j,count;
	count=0;
	for(i=0;i<node_add;i++)
	{
		if(node2_total[i]>0)
		{
			node2_array[count]=node2_total[i];
			count++;
		}
	}

	for(i=0;i<node_add;i++)
	{
		node2_total[i]=0;
	}

	for(i=0;i<node_num;i++)
	{
		j=node2_array[i]-1;
		node2_total[j]=i+1;
	}
}
///////////////////////////////////////////CSR Initial//////////////////////////////
void Stiff_initial(int Node2,short max,long node2_num,long elem_num,long *Elem1_count,long **Elem_array,long *node2_total,long *Istiff,long *Jstiff)
{
	long i,i1,j1,j2,k,*id,index1,index2,index3,index4,*Jtmp,Max;
	
	Max=2*max;
	for(i=0;i<Istiff[2*node2_num];i++)
	{
		Jstiff[i]=-1;
	}

	id=Arraylong1(Node2);
	Jtmp=Arraylong1(Max);

	for(i=0;i<elem_num;i++)
	{
		k=Elem1_count[i]-1;

		for(i1=0;i1<Node2;i1++)
		{
			id[i1]=Elem_array[k][i1]-1;
			id[i1]=node2_total[id[i1]]-1;
		}

		for(j1=0;j1<Node2;j1++)
		{
			index1=2*id[j1];
			index2=index1+1;

			for(j2=0;j2<Node2;j2++)
			{
				index3=2*id[j2];
				index4=index3+1;

				Ini_assemble(Max,index1,index3,Jtmp,Istiff,Jstiff);
				Ini_assemble(Max,index1,index4,Jtmp,Istiff,Jstiff);
				Ini_assemble(Max,index2,index3,Jtmp,Istiff,Jstiff);
				Ini_assemble(Max,index2,index4,Jtmp,Istiff,Jstiff);
			}
		}  
	}
	free(Jtmp);
	free(id);
}

void Ini_assemble(int Max,long n1,long n2,long *Jtmp,long *Istiff,long *Jstiff)
{
	int flag,flag1,flag2;
	long i,n11,k;

	flag=0;
	flag1=0;
	flag2=0;

	n11=n1+1;

	for(i=0;i<Max;i++)
		Jtmp[i]=-1;

	for(i=Istiff[n1];i<Istiff[n11];i++)
		if(Jstiff[i]>-1)
		{
			flag=1;
			break;
		}
	
	if(flag==1)
	{
		for(i=Istiff[n1];i<Istiff[n11];i++)
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
			for(i=Istiff[n1];i<Istiff[n11];i++)
			{
				if(n2<Jstiff[i])
				{
					flag2=1;
					for(k=i;k<Istiff[n11]-1;k++)
					{
						Jtmp[k-i]=Jstiff[k];
					}
					Jstiff[i]=n2;
					for(k=i;k<Istiff[n11]-1;k++)
					{
						Jstiff[k+1]=Jtmp[k-i];
					}
					break;
				}
			}		
			if(flag2==0)
			{
				for(i=Istiff[n1];i<Istiff[n11];i++)
					if(Jstiff[i]==-1)
					{	
						Jstiff[i]=n2;
						break;
					}
			}
		}
	}
	else
	{
		Jstiff[Istiff[n1]]=n2;
	}	
}

void Stiff_count(long node2_num,long *Istiff_ini,long *Jstiff_ini,long *Istiff)
{
	long i,j,i1,node_add2,count;
	
	node_add2=2*node2_num;

	count=0;
	for(i=0;i<node_add2;i++)
	{
		i1=i+1;
		for(j=Istiff_ini[i];j<Istiff_ini[i1];j++)
		{
			if(Jstiff_ini[j]>=0)
				count++;
		}
		Istiff[i1]=count;
	}
}
