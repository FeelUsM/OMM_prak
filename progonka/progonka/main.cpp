#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _MSC_VER
#pragma warning(disable:4996)
#endif

double alp_l, bet_l, alp_r, bet_r;
double * A, * B, * C, * F, * d, * dd, * y;
int N;

void printmas(char * s, double * m)
{
	printf("%s\n",s);
	for(int i=0; i<=N; i++)
		printf("%f ",m[i]);
	printf("\n");
}

int main(int argc, char * argv[])
{
	char path[10000];
	if(argc==1)
		strcpy(path,"../");
	else
		strcpy(path,argv[1]);
	char * filename=path+strlen(path);
	strcpy(filename,"input.txt");
	FILE * in = fopen(path,"r");
	assert(in!=0);

	fscanf(in,"%d %lf %lf %lf %lf",&N,&alp_l,&bet_l,&alp_r,&bet_r);
	A=new double[N+1];
	B=new double[N+1];
	C=new double[N+1];
	F=new double[N+1];
	d=new double[N+1];
	dd=new double[N+1];
	y=new double[N+1];
	int count;
	for(int i=0; i<=N; i++)
		assert(fscanf(in,"%lf",A+i)==1);
	 //printmas("A:",A);
	for(int i=0; i<=N; i++)
		assert(fscanf(in,"%lf",B+i)==1);
	for(int i=0; i<=N; i++)
		assert(fscanf(in,"%lf",C+i)==1);
	 //printmas("C:",C);
	 //printmas("B:",B);
	for(int i=0; i<=N; i++)
		assert(fscanf(in,"%lf",F+i)==1);
	 //printmas("F:",F);

	printf("---------\n");
	d[1]=alp_l; dd[1]=bet_l;
	for(int i=1; i<N; i++)
	{
		double zn=C[i]-A[i]*d[i];
		d[i+1]=B[i]/zn;
		dd[i+1]=(A[i]*dd[i]-F[i])/zn;
		printf("%d: %f, %f %f %f --> %d: %f %f\n",i,A[i],-C[i],B[i],F[i],i+1,d[i+1],dd[i+1]);
	}
	printmas("===========\nd",d);
	printmas("dd",dd);

	printf("---------\n");
	y[N]=(bet_r+dd[N]*alp_r)/(1-d[N]*alp_r);
	for(int i=N; i>0; i--)
	{
		y[i-1]=d[i]*y[i]+dd[i];
		printf("%d: %f %f %f --> %d: %f\n",i,y[i],d[i],dd[i],i-1,y[i-1]);
	}

	strcpy(filename,"output.txt");
	FILE * out= fopen(path,"w");
	for(int i=0; i<=N; i++)
		fprintf(out,"%lf ",y[i]);
	fclose(out);
	//system("pause");
	return 0;
}