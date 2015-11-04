#include <stdio.h>
#include <math.h>
#include <string.h>

#ifdef _MSC_VER
#pragma warning(disable:4996)
#endif
#define ifnot(expr) if(!(expr))

class array2d
{
	double * p;
	int Nx;
public:
	array2d()	{}
	array2d(double * P, int N)	:p(P),Nx(N)	{}
	double * operator[](int y)
	{	return p+Nx*y;	}
};

double X0, T0, hx, tau, epsilon;
int Nx, Nt, maxiter;
int Niter, Ncycle;

/*
du/dt+F(u,x,t)*du/dx=0
(0)
	(!!^v[i+1]!!-v[i+1])/tau
+	(^v[i]-v[i])/tau
+	(!!^v[i+1]!!-^v[i])/tau	*F((!!^v[i+1]!!-^v[i])/2,x[i+1/2],^t)
+	(v[i+1]-v[i])/tau		*F((v[i+1]-v[i])/2,x[i+1/2],t)
=0
(1)
*/
double F(double u, double x, double t)
{	return 2*u + t;	}
double G(double u, double x, double t)
{	return u*(t+u);	}
/*
double (*Tsolver)(double x, double t, double yn, double ynp1, double yn_)
x,t - координаты y[i]
yn		y[i]
ynp1	y[i+1]
yn_		^y[i]
return 	!!^y[i+1]!!
*/
typedef double (*Tsolver)(double,double,double,double,double);
/*
(1)преобразуется как
!!y!!=expr(!!y!!,...)
и решается получившееся ур-е
//зацикливается, если d expr/d!!y!! <=-1
*/
double T2solver(double x, double t, double yn, double ynp1, double yn_)
{
	int i=0;
	double yold , ynew = (ynp1+yn_)/2;
	do
	{
		yold = ynew;
		ynew = ynp1+yn-yn_
			+tau/hx*F((yold+yn_)/2,x+hx/2,t+tau)*(yold-yn_)
			+tau/hx*F((ynp1+yn )/2,x+hx/2,t    )*(ynp1-yn );//expr(!!y!!,...)
		if(++i>maxiter)	throw "zaciklivanie";
		Niter++;
	}while(fabs(ynew-yold)>epsilon);
	Ncycle++;
	return ynew;
}
/*
(0) преобразуется как 
du/dt-dG(u,x,t)/dx=0
Q(!!y!!) = 	(^y[i]-y[i]+!!^y[i+1]!!-y[i+1])/tau -
		-	(G(!!^y[i+1]!!,x[i+1],^t)-^G[i]+G[i+1]-G[i])/h = 0
*/
double Newton(double x, double t, double yn, double ynp1, double yn_)
{
	int i=0;
	double yold , ynew = (ynp1+yn_)/2;
	do
	{
		yold = ynew;
		double Gold	=G(yold,t+tau,x+hx);
		double Gn_	=G(yn_ ,t+tau,x   );
		double Gnp1	=G(ynp1,t    ,x+hx);
		double Gn	=G(yn  ,t    ,x   );
		double Q	=(yn_-yn+yold-ynp1)/tau-(Gold-Gn_+Gnp1-Gn)/hx;
		double Q_	=1/tau-F(yold,x+hx,t+tau)/hx;
		ynew = yold - Q/Q_;
		if(++i>maxiter)	throw "zaciklivanie";
		Niter++;
	}while(fabs(ynew-yold)>epsilon);
	Ncycle++;
	return ynew;
}
void main_cycle(array2d U, Tsolver solver)
{
	Niter=Ncycle=0;
	for(int it=0; it<Nt; it++)
		for(int ix=0; ix<Nx; ix++)
			U[it+1][ix+1] = solver(X0+ix*hx,T0+it*tau, U[it][ix],U[it][ix+1],U[it+1][ix]);
}

int main(int argc, char * argv[])
{
try{
	//printf("hello\n");
	if(argc>1)
		for(int i=0; i<argc; i++)
			printf("%d: '%s'\n",i,argv[i]);

	char * nameinput, * nameoutput;
	nameinput= (argc>=2) ? argv[1] : "../input.txt";
	nameoutput=(argc>=3) ? argv[2] : "../output.txt";

	FILE * fin=fopen(nameinput,"r");
	ifnot(fin) throw "cannot open input file";
	ifnot(8==fscanf(fin,"%lg %lg %d %lg %lg %d %lg %d",&X0,&hx,&Nx,&T0,&tau,&Nt,&epsilon,&maxiter))
		throw "cannot read start parameters from input";

	double * array=new double[(Nx+1)*(Nt+1)];

	for(int it=0; it<=Nt; it++)
		ifnot(1==fscanf(fin,"%lg",array+it*(Nx+1)))
			throw "cannot read time gran usl";
	double tmp;

	ifnot(1==fscanf(fin,"%lg",&tmp))
		throw "cannot start read space gran usl";
	ifnot(tmp==array[0])
		throw "gran usl must be agree in (T0,X0)";

	for(int ix=1; ix<=Nx; ix++)
		ifnot(1==fscanf(fin,"%lg",array+ix))
			throw "cannot read space gran usl";

	if(argc==4 && 0==strcmp("--t2",argv[3]))
		main_cycle(array2d(array,Nx+1),T2solver);
	else
		main_cycle(array2d(array,Nx+1),Newton);

	FILE * fout=fopen(nameoutput,"w");
	ifnot(fout)	throw "cannot open output file";
	double sriters=double(Niter)/Ncycle;
	fprintf(fout,"%f\n",sriters);
	for(int it=0; it<=Nt; it++)
	{
		for(int ix=0; ix<=Nx; ix++)
			fprintf(fout,"%g ",*array++);
		fprintf(fout,"\n");
	}
}
catch(const char * s)
{
	printf("ERROR:\n");
	printf("%s\n",s);
	if(argc==1)	getchar();//"pause");
	return 16;
}
catch(...)
{
	printf("UNKNOWN ERROR\n");
	if(argc==1)	getchar();//"pause");
	return 9;
}

	return 0;
}
