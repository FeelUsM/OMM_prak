#include <stdio.h>
#include <stdlib.h>
//#include <assert.h>
#include <math.h>
#include <string.h>
#include <functional>//чтобы из == и < получить !=, <=, >=, >, хотя это все равно не используется
using namespace std;
#define ifnot(expr) if(!(expr))

#ifdef _MSC_VER
#pragma warning(disable:4996)
#endif

template <class T>
class superpointer
{
	typedef superpointer<T> my_t;
	T * p;
	int dN;
public:
	superpointer()	{}
	superpointer(T * P, int N)	:p(P),dN(N)	{}
	T & operator*()	{return *p;}
	my_t & operator++()	{p+=dN; return *this;}
	my_t & operator--() {p-=dN; return *this;}
	my_t & operator++(int) {my_t tmp=*this; (*this)++; return tmp;}
	my_t & operator--(int) {my_t tmp=*this; (*this)--; return tmp;}
	my_t & operator+=(int x)	{p+=x*dN; return *this;}
	my_t & operator-=(int x)	{p-=x*dN; return *this;}
	T & operator[](int x)	{return *(p+x*dN);}//вообще из этого класса используется тольк этот метод
	bool operator==(const my_t & r)const
	{	assert(dN==r.dN);	return p==r.p;	}
	bool operator<(const my_t & r)const
	{	assert(dN==r.dN);	return p<r.p;	}
	int operator-(const my_t & r)const
	{	assert(dN==r.dN);	return (r.p-p)/dN;	}
};
template <class T> inline
superpointer<T> operator+(superpointer<T> l, int r)
{	l+=r; return l;	}
template <class T> inline
superpointer<T> operator+(int r, const superpointer<T> & l)
{	l+=r; return l;	}
template <class T> inline
superpointer<T> operator-(const superpointer<T> & l, int r)
{	l-=r; return l;	}

class superarrayYX
{
	double * p;
public:
	double * pointer()	{	return p;	}
	int _Nx, _Ny;
	superarrayYX()	{}
	superarrayYX(double * P, int x, int y)	:p(P),_Nx(x),_Ny(y)	{}
	double * operator[](int i)
	{	return p+i*(_Nx+1);	}
};
class superarrayXY
{
	double * p;
public:
	double * pointer()	{	return p;	}
	int _Nx, _Ny;
	superarrayXY()	{}
	superarrayXY(double * P, int x, int y)	:p(P),_Ny(x),_Nx(y)	{}
	superpointer<double> operator[](int i)
	{	return superpointer<double>(p+i,_Ny+1);	}
};
//==========================================================
double X0, Y0, T0;//координаты левого нижнего ближайшего угла
double hx, hy, tau;//шаги
int _Nx_, _Ny_, Nt;//размеры
double koef;//см.ниже

//коэффициенты гран-условий
double	alp_xl, bet_xl,		//alp_xl*du/dx - bet_xl*u = mu_xl
		alp_xh, bet_xh,		//alp_xh*du/dx + bet_xh*u = mu_xh
		alp_yl, bet_yl,		//alp_yl*du/dy - bet_yl*u = mu_yl
		alp_yh, bet_yh;		//alp_yh*du/dy + bet_yh*u = mu_yh
//функции с гранусловиями
double * mu_xl, * mu_xh, * mu_yl, * mu_yh;
//функции с гранусловиями, если задан флаг --internal-gu
double fmu_xl(double x, double y, double t)
{	return 0;	}
double fmu_xh(double x, double y, double t)
{	return 0;	}
double fmu_yl(double x, double y, double t)
{	return 0;	}
double fmu_yh(double x, double y, double t)
{	return 0;	}

//d/dt = d/dx(P*d/dx u)+d/dy(Q*d/dy u) + F
double P(double x, double y, double t)
{	return koef;	}
double Q(double x, double y, double t)
{	return koef;	}
double F(double x, double y, double t)
{	return exp(t)*x*y*(1-x);	}
//{	return x*y*(1-x);	}
//================================================
//функции отладки: (выключенные)
double * mainpointer;
template <class P>
void printoff(FILE * f, int N, P p, const char * mes=0, bool stop=false)
{return;
	if(mes)	fprintf(f,"%s:",mes);

	for(int i=0; i<=N; i++)
		fprintf(f,"%d ",&(p[i])-mainpointer);

	if(f==stdout && stop)
		getchar();
	else
		fprintf(f,"\n");
	fflush(f);
}
template <class P>
void printarray(FILE * f, int N, P p, const char * mes=0, bool stop=false)
{return;
	if(mes)	fprintf(f,"%s:",mes);

	for(int i=0; i<=N; i++)
		fprintf(f,"%f ",p[i]);

	if(f==stdout && stop)
		getchar();
	else
		fprintf(f,"\n");
	fflush(f);
}
void printlayer(FILE * f, double * p, const char * mes=0, bool stop=false)
{return;
	if(mes)	fprintf(f,"%s\n",mes);

	for(int iy=0; iy<=_Ny_; iy++)
	{
		for(int ix=0; ix<=_Nx_; ix++)
			fprintf(f,"%f ",*p++);
		fprintf(f,"\n");
	}

	if(f==stdout && stop)
		getchar();
	else
		fprintf(f,"\n");
	fflush(f);
}
template <class sa>
void fillfield(sa arr)
{
	printf("Nx=%d; Ny=%d\n",arr._Nx,arr._Ny);
	int n=0;
	for(int iy=0; iy<=arr._Ny; iy++)
		for(int ix=0; ix<=arr._Nx; ix++)
			arr[iy][ix]=n++;
}
void debugXY()
{
	_Nx_=7;	_Ny_=10;
	double * p=new double[(_Nx_+1)*(_Ny_+1)];
	fillfield(superarrayXY(p,_Nx_,_Ny_));
	printlayer(stdout,p,"debugXY",true);
}
//================================================
/*
A[m]*y[m-1]-C[m]*y[m]+B[m]*y[m+1]=F[m]//m in 1..N-1
y[0]=alp_l*y[1]+bet_l
y[N]=alp_r*y[N-1]+bet_r
*/
template<class PA, class PB, class PC, class PF, class PY>
void progonka(int N, PY y, PA A, PB B, PC C, PF F, double alp_l, double bet_l, double alp_r, double bet_r)
{
	printarray(stdout,N-2,A+1,"A");
	printarray(stdout,N-2,B+1,"B");
	printarray(stdout,N-2,C+1,"C");
	printarray(stdout,N-2,F+1,"F");
	printarray(stdout,N,y,"in y");
	printoff(stdout,N,y,"off y");
	double * d, * dd;
	d=new double[N+1];//в зависимости от направления N меняется
	dd=new double[N+1];//однако здесь все равно еще это можно оптимизировать
	d[1]=alp_l; dd[1]=bet_l;
	for(int i=1; i<N; i++)
	{
		double zn=C[i]-A[i]*d[i];
		d[i+1]=B[i]/zn;
		dd[i+1]=(A[i]*dd[i]-F[i])/zn;
	}
	y[N]=(bet_r+dd[N]*alp_r)/(1-d[N]*alp_r);
	for(int i=N; i>0; i--)
		y[i-1]=d[i]*y[i]+dd[i];
	printarray(stdout,N,y,"out y");
	//printoff(stdout,N,y,"out y off");
	delete d[];
	delete dd[];
}
typedef double (*TF)(double,double,double);
/*
неявная по X, явная по Y
*/
template<class superarray>
void solve_layer(superarray out, superarray in, int it/*для гран-условия*/, double t/*для неоднородности*/,
				 double X0, double Y0, double hx, double hy, TF P, TF Q, TF f,
				 double alp_xl, double bet_xl, double alp_xh, double bet_xh, double * mu_xl, double * mu_xh,
				 double alp_yl, double bet_yl, double alp_yh, double bet_yh, double * mu_yl, double * mu_yh
				 )
{
	mainpointer=out.pointer();//для отладки
	int Nx=in._Nx, Ny=in._Ny;
	//прогонка вдоль x
	double * A, * B, * C, * F;
	A=new double[Nx+1];
	B=new double[Nx+1];
	C=new double[Nx+1];
	F=new double[Nx+1];
	for(int iy=1; iy<Ny; iy++)
	{
		for(int ix=1; ix<Nx; ix++)
		{
			double x=X0+ix*hx, y=Y0+iy*hy;
			A[ix] = P(x-hx/2,y,t) /hx/hx;
			B[ix] = P(x+hx/2,y,t) /hx/hx;
			C[ix] = A[ix] + B[ix] + 2/tau;
			F[ix] = f(x,y,t) + in[iy][ix]*2/tau +
				( Q(x,y-hy/2,t)*in[iy-1][ix] + Q(x,y+hy/2,t)*in[iy+1][ix]
					-(Q(x,y-hy/2,t)+Q(x,y+hy/2,t))*in[iy][ix] )/hy/hy;
				F[ix]=-F[ix];
		}
		progonka(Nx,out[iy],A,B,C,F,
			alp_xl						/(alp_xl+bet_xl*hx),
			-mu_xl[it*(Ny+1)+iy]*hx		/(alp_xl+bet_xl*hx),
			alp_xh						/(alp_xh+bet_xh*hx),
			mu_xh[it*(Ny+1)+iy]*hx		/(alp_xh+bet_xh*hx)
			);
		printlayer(stdout,out.pointer(),"step progonka",true);
	}
	//гран. условия по краям вдоль x
	double zn = alp_yl+bet_yl*hy;
	for(int ix=0; ix<=Nx; ix++)
		out[0][ix] = alp_yl*out[1][ix]/zn - mu_yl[it*(Nx+1)+ix]*hx/zn;
	zn = alp_yh+bet_yh*hy;
	for(int ix=0; ix<=Nx; ix++)
		out[Ny][ix] = alp_yh*out[Ny-1][ix]/zn - mu_yh[it*(Nx+1)+ix]*hx/zn;
	delete A[];
	delete B[];
	delete C[];
	delete F[];
}
double PP(double x, double y, double t)
{	return P(y,x,t);	}
double QQ(double x, double y, double t)
{	return Q(y,x,t);	}
double FF(double x, double y, double t)
{	return F(y,x,t);	}
/*
du/dt=Lu+F
Lu=d/dx(P*du/dx)+d/dy(Q*du/dy)
alp_xl*du/dx-bet_xl|(x=X0)=mu_xl
alp_yl*du/dy-bet_yl|(y=Y0)=mu_yl
alp_xh*du/dx+bet_xh|(x=X0+Nx*hx)=mu_xh
alp_yh*du/dy+bet_yh|(y=Y0+Ny*hy)=mu_yh
u|(t=0)-задано
*/
void main_cycle(double * sloy)
{
	double * promsloy=new double[(_Nx_+1)*(_Ny_+1)];
	for(int i=0; i<Nt; i++)
	{
		solve_layer(superarrayYX(promsloy,_Nx_,_Ny_),superarrayYX(sloy,_Nx_,_Ny_), 2*i, T0+i*tau+tau/2,
			X0, Y0, hx, hy, P, Q, F,
			alp_xl, bet_xl, alp_xh, bet_xh, mu_xl, mu_xh,
			alp_yl, bet_yl, alp_yh, bet_yh, mu_yl, mu_yh);
		printlayer(stdout,promsloy,"after first half layer",true);
		solve_layer(superarrayXY(sloy,_Nx_,_Ny_),superarrayXY(promsloy,_Nx_,_Ny_), 2*i+1, T0+i*tau+tau/2,
			Y0,X0,hy,hx, QQ, PP, FF, 
			alp_yl, bet_yl, alp_yh, bet_yh, mu_yl, mu_yh,
			alp_xl, bet_xl, alp_xh, bet_xh, mu_xl, mu_xh);
		printlayer(stdout,sloy,"after first full layer",true);
	}
	delete promsloy[];
}
//================================================
#define Nx _Nx_
#define Ny _Ny_
int main(int argc, char * argv[])
{
	//debugXY();
try{
	//printf("hello\n");
	if(argc>1)
		for(int i=0; i<argc; i++)
			printf("%d: '%s'\n",i,argv[i]);

	char * namegranusl, * nameinput, * nameoutput;
	namegranusl=(argc>=2) ? argv[1] : "../gran_usl.txt";
	nameinput=	(argc>=3) ? argv[2] : "../input_layer.txt";
	nameoutput=	(argc>=4) ? argv[3] : "../output_layer.txt";

	FILE * fgu=fopen(namegranusl,"r");
	ifnot(fgu!=0) throw "cannot open gran_usl file";

	ifnot(18==fscanf(fgu,"%lf %lf %lf %lf %lf %d %d %lf %lf %d %lf %lf %lf %lf %lf %lf %lf %lf",
		&koef, &X0, &Y0, &hx, &hy, &Nx, &Ny,	//коэффициент и размеры поля
		&T0, &tau, &Nt,							//задание (кол-во циклов)
		&alp_xl, &bet_xl, 
		&alp_xh, &bet_xh,
		&alp_yl, &bet_yl,
		&alp_yh, &bet_yh)) throw "cannot read formating data from gran_usl file";
	mu_xl=new double[(Ny+1)*2*Nt];
	mu_xh=new double[(Ny+1)*2*Nt];
	mu_yl=new double[(Nx+1)*2*Nt];
	mu_yh=new double[(Nx+1)*2*Nt];
	if(argc==5 && 0==strcmp(argv[4],"--interface-gu") || argc==1)
	{
		double * p=mu_xl;
		for(int it=0; it<2*Nt; it++)
			for(int iy=0; iy<=Ny; iy++)
				ifnot(1==fscanf(fgu,"%lf",p++)) throw "cannot start read 'gran usl's in gran_usl file";
		p=mu_xh;
		for(int it=0; it<2*Nt; it++)
			for(int iy=0; iy<=Ny; iy++)
				ifnot(1==fscanf(fgu,"%lf",p++)) throw "cannot read 'gran usl's in gran_usl file";
		p=mu_yl;
		for(int it=0; it<2*Nt; it++)
			for(int ix=0; ix<=Nx; ix++)
				ifnot(1==fscanf(fgu,"%lf",p++)) throw "cannot read 'gran usl's in gran_usl file";
		p=mu_yh;
		for(int it=0; it<2*Nt; it++)
			for(int ix=0; ix<=Nx; ix++)
				ifnot(1==fscanf(fgu,"%lf",p++)) throw "cannot finish read 'gran usl's in gran_usl file";
	}
	else
	{
		double * p=mu_xl;
		for(int it=1; it<=2*Nt; it++)
			for(int iy=0; iy<=Ny; iy++)
				*p++= fmu_xl(X0,Y0+iy*hy,T0+it*tau/2);
		p=mu_xh;
		for(int it=1; it<=2*Nt; it++)
			for(int iy=0; iy<=Ny; iy++)
				*p++= fmu_xh(X0+Nx*hx,Y0+iy*hy,T0+it*tau/2);
		p=mu_yl;
		for(int it=1; it<=2*Nt; it++)
			for(int ix=0; ix<=Nx; ix++)
				*p++= fmu_yl(X0+ix*hx,Y0,T0+it*tau/2);
		p=mu_yh;
		for(int it=1; it<=2*Nt; it++)
			for(int ix=0; ix<=Nx; ix++)
				*p++= fmu_yh(X0+ix*hx,Y0+Ny*hy,T0+it*tau/2);
	}

	double * sloy=new double[(Nx+1)*(Ny+1)];

	FILE * fin=fopen(nameinput,"r");
	ifnot(fin!=0) throw "cannot open input_layer file";
	double * p=sloy;
	for(int iy=0; iy<=Ny; iy++)
		for(int ix=0; ix<=Nx; ix++)
			ifnot(1==fscanf(fin,"%lf",p++)) throw "cannot read from input_layer file";

	main_cycle(sloy);

	FILE * fout=fopen(nameoutput,"w");
	ifnot(fout!=0) throw "cannot open to write output_layer file";
	p=sloy;
	for(int iy=0; iy<=Ny; iy++)
	{
		for(int ix=0; ix<=Nx; ix++)
			fprintf(fout,"%lg ",*p++);
		fprintf(fout,"\n");
	}
	fclose(fout);
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
