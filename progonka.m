clc
l=-1;
N=100;
r=10;
h=(r-l)/N
x=l:h:r;
%difur: a(x)y''+b(x)y'+c(x)y=F(x)
a=ones(size(x));
b=zeros(size(x));
c=ones(size(x));
F=zeros(size(x));
%gran usl: 
%alp__l*y-bet__l*y'=mu__l
%alp__r*y+bet__r*y'=mu__r
alp__l=0;
bet__l=-1;
mu__l=1;
alp__r=0;
bet__r=1;
mu__r=0;

A=a./h./h+b./h./2
B=a./h./h-b./h./2
C=2*a./h./h-c
F
alp_l=   alp__l/(alp__l+bet__l*h)
bet_l=-h*mu__l/(alp__l+bet__l*h)
alp_r=   alp__r/(alp__r+bet__r*h)
bet_r= h*mu__r/(alp__r+bet__r*h)

%todo
%написать проверки устойчивости прогонки

fout=fopen('progonka/input.txt','w');
fprintf(fout,'%d %g %g %g %g\n',length(x)-1,alp_l,bet_l,alp_r,bet_r);
for q=A
    fprintf(fout,'%g ',q);
end;
  fprintf(fout,'\n');
for q=B
    fprintf(fout,'%g ',q);
end;
  fprintf(fout,'\n');
for q=C
    fprintf(fout,'%g ',q);
end;
  fprintf(fout,'\n');
for q=F
    fprintf(fout,'%g ',q);
end;
  fprintf(fout,'\n');
fclose(fout);

err=system('progonka\Debug\progonka.exe progonka\');
if(err)
    error('error in call progonka');
end;

fin=fopen('progonka/output.txt');
y=fscanf(fin,'%f',length(x));
plot(x,y);

system('del progonka\input.txt progonka\output.txt');

