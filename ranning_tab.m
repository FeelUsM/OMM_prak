clc
x0=0;
xr=-1;
t0=0;
tr=5;
Nx=100;
Nt=100;
hx=(xr-x0)/Nx;
tau=(tr-t0)/Nt;
x=x0:hx:xr;
t=t0:tau:tr;
[X T]=meshgrid(x,t);

%todo
%рисовать график не по всем точкам

epsilon=0.001;
maxiter=10000000;

u_x0=1-x;
u_0t=(2-t.^2)./(4*t+2);

fin=fopen('ranning_tab\input.txt','w');
fprintf(fin,'%g %g %d \n%g %g %d \n%g %d\n',x0,hx,Nx,t0,tau,Nt,epsilon,maxiter);
for q=u_0t
    fprintf(fin,'%g ',q);
end;
fprintf(fin,'\n');
for q=u_x0
    fprintf(fin,'%g ',q);
end;
fclose(fin);

err=system('ranning_tab\Debug\ranning_tab.exe ranning_tab/input.txt ranning_tab/output.txt');
if(err)
    error('error in call ranning_tab');
end;

u=[];
fout=fopen('ranning_tab\output.txt');
sriters=fscanf(fout,'%f',1)
for i=0:Nt
    %disp(i)
    q=fscanf(fout,'%f ',Nx+1);
    u=[u; q.'];
end;
fclose(fout);

%window=figure;
surf(X,T,u);
%print(window,'-dpng','ranning_tab');
%close(window);


