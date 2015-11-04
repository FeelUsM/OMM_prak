clc
x0=0;
xr=1;
y0=0;
yr=1;
Nx=100;
Ny=100;
hx=(xr-x0)/Nx;
hy=(yr-y0)/Ny;
x=x0:hx:xr;
y=y0:hy:yr;
[X Y]=meshgrid(x,y);

%todo
% в print/scan заменить %f на %g
% научиться писать часла в матлабовские строки

if 1 %начальные условия
indata=zeros(size(X));%функция начальных условий
fin=fopen('alternating_direction/input_layer.txt','w');
for iy=1:Ny+1
    for i=indata(iy,:)
        fprintf(fin,'%g ',i);
    end;
    fprintf(fin,'\n');
end;
fclose(fin);
end;

T=0;
dT=0.01;
Nt=100;
tau=dT/Nt;

window=figure;
names=[
    'alternating_direction_1'
    'alternating_direction_2'
    'alternating_direction_3'
    'alternating_direction_4'
    'alternating_direction_5'
    'alternating_direction_6'
    'alternating_direction_7'
    ];
for name=names.' %несколько промежутков по времени, после которых рисуем графики
name=name.';
    
%параметры всего задания и гран. условия
fgu=fopen('alternating_direction/gran_usl.txt','w');
fprintf(fgu,'%g \n%g %g %g %g %d %d\n %g %g %d\n',9., x0, y0, hx, hy, Nx, Ny, T, tau, Nt);
fprintf(fgu,' %g %g %g %g %g %g %g %g\n', 0, -1, 0, 1, 0, -1, 1, 0);
fclose(fgu);

disp('from')
disp(T)
disp('to')
T=T+dT;
disp(T)

% ======= вызов ========
err=system('alternating_direction\Release\alternating_direction.exe alternating_direction/gran_usl.txt alternating_direction/input_layer.txt alternating_direction/output_layer.txt');
if(err)
    error('error in call alternating_direction');
end;
%очистка и переименование файлов
system('del alternating_direction\input_layer.txt');
system('ren alternating_direction\output_layer.txt input_layer.txt');

u=[]; % считываем результат
fout=fopen('alternating_direction\input_layer.txt');
%fout=fopen('alternating_direction\output_layer.txt')
for i=0:Ny
    %disp(i)
    q=fscanf(fout,'%f ',Nx+1);
    u=[u; q.'];
end;
fclose(fout);

surf(X,Y,u);
print(window,'-dpng',name);

end; %for name=names

close(window);
%system('del alternating_direction\input_layer.txt');
system('del alternating_direction\gran_usl.txt');

%------------------------------------------------------
%неоднородность моего ур-я в t=0
window=figure;
surf(X,Y,X.*Y.*(1-X));
print(window,'-dpng','alternating_direction_neodnorodnost');
close(window);



