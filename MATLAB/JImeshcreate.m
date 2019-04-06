function [x,z]=JImeshcreate(min_x,max_x,cell_quantity_x,min_z,max_z,cell_quantity_z)

x.min=min_x;
x.max=max_x;
x.cell_quantity=cell_quantity_x;

z.min=min_z;
z.max=max_z;
z.cell_quantity=cell_quantity_z;


%persiapan plot,
x.coor_plot=linspace(x.min,x.max,x.cell_quantity+1);
z.coor_plot=linspace(z.min,z.max,z.cell_quantity+1);

x.dx=x.coor_plot(2)-x.coor_plot(1);
z.dz=z.coor_plot(2)-z.coor_plot(1);
%persiapan plot value dari setiap cell


%ketika kita ingin plot 5x5 cell, maka harus dibuat
%boundary sebanyak x=1:6 dan z=1:6 dan value modelnya harus 6x6, dimana
%(:,last)=0 dan (last,:)=0

%nilai koordinat dari model, dimana nilai koordinatnya adalah titik tengah
%dari cell alias dx/2 dan dz/2 untuk cell pertama.
x.coor_m=x.min+x.dx/2:x.dx:(x.max-x.dx/2);
z.coor_m=z.min+z.dz/2:z.dz:(z.max-z.dz/2);
