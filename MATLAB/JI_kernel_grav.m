function a=JI_kernel_grav(x,z,data)

x_m=x.coor_m;
z_m=z.coor_m;
dx=x.dx;
dz=z.dz;
xobs=data.grav_coor_x;
zobs=data.grav_coor_z;


Gg=6.67e-11; %G
a=zeros(numel(xobs),length(x_m)*length(z_m));
    for i=1:length(xobs)
    k=1;
        for ii=1:length(x_m)
            for jj=1:length(z_m)
                A=xobs(i)-x_m(ii)+dx/2;
                B=xobs(i)-(x_m(ii)+dx/2);
                C=zobs(i)-z_m(jj)+dz/2;
                D=zobs(i)-(z_m(jj)+dz/2);

                a(i,k)=Gg*(A*log((A^2+D^2)/(A^2+C^2)+10e-10)...
                -B*log((B^2+D^2)/((B^2+C^2)+10e-10))...
                +2*D*(atan(A/(D+10e-10))-atan(B/(D+10e-10)))...
                -2*C*(atan(A/(C+10e-10))-atan(B/(C+10e-10))));
%             if k==1771 && i==length(xobs)
%                 
%                y=asdfghj; 
%             end
            k=k+1;
            end
        end
    end
    
a=a*1e5;
end
