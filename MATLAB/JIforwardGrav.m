function [g]=JIforwardGrav(model,dx,coor_m_x,coor_obs_x,dz,coor_m_z,coor_obs_z)

[m,n]=size(model);

g.const=6.67e-11; %G


for i=1:length(coor_obs_x);
    kk=1;
    clc
    progress=i*100/length(coor_obs_x)
        for ii=1:length(coor_m_x);
            for jj=1:length(coor_m_z);
                A=coor_obs_x(i)-coor_m_x(ii)+dx/2;
                B=coor_obs_x(i)-(coor_m_x(ii)+dx/2);
                C=coor_obs_z(i)-coor_m_z(jj)+dz/2;
                D=coor_obs_z(i)-(coor_m_z(jj)+dz/2);
%                 plot(A,zi(jj),'-bo');
%                 hold on
%                 plot(B,zi(jj),'-bo');
%                 plot(C,xi(ii),'-bo');
%                 plot(D,xi(ii),'-bo');
                g.kernel(i,kk)=g.const*(A*log((A^2+D^2)/(A^2+C^2))-B*log((B^2+D^2)/(B^2+C^2))+2*D*(atan(A/D)-atan(B/D))-2*C*(atan(A/C)-atan(B/C)));
            kk=kk+1;
            end
        end
end
[m,n]=size(model);
g.rho=reshape(model,[m*n,1]);
g.obs=g.kernel*g.rho;

g.obs=g.obs*1e8;
disp('==========================================================================')
disp('==========================================================================')
disp('====                                                                  ====')
disp('====                           Forward Gravity                        ====')
disp('====                    Based  on Last & Kubik 1983                   ====')
disp('====                               By                                 ====')
disp('====                        Teknik Geofisika ITS                      ====') 
disp('====                               ---                                ====')
disp('====                        M.Irsyad Hibatullah                       ====')
disp('====                         Jeremy Gohitua M.                        ====')
disp('====                           Nuha Malihati                          ====')
disp('====                         Firman Syaifuddin                        ====')
disp('====                          Dwa Desa Warnana                        ====')
disp('====                          Juan Pandu G.N.A                        ====')
disp('====                                                                  ====')
disp('==========================================================================')
disp('==========================================================================')

end




%             r1=sqrt(((coor_m_z(k)-dz/2)^2)+((coor_obs_x(i)-coor_m_x(j)+dx/2)^2));
%             r2=sqrt(((coor_m_z(k)+dz/2)^2)+((coor_obs_x(i)-coor_m_x(j)+dx/2)^2));
%             r3=sqrt(((coor_m_z(k)-dz/2)^2)+((coor_obs_x(i)-coor_m_x(j)-dx/2)^2));
%             r4=sqrt(((coor_m_z(k)+dz/2)^2)+((coor_obs_x(i)-coor_m_x(j)-dx/2)^2));
%             
%             theta1=atan((coor_obs_x(i)-coor_m_x(j)-dx/2)/(coor_m_z(k)-dz/2));
%             theta2=atan((coor_obs_x(i)-coor_m_x(j)-dx/2)/(coor_m_z(k)-dz/2));
%             theta3=atan((coor_obs_x(i)-coor_m_x(j)-dx/2)/(coor_m_z(k)-dz/2));
%             theta4=atan((coor_obs_x(i)-coor_m_x(j)-dx/2)/(coor_m_z(k)-dz/2));
%    g.kernel(i,kk)=2*g.const*((coor_obs_x(i)-coor_m_x(j)+dx/2)*log(r2*r3/r1*r4)+dx*log(r3/r4)-(coor_m_z(k)+dz/2)*(theta4-theta2)+(coor_m_z(k)-dz/2)*(theta3-theta1));
%       kk=kk+1;

%                 r=((coor_m_x(j)-coor_obs_x(i))^2)+((coor_m_z(k)-coor_obs_z(i))^2);
% %         (2*g.const*dx*dz*((coor_m_z(j)-coor_obs_z(i))))/r;
%         g.kernel(i,kk)=(2*g.const*(coor_m_z(k)-coor_obs_z(i))*dx*dz)/r;
