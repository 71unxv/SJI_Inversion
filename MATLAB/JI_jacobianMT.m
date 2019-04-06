% Function Jacobian %
function[A] = JI_jacobianMT(frequency_vector,m0_res,thicknesses,data_obs,tipe_inversi);
if tipe_inversi=='1D';
    roa=data_obs;
    par = 0.01;
    lr=length(m0_res);
    lt=length(thicknesses);
        for i=1:length(frequency_vector);
            [roa1]=JIforwardMT(m0_res, thicknesses,frequency_vector(i));
        end 
    clear i
     r = m0_res;
    lr=length(m0_res);
    r2=r;
        for i2 = lr:-1:1;
         r2(i2) = (r(i2)*par)+r(i2);
               for ii = 1:length(frequency_vector);
                  s = frequency_vector(ii);
                  [g] = JIforwardMT(r2,thicknesses,s);
                 roa2(ii,:) = g;
                end
           A1(:,i2) = [(roa2-roa1)/(r(i2)*par)]*r(i2)./roa';
           r2 = r;
        end
    A = [A1];
end
% end
% 
% if tipe_inversi=='1.5D';
% 
%     roa=data_obs;
%     par = 0.1;
%     [lr,lrr]=size(m0_res);
% 
%     for i=1:length(frequency_vector);
%         [roa1]=JIforwardMT(m0_res, thicknesses,frequency_vector(i));
%     end 
%         for j=1:lrr;
%         if j==1;
%                clear i
%         r = m0_res;
%         lr=length(m0_res);
%         r2=r;
%             for i2 = 1:lr;
%              r2(i2) = (r(i2)*par)+r(i2);
%                    for ii = 1:length(frequency_vector);
%                       s = frequency_vector(ii);
%                       [g] = JIforwardMT(r2,thicknesses,s);
%                      roa2(ii,:) = g;
%                     end
%                A1(:,i2) = [(roa2-roa1)/(r(i2)*par)]*r(i2)./roa';
%                r2 = r;
%             end
%         A = [A1];
%         else
%         clear i
%         r = m0_res;
%         lr=length(m0_res);
%         r2=r;
%             for i2 = 1:lr;
%              r2(i2) = (r(i2)*par)+r(i2);
%                    for ii = 1:length(frequency_vector);
%                       s = frequency_vector(ii);
%                       [g] = JIforwardMT(r2,thicknesses,s);
%                      roa2(ii,:) = g;
%                     end
%                A1(:,i2) = [(roa2-roa1)/(r(i2)*par)]*r(i2)./roa';
%                r2 = r;
%             end
%         A = [A A1];
%         end


        end
% end
% 
% end