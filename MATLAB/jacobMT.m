function [J_appres]=jacobMT(par,init_model,thicknesses,frequency)



     J_appres=zeros(length(frequency),length(init_model));
     J_phase=zeros(length(frequency),length(init_model));
for i=1:length(frequency)
    [app_rho_m,phase_m]=JIforwardMT(exp(init_model),thicknesses,frequency(i));
    for j=1:length(init_model)
        mdm=init_model;
        dm=((init_model(j))*par);
        mdm(j)=(init_model(j))+dm;
    [app_rho_mdm,phase_mdm]=JIforwardMT(exp(mdm),thicknesses,frequency(i));
    J_appres(i,j)=((log(app_rho_mdm)-log(app_rho_m))./(dm));
    J_phase(i,j)=((log(phase_mdm)-log(phase_m))./(dm));
%     df/dm=(f(mdm)-f(m))/dm
    end
end
%dilakukan di skala log

end
