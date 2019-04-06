
clear
%load data
Grav=load('20_point_Gravity_Inv_Result.mat');
MT=load('20_point_15D_MT_Inv_Result_2.mat');


data_MT=MT.data(:);
model_MT=MT.model(:);

data_Grav=Grav.data(:);
model_Grav=Grav.model(:);
clear Grav MT

cell_num_x=40;
cell_num_z=20;

min_x=0;
max_x=20000;

min_z=0;
max_z=10000;

[x,z]=JImeshcreate(min_x,max_x,cell_num_x,min_z,max_z,cell_num_z);


clear cell_num_x cell_num_z min_x min_z max_z max_x


%% plot data observasi
plot_param.MT_limitphase=[min(data_MT.obs_phase(:)) max(data_MT.obs_phase(:))];
plot_param.MT_limitappres=log([min(data_MT.obs_apparent_res(:)) max(data_MT.obs_apparent_res(:))]);

plot_param.Grav_model=[min(model_Grav.True_Dens(:)) max(model_Grav.True_Dens(:))];
plot_param.MT_model=log([min(model_MT.True_Res(:)) max(model_MT.True_Res(:))]);


% JIplot_data_mt(x,data_MT,'obs','figure(1)','subplot(2,3,2)','subplot(2,3,3)',plot_param.MT_limitphase,plot_param.MT_limitappres)
% JI_plot_Gobs(data_Grav,'obs','figure(1)','subplot(2,3,1)')


%% pendefinisian parameter inversi
  inv_param_Grav.eps=0.2;
  inv_param_Grav.eps_min=0.0002;
  inv_param_Grav.eps_max=0.0002;
  
  inv_param_Grav.RMSE=100;
  inv_param_Grav.min_RMSE=0.1;
  
  inv_param_Grav.max_Iter=50;
  
  inv_param_Grav.optim='pcg';
  
  inv_param_Grav.flat_degree=0.5;
  inv_param_Grav.smoothness_order=1;
  
  inv_param_Grav.smoothness_order_reference=4;
  inv_param_Grav.tau=0.5; %reference weight
  
  inv_param_Grav.compactness='off';
  inv_param_Grav.beta=0.3;
  


  inv_param_MT.eps=0.2;
  inv_param_MT.eps_min=0.0000000002;
  inv_param_MT.eps_max=0.8;

  inv_param_MT.RMSE=100;
  inv_param_MT.min_RMSE=0.001;
  
  inv_param_MT.max_Iter=20;
  inv_param_MT.optim='pcg';
  
  inv_param_MT.flat_degree=0.5;
  inv_param_MT.smoothness_order=2;
  inv_param_MT.smoothness_order_reference=4;
  inv_param_MT.tau=0.5; %reference weight
  
  inv_param_MT.FDpar=1/100;
  
  inv_param_MT.min_m=log(50);
  inv_param_MT.max_m=log(800);
  
  inv_param_Grav.min_m=500;
  inv_param_Grav.max_m=3500;
 

a_dens=0.31;
a_res=240;
b=0.25;
c=0.24;
  
%% pendefinisian Initial Model

  model_MT.initial   = JI_init_model_increment(5,1000,data_MT.point_number,z.cell_quantity); %inisial model resistivitas
  model_MT.initial   = log(model_MT.initial);
  model_Grav.initial = JI_init_model_increment(1000,3000,x.cell_quantity,z.cell_quantity); %inisial model densitas

% model_MT.initial   =abs(imresize(model_MT.initial,[z.cell_quantity data_MT.point_number]));
% model_Grav.initial =abs(imresize(model_Grav.initial,[z.cell_quantity x.cell_quantity]));


%% Pre-Inversi
  %Grav
  inv_param_Grav.G = JI_kernel_grav(x,z,data_Grav);
  inv_param_Grav.Wm = JI_wm_2D(size(model_Grav.initial),inv_param_Grav,0);
  inv_param_Grav.W = JI_depth_weight_grav(x,z,inv_param_Grav,model_Grav);
  inv_param_Grav.Wref1 = JI_wm_2D([z.cell_quantity x.cell_quantity],inv_param_Grav,1);
  
  %MT
  inv_param_MT.Wm = JI_wm_2D(size(model_MT.initial),inv_param_MT,0);
  inv_param_MT.Wref1   = JI_wm_2D([z.cell_quantity data_MT.point_number],inv_param_MT,1);
  
  tau=tau_matrix(model_Grav,model_MT,inv_param_Grav,inv_param_MT);

  
Wm=blkdiag(inv_param_Grav.Wm,inv_param_MT.Wm);
Wref1=blkdiag(inv_param_Grav.Wref1,inv_param_MT.Wref1);

subplot(2,1,1)
pcolor(extend(model_Grav.initial));axis ij;colorbar ;caxis(plot_param.Grav_model)


subplot(2,1,2)
pcolor(extend(model_MT.initial));axis ij;colorbar ;caxis(plot_param.MT_model)


%% Run Inversi
Iter=1;
while Iter < inv_param_MT.max_Iter

    J = JI_jacobian_MT_15D_2(z, model_MT, data_MT, inv_param_MT);
    G = blkdiag(inv_param_Grav.G*inv_param_Grav.W,J);

    model_Grav.ref_MT = JI_res2dens(z.coor_m,a_dens,a_res,b,c,exp(JI_interp_MT15D (x,z,data_MT,model_MT)));
    model_MT.ref_Grav = log(JI_get_point_model(x,data_MT,(JI_dens2res(z.coor_m,a_dens,a_res,b,c,(model_Grav.initial)))));

    Mref1=[(model_Grav.initial(:)-model_Grav.ref_MT(:));(model_MT.initial(:)-model_MT.ref_Grav(:))];

    [data_MT.cal_apparent_res , data_MT.cal_phase] = JIforwardMT_15D(z , data_MT , model_MT);

    dg = data_Grav.Gobs(:)-(inv_param_Grav.G*model_Grav.initial(:));

    inv_param_Grav.RMSE = sqrt((dg'*dg)./length(dg));
    [df , inv_param_MT] = JI_df_MT_2(data_MT, inv_param_MT);

    disp(['Iter:  ' num2str(Iter)])

    disp(['RMS Grav: ' num2str(inv_param_Grav.RMSE) ' RMS MT: ' num2str(inv_param_MT.RMSE)])

    dd=[dg(:);df(:)];

    eps=eps_matrix(model_Grav,model_MT,inv_param_Grav,inv_param_MT);
% N=(G'*G)+(Wm.*eps)+(Wref1);
% n=(G'*dd)+Wref1*Mref1;

N=(G'*G)+(Wm.*eps);
n=(G'*dd);



dm=bicgstab(N,n);
clc
disp(['Iter:  ' num2str(Iter)])
disp(['RMS Grav: ' num2str(inv_param_Grav.RMSE) ' RMS MT: ' num2str(inv_param_MT.RMSE)])
dm_grav = dm(1:numel(model_Grav.initial));
dm_mt   =dm(numel(model_Grav.initial)+1:numel(model_Grav.initial)+numel(model_MT.initial));

model_Grav.initial=model_Grav.initial+reshape(dm_grav,size(model_Grav.initial));
model_MT.initial=model_MT.initial+reshape(dm_mt,size(model_MT.initial));

subplot(2,1,1)
pcolor(extend(model_Grav.initial));axis ij;colorbar ;caxis(plot_param.Grav_model)

subplot(2,1,2)
pcolor(extend(model_MT.initial));axis ij;colorbar ;caxis(plot_param.MT_model)

Iter=Iter+1;
pause(0.01)
end








%% fungsi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Joint MT-Grav %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eps=eps_matrix(model_Grav,model_MT,inv_param_Grav,inv_param_MT)
epsgrav=ones(numel(model_Grav.initial)).*inv_param_Grav.eps;
epsmt  = ones(numel(model_MT.initial)).*inv_param_MT.eps;
eps=blkdiag(epsgrav,epsmt);
end

function tau=tau_matrix(model_Grav,model_MT,inv_param_Grav,inv_param_MT)
taugrav=ones(numel(model_Grav.initial)).*inv_param_Grav.tau;
taumt  = ones(numel(model_MT.initial)).*inv_param_MT.tau;
tau=blkdiag(taugrav,taumt);
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Grav %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Grav engine
  function [model,data,inv_param_grav] = JI_inv_grav(x,z,model,data,inv_param_grav)
  G=inv_param_grav.G;
  H=inv_param_grav.H;
  W=inv_param_grav.W;
  data.Gcal=G*model.initial(:);
  dg=data.Gobs-data.Gcal;
  inv_param_grav.RMSE=sqrt((dg'*dg)/length(dg));
  if strcmp(inv_param_grav.optim,'pinv')
  dm=(pinv((((G*W)'*(G*W))+W.*inv_param_grav.eps)))*((G*W)'*dg);
  else
     dm=eval(strcat(string(inv_param_grav.optim),"((((G*W)'*(G*W))+W.*inv_param_grav.eps),(G*W)'*dg);"));clc
  end
  model.initial=model.initial+reshape(dm,[z.cell_quantity x.cell_quantity]);
  
  end
  function [W]                         = JI_depth_weight_grav(x,z,inv_param_grav,model)
  P=eye(x.cell_quantity*z.cell_quantity).*0.001; % if geology is known, it became 0.001;
  V=P; % first iteration V= identitiy matrix
  Q=P; %pre allocating
for j=1:(x.cell_quantity*z.cell_quantity)
    [row,~]=ind2sub([z.cell_quantity x.cell_quantity],j);    
    Q(j,j)=1/((z.coor_m(row)+10e-10)^inv_param_grav.beta);
    if strcmp(inv_param_grav.compactness,'on')
    V(j,j)=(1./((model.initial(j)^2)+inv_param_grav.eps));
    else
    V(j,j)=1;
    end
end
W=(P)*Q*V;
W=inv(W);
  end
  function JIplot_model_grav(x,z,model,data,datatype,Iter,inv_param_grav,fignum,subnum)
 
RMSE=inv_param_grav.RMSE;

eval(fignum)
eval(subnum)

if strcmp('true',datatype)
    pcolor(x.coor_plot,z.coor_plot,extend(imresize(model.True_Dens,[z.cell_quantity x.cell_quantity])));
    
    axis ij; c=colorbar; caxis([min(min(model.True_Dens)) max(max(model.True_Dens))]) ; c.Label.String = 'Density';
    
else
    pcolor(x.coor_plot,z.coor_plot,extend(model.initial));
    axis ij; c=colorbar; caxis([min(min(model.True_Dens)) max(max(model.True_Dens))]) ; c.Label.String = 'Density';
end


hold on
plot(data.grav_coor_x,data.grav_coor_z,'o','MarkerSize',5,'MarkerFaceColor','r','MarkerEdgeColor','none')
xlabel('Distance (m)')
ylabel('Depth (m)')
hold off



% if strcmp('true',datatype)
%     title('True Model')
% else
%     title(['2D Gravity Inv Result | Iteration : ' num2str(Iter) ' | RMS Error : ' num2str(RMSE(Iter))])
% end

  end
  function JI_plot_Gobs(data,datatype,fignum,subnum)
 eval(fignum)
 eval(subnum)
 if strcmp(datatype,'cal')
    plot(data.grav_coor_x,data.Gcal,'b-');
 elseif strcmp(datatype,'obs')
    plot(data.grav_coor_x,data.Gobs,'b.','MarkerSize',10);
 elseif strcmp(datatype,'cal & obs') || strcmp(datatype,'obs & cal')
    plot(data.grav_coor_x,data.Gobs,'b.','MarkerSize',10);
        hold on
    plot(data.grav_coor_x,data.Gcal,'b-');
    hold off
 else
     disp('data type error')
     pause
 end
 ylabel('miliGal')
 xlabel('distance(m)')
  if strcmp('obs',datatype)
    title('G Obs')
 elseif strcmp('cal',datatype)
    title('G cal')
 else
     title('G obs & cal')
 end
  end
 
 %% MT Engine
  function [model,data , inv_param]  = JI_inv_MT15D(z,model, data, inv_param)
     J  =  JI_jacobian_MT_15D_2(z,model,data,inv_param);
        [data.cal_apparent_res,data.cal_phase] =  JIforwardMT_15D(z,data,model);
        [df,inv_param]                         =  JI_df_MT_2(data,inv_param);
        if inv_param.RMSE<inv_param.min_RMSE
            disp('iteration has stopped because minimum RMSE has reached')
        else
        [U,s,~]  =  JI_csvd(J);

        inv_param.eps = l_curve(U,s,df);
            if inv_param.eps_max<inv_param.eps
               inv_param.eps=inv_param.eps_max;
            elseif inv_param.eps_min>inv_param.eps
               inv_param.eps=inv_param.eps_min;
            end
        inv_com="(((J'*J)+(inv_param.H.*inv_param.eps)),J'*df);";
        inv_com=strcat(inv_param.optim,inv_com);
        dm=eval(string(inv_com)); clc
        model.initial=model.initial+reshape(dm,size(model.initial));
        end
    end
  function [df , inv_param]          = JI_df_MT_2(data, inv_param)
    appres_idx = 1:2:(data.point_number*2);
    phase_idx  = 2:2:(data.point_number*2);
    cal(:,appres_idx)=data.cal_apparent_res;
    cal(:,phase_idx)=data.cal_phase;

    obs(:,appres_idx)=data.obs_apparent_res;
    obs(:,phase_idx)=data.obs_phase;

    df=log(obs(:))-log(cal(:));
    inv_param.RMSE=sqrt((df'*df)/length(df));
    end
  function [J]                       = JI_jacobian_MT_15D_2(z, model, data, inv_param)
    thickness=ones(z.cell_quantity-1,1).*z.dz;
    %pre allocating
    J_appres=zeros(length(data.frequency(:,1)),length(model.initial(:,1)),data.point_number);
    J_phase=J_appres;

    for k=1:data.point_number
    [J_appres(:,:,k),J_phase(:,:,k)]=devel_JI_jacobianMT_02(inv_param.FDpar,model.initial(:,k),thickness,data.frequency(:,k));

    if k==1
    comJ=strcat("([J_appres(:,:,",num2str(k),");J_phase(:,:,",num2str(k),")])");
    else
        comJ=strcat(comJ,",([J_appres(:,:,",num2str(k),");J_phase(:,:,",num2str(k),")])");
    end

    end
    J=eval(strcat("blkdiag(",comJ,");"));

    end
  function [d_cal_app , d_cal_phase] = JIforwardMT_15D(z , data , model)
    thickness=ones(z.cell_quantity-1,1).*z.dz;
    %pre allocating
    d_cal_app=zeros(size(data.frequency,1),data.point_number);
    d_cal_phase=d_cal_app;

    for j=1:data.point_number
        for i=1:size(data.frequency,1)
            [d_cal_app(i,j),d_cal_phase(i,j)]=JIforwardMT(exp(model.initial(:,j)),thickness,data.frequency(i,j));
        end
    end

    end

%% Grav plot data & model


%% Grav MISC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=JI_get_point_model(x,data,model_initial)
model_col_ind = dsearchn(x.coor_m',data.coor_mt_x);
out=model_initial(:,model_col_ind);
end

 %% MT Plot Function
function JI_plot_1D_model(model,point,z,Iter,fig_num,submodel,limit)
%plot 1D curve model sebagai tools untuk analisa hasil inversi
%sementara kecurigaan ada pada 
%   -augmentasi matriksnya masih salah, mungkin bisa dicoba cari referensi
%       lain
%   - memang solusi 1Dnya kalau di tumpuk seperti itu, cara nge ceknya
%       yha... di plot satu-satu

eval(fig_num)
eval(submodel)
temp=(linspace(z.min,z.max,size(model.True_Res,1)));
temp=temp(2)-temp(1);

stairs([0;model.True_Res(:,point)],...
    cumsum([0;ones(size(model.True_Res,1),1).*temp]),':k','LineWidth',1.5) ;axis ij ;xlim(limit)
hold on
stairs([0;exp(model.inv_result(:,point,Iter))],...
    cumsum([0; diff(z.coor_plot)']),'b','LineWidth',1.5) ;axis ij ;xlim(limit)
hold off

xlabel('Resistivity(rho)')
% ylabel('Depth(m)')
 end
function JIplot_model_mt(x,z,model,data,datatype,inv_param,Iter,fig_num,submodel,limit)

RMSE=inv_param.RMSE;

if strcmp('true',datatype)
    plotmodel=imresize(log(model.True_Res),[z.cell_quantity x.cell_quantity]);
else
    plotmodel=JI_interp_MT15D (x,z,data,model);
end

eval(fig_num)
eval(submodel)
pcolor(x.coor_plot,z.coor_plot,extend((plotmodel)));axis ij; colorbar ; caxis(limit);

xlabel('Distance (m)')
ylabel('Depth (m)')
hold on
plot(data.coor_mt_x,data.coor_mt_z,'v','MarkerSize',7,'MarkerFaceColor','k','MarkerEdgeColor','none')
hold off
% legend('Log ( Resistivity (Ohm/m) )','MT point','Location','southoutside')

if strcmp('true',datatype)
    title('True Model')
else
    title(['1.5D Inv Result | Iteration : ' num2str(Iter) ' | RMS Error : ' num2str(RMSE)])
end

end
function JIplot_data_mt(x,data,datatype,fig_num,subphase,subappres,limitphase,limitappres)
model_col_ind = dsearchn(x.coor_m',data.coor_mt_x);

plotdata.apparent_res=nan(size(data.frequency,1),x.cell_quantity);
plotdata.phase=plotdata.apparent_res;
% % create temporary plot model
if strcmp('obs',datatype)
    plotdata.apparent_res(:,model_col_ind)=log(data.obs_apparent_res);
    plotdata.phase(:,model_col_ind)=(data.obs_phase);
else % 'cal'
    plotdata.apparent_res(:,model_col_ind)=log(data.cal_apparent_res);
    plotdata.phase(:,model_col_ind)=(data.cal_phase);
end

[plotdata.mt_point_x1,plotdata.mt_point_z1]=meshgrid(data.coor_mt_x,min(log([data.frequency(1,1);data.frequency(:,1)])));
[plotdata.mt_point_x2,plotdata.mt_point_z2]=meshgrid(data.coor_mt_x,max(log([data.frequency(1,1);data.frequency(:,1)])));
% % plot Phase
eval(fig_num)
eval(subphase)
pcolor(x.coor_plot,log([data.frequency(1,1);data.frequency(:,1)]),extend(inpaintn(plotdata.phase))); colorbar ;caxis(limitphase)
xlabel('distance(m)')
ylabel('Log(  Frequency(Hz)  )')
hold on
plot([plotdata.mt_point_x1],plotdata.mt_point_z2,'v','MarkerSize',7,'MarkerFaceColor','k','MarkerEdgeColor','none')
hold off
% legend('Phase(Degree)','MT point','Location','southoutside')

if strcmp('obs',datatype)
    title('Obs.Phase')
else
    title('Cal.Phase')
end
 % % plot apparent Resistivity
eval(subappres)
pcolor(x.coor_plot,log([data.frequency(1,1);data.frequency(:,1)]),extend(inpaintn(plotdata.apparent_res))); colorbar ;caxis(limitappres)
xlabel('distance(m)')
ylabel('Log(  Frequency(Hz)  )')
hold on
plot([plotdata.mt_point_x1],plotdata.mt_point_z2,'v','MarkerSize',7,'MarkerFaceColor','k','MarkerEdgeColor','none')
% legend('Apparent Resistivity(Log(Rho))','MT point','Location','southoutside')

if strcmp('obs',datatype)
    title('Obs.Apparent Resistivity')
else
    title('Cal.Apparent Resistivity')
end

end

function plotmodel=JI_interp_MT15D (x,z,data,model)
%interpolasi NaN value pada MT 1.5D
%mungkin bisa di cek kecepatannya apabila dibandingkan dengan fungsi
%interp/griddata

model_col_ind = dsearchn(x.coor_m',data.coor_mt_x);
plotmodel=nan(z.cell_quantity,x.cell_quantity);

    plotmodel(:,model_col_ind)=model.initial;
    plotmodel=inpaintn(plotmodel);
end

% function JI_MT_match_point
%%  MT MISC
function A=extend(A)
A=[A zeros(size(A,1),1)
    zeros(1,size(A,2)+1)];
end
function model_initial=JI_init_model_increment(min_mod,max_mod,x_num,z_num)
    disp('create initial model')
    mod=linspace(min_mod,max_mod,z_num); 
    model_initial=zeros(z_num,x_num);
         for i=1:z_num
          model_initial(i,1:x_num)=mod(i);
         end
end

%% catatan konversi

% a_dens=0.31;
% a_res=240;
% b=0.25;
% c=0.24;
% 
% JI_res2dens(z_depth,a_dens,a_res,b,c,res)
% JI_dens2res(z_depth,a_dens,a_res,b,c,dens)
% JI_dens2vp(a_dens,b,dens)
% JI_vp2dens(a_dens,b,vp)
% JI_vp2res(a_res,c,Vp,z)
% JI_res2vp(a_res,c,Res,z)



