%JIforwardMT(resistivities, thicknesses,frequency)

function [apparentResistivity, phase] = JIforwardMT(resistivities, thicknesses,frequency)

mu = 4*pi*1E-7; %Magnetic Permeability (H/m)  
w = 2 * pi * frequency; %Angular Frequency (Radians);
nn=length(resistivities); %Number of Layers

impedances = zeros(nn,1);
%www.digitalearthlab.com
%Layering in this format
% Layer      j
% Layer 1    1
% Layer 2    2
% Layer 3    3
% Layer 4    4
% Basement   5

% Steps for modelling (for each geoelectric model and frequency)
% 1. Compute basement impedance Zn using sqrt((w * mu * resistivity))
% 2. Iterate from bottom layer to top(not the basement)
    % 2.1. Calculate induction parameters
    % 2.2. Calculate Exponential factor from intrinsic impedance
    % 2.3 Calculate reflection coeficient using current layer
    %          intrinsic impedance and the below layer impedance
    
% 3. Compute apparent resistivity from top layer impedance
        %   apparent resistivity = (Zn^2)/(mu * w)

%Symbols
% Zn - Basement Impedance
% Zi - Layer Impedance
% wi - Intrinsic Impedance
% di - Induction parameter
% ei - Exponential Factor
% ri - Reflection coeficient
% re - Earth R.C.
        
%Step 1 : Calculate basement impedance  
Zn = sqrt(sqrt(-1)*w*mu*resistivities(nn)); 
impedances(nn) = Zn; 

%Iterate through layers starting from layer j=n-1 (i.e. the layer above the basement)        
for j = nn-1:-1:1
%     resistivities(j) = resistivities(j);
%     thicknesses = thicknesses(j);
                
    % 3. Compute apparent resistivity from top layer impedance
    %Step 2. Iterate from bottom layer to top(not the basement) 
    % Step 2.1 Calculate the intrinsic impedance of current layer
    dj = sqrt(sqrt(-1)* (w * mu * (1/resistivities(j))));% dj - Induction parameter
    wj = dj * resistivities(j); %wj - Intrinsic Impedance
    
        
    % Step 2.2 Calculate Exponential factor from intrinsic impedance
    ej = exp(-2*thicknesses(j)*dj); %ej - Exponential Factor              

    % Step 2.3 Calculate reflection coeficient using current layer
    %          intrinsic impedance and the below layer impedance
    belowImpedance = impedances(j + 1);
    rj = (wj - belowImpedance)/(wj + belowImpedance); % rj - Reflection coeficient
    re = rj*ej; % re - Earth R.C.
    Zj = wj * ((1 - re)/(1 + re)); % Zj - Layer Impedance
    
    impedances(j) = Zj;               
end
% Step 3. Compute apparent resistivity from top layer impedance
Z = impedances(1);
absZ = abs(Z); 

apparentResistivity = (absZ * absZ)/(mu * w);
phase = atan2(imag(Z),real(Z));
% disp('==========================================================================')
% disp('==========================================================================')
% disp('====                                                                  ====')
% disp('====                      Forward Magnetotelluric                     ====')
% disp('====                    Based  on .................                   ====')
% disp('====                               By                                 ====')
% disp('====                        Teknik Geofisika ITS                      ====') 
% disp('====                               ---                                ====')
% disp('====                        M.Irsyad Hibatullah                       ====')
% disp('====                         Jeremy Gohitua M.                        ====')
% disp('====                           Nuha Malihati                          ====')
% disp('====                         Firman Syaifuddin                        ====')
% disp('====                          Dwa Desa Warnana                        ====')
% disp('====                          Juan Pandu G.N.A                        ====')
% disp('====                                                                  ====')
% disp('==========================================================================')
% disp('==========================================================================')

    