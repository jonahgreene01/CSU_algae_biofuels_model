function [m_evap,Q_evap,Re_L,Gr] = EvaporationV2(area, perimeter,TR,Tamb,RH,WNDSPD,length_orp,dyn_visc_air,rho_air)
% Combined evaporation function including free and forced evaporation
% Free evaporation is modeled following the methods explained in Adams et
% al. 1990
% Forced evaporation is modeled following the methods explained in
% Sartori et al. 1987 and 1999
% When both forced and free convection are equal, the evaporation rate is
% calculated squaring them according to Adams et al. 1990


% Inputs:
% pond temperature (K)
% ambient temperature (K)
% relative humidity (%)
% wind speed (m/s)
% pond length (m)


% Declare constants
g = 9.81;% gravity m/s6
L_f = area/perimeter; % characteristic length
L_n = area/perimeter;  % characteristic length
D_w_a = 2.4*10^-5;
M_water = 0.018; %kg/mol
R = 8.314; %Universal gas constant Pa*m3/mol*K
hfg_water = 2.45*10^6;
n = 3; % Can be changed to 4 and 7/2 according to Incropera 

% correlation to adjust for height from prevailing windspeed
if length_orp < 100
    z =  0.5 ; 
    zo = 10;% 
else
    z = 3;
    zo = 10;
end

beta = 1./((TR+Tamb)/2);  % coefficient of expansion
kin_visc = dyn_visc_air./rho_air;  % kinematic viscosity


pw = 3385.5.*exp(-8.0929 + 0.97608.*(TR + 42.607 - 273.15).^0.5); % vapor pressure at water temp Pa
pa = 3385.5.*exp(-8.0929 + 0.97608.*(Tamb + 42.607 - 273.15).^0.5); % vapor pressure at air temperature Pa

%Calculate the Schmidt number
Sch_L = kin_visc./D_w_a;
% Correct for wind speed height
if WNDSPD == 0
    na = 0;
else
    na = (.37-.0881*log10(abs(WNDSPD)))/(1-.0881*log10(z/zo));
end

V = WNDSPD*(z/zo)^na;

%Calculate the Reynold's and Grasshoff Number
Gr = (g*beta*abs(TR-Tamb)*(L_n^(3)))/(kin_visc^2);
Ra_L =  (g*beta*abs(TR-Tamb)*(L_n^3))/(dyn_visc_air*kin_visc);
Re_L = L_f*V/kin_visc;

if order_mag(Gr) < order_mag(Re_L^2)   % forced convection dominates (Gr(i,1)/Re_L(i,1)^2) < 1E-1
    
    if Re_L < (3*10^5) %then the flow is laminar and use:
        Sh_L = 0.65*(Re_L^0.5)*(Sch_L^(1/3)); % =0.628
    elseif Re_L > (5*10^5) %then the flow is turbulent and use:
        Sh_L = 0.045*(Re_L.^0.8).*(Sch_L^(1/3)); % =.035
    else   %Average the values sherwood values for laminar and turbulent
        Sh_L = ((0.65*(Re_L^0.5).*(Sch_L^(1/3)))+(0.045*(Re_L^0.8).*(Sch_L^(1/3))))/2;
    end
    
    Kf = Sh_L.*D_w_a./L_f;
    
    %Now we can calculate the evaporation rate in kg/m2*s
    m_evap = Kf*((pw/TR)-((RH/100)*pa/Tamb))*M_water/R; %Evaporation rate in kg/s*m2
    Q_evap = -m_evap*hfg_water; %Evaporative heat transfer in W/m2,
    
elseif   order_mag(Gr) > order_mag(Re_L^2) % natural convection dominates  (Gr(i,1)/Re_L(i,1)^2) > 10
    if Ra_L < 8E6 && Ra_L > 2.6E4 % Laminar natural Convection
        Sh_L = 0.54*Ra_L^0.25; % c= 0.54
    elseif Ra_L > 8E6
        Sh_L = 0.15*Ra_L^(1/3); % d = 0.15
    else
        Sh_L = ((0.54*Ra_L^0.255)+(0.15*Ra_L^0.327))/2;
    end
    
    Kn = Sh_L*D_w_a/L_n;
    
    %Now we can calculate the evaporation rate in kg/m2*s
    m_evap = Kn*((pw/TR)-((RH/100)*pa/Tamb))*M_water/R; %Evaporation rate in kg/s*m2
    Q_evap = -m_evap*hfg_water; %Evaporative heat transfer in W/m2,
    
elseif   order_mag(Gr) == order_mag(Re_L^2)  %(Gr(i,1)/Re_L(i,1)^2 )== 1 % both natural and forced are considered
    
    if Re_L < (3*10^5) %then the flow is laminar and use:
        Sh_f = 0.65*(Re_L^0.5)*(Sch_L^(1/3));
    elseif Re_L > (5*10^5) %then the flow is turbulent and use:
        Sh_f = 0.045*(Re_L^0.8)*(Sch_L^(1/3));
    else   %Average the values sherwood values for laminar and turbulent
        Sh_f = ((0.65*(Re_L^0.5).*(Sch_L^(1/3)))+(0.045*(Re_L^0.8).*(Sch_L^(1/3))))/2;
    end
      
    if Ra_L < 8E6 && Ra_L > 2.6E4 % Laminar natural mass transfer
        Sh_n = 0.54*Ra_L^0.25;
    elseif Ra_L > 8E6
        Sh_n = 0.15*Ra_L^(1/3); % Turbulent natural mass transfer
    else
        Sh_n = ((0.54*Ra_L^0.255)+(0.15*Ra_L^0.327))/2;
    end
    % Sherwood number for mixed 
    Sh_L =  ((Sh_f^n)+ (Sh_n^n))^(1/n);
    K = Sh_L.*D_w_a./L_f;
    m_evap = K*((pw/TR)-((RH/100)*pa/Tamb))*M_water/R; %Evaporation rate in kg/s*m2
    Q_evap = -m_evap*hfg_water; %Evaporative heat transfer in W/m2,
    
end

if TR < 273.15
    m_evap = 0;
    Q_evap = 0;
end
end

function [n] = order_mag( val, base )
%Order of magnitude of number for specified base. Default base is 10.
%order(0.002) will return -3., order(1.3e6) will return 6.
%Author Ivar Smith
if nargin < 2
    base = 10;
end

for i = 1:length(val)
    if val(i,1) == 0 
    n(i,1) =0;
    else 
    n(i,1) = round(floor(log(abs(val(i,1)))./log(base)));
    end
end
end

