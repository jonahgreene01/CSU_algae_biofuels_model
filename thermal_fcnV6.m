function [TR,evap_matrix, Water_con, M_Evap] = thermal_fcnV6(weather_matrix,geom_vec)
% Runs thermal model, output is pond temperature over time
% inputs are dimensions of pond and weather data
% inputs_vec = [length,width,depth]';
% weather_matrix = [Tamb, RH,GHI,WNDSPD]';

% Decompose Weather Matrix
T_amb = weather_matrix(:,1)+273.15;
RH =  weather_matrix(:,2);
GHI= weather_matrix(:,3);
WNDSPD= weather_matrix(:,4);

% Decomopose inputs vector
length_orp = geom_vec(1);
width = geom_vec(2);
depth = geom_vec(3);

%Define unchanging parameters
h = 3600; %Time step is 1 hour
pi = 3.14159;
rho_algae = 1000.0;  %kg/m3
cp_algae = 4184.0; % specific heat 
area = pi*0.25*(width^2) + width*(length_orp-width); % Pond area
volume = area*depth; % volume in m3
perimeter = (3.14*width)+(2*(length_orp-width)); % Pond perimeter
Tg = mean(T_amb); % Ground temperature K

n1 = 1;
n2 = length(T_amb);

TR = zeros(length(T_amb),1);
mass_flux = zeros(length(T_amb),1);
solar = zeros(length(T_amb),1);
convection = zeros(length(T_amb),1);
atmospheric = zeros(length(T_amb),1);
conduction = zeros(length(T_amb),1);
evaporation = zeros(length(T_amb),1);
inflow = zeros(length(T_amb),1);
reradiation = zeros(length(T_amb),1);
Water_consump = zeros(length(T_amb),1);

for i = n1:(n2-1) %Only go to n-1 because program reads to i+1
    
    if (i == n1)
        TR(i,1) = T_amb(i,1);
    else
        TR(i,1) = TR(i,1);
    end
    
    
    %Thermal Model - 4th Order Runge-Kutta Scheme
    %Compute Air Properties at this time step using the film temperature
    [rho_air, ~, ~, ~, dyn_visc_air, ~, ~]=Air_Properties((TR(i,1)+T_amb(i,1))/2);
    
    %Compute thermal mass based on algal properties
    Therm_mass = volume*rho_algae*cp_algae;
    
    %Fist RK4 Fluxes
    [Q_Solar] = Direct_Solar_HT_ORP(GHI(i,1));
    [Q_Convection,~] = Convection_HT_ORPV3(area,perimeter, TR(i,1), T_amb(i,1), WNDSPD(i,1), dyn_visc_air, rho_air);
    [Q_Longwave_Atmo] = Longwave_Atmo_HT_ORPV4(T_amb(i,1));
    [Q_Ground] = Ground_ConductionV2(TR(i,1),Tg,h);
    [M_Evap, Q_Evap,~,~] = EvaporationV2(area, perimeter, TR(i,1), T_amb(i,1), RH(i,1), WNDSPD(i,1),length_orp, dyn_visc_air, rho_air);
    [Q_Inflow] = Inflow_HT_ORP(M_Evap, cp_algae, TR(i,1), T_amb(i,1));
    [Q_Rerad] = Reradiation_HT_ORP(TR(i,1));
    
    %RK4_1 time derivative of temperature - K1=dT1
    dT1 =(Q_Evap + Q_Convection +Q_Inflow + Q_Longwave_Atmo + Q_Solar + Q_Ground + Q_Rerad)*area/Therm_mass;
    T1 = TR(i,1) + 0.5*h*dT1;
    
    %Air Properties - Thermal mass already defined
    [rho_air, ~, ~, ~, dyn_visc_air, ~, ~]=Air_Properties((T1+((T_amb(i,1)+T_amb(i+1,1))/2))/2);
    
    %Second RK4 Fluxes
    [Q_Solar] = Direct_Solar_HT_ORP((GHI(i,1)+GHI(i+1,1))/2);
    [Q_Convection,~] = Convection_HT_ORPV3(area,perimeter, T1, ((T_amb(i,1)+T_amb(i+1,1))/2), (WNDSPD(i,1)+WNDSPD(i+1,1))/2, dyn_visc_air, rho_air);
    [Q_Longwave_Atmo] = Longwave_Atmo_HT_ORPV4((T_amb(i,1)+T_amb(i+1,1))/2);
    [Q_Ground] = Ground_ConductionV2(T1,Tg,h);
    [M_Evap, Q_Evap,~,~] = EvaporationV2(area, perimeter,T1, (T_amb(i,1)+T_amb(i+1))/2, (RH(i,1)+RH(i+1,1))/2, (WNDSPD(i,1)+WNDSPD(i+1,1))/2,length_orp, dyn_visc_air, rho_air);
    [Q_Inflow] = Inflow_HT_ORP(M_Evap, cp_algae, T1, (T_amb(i,1)+T_amb(i+1, 1))/2);
    [Q_Rerad] = Reradiation_HT_ORP(T1);
    
    % RK4_2 time derivative - K2 = dT2
    dT2 = (Q_Evap +  Q_Convection + Q_Inflow + Q_Longwave_Atmo + Q_Solar + Q_Ground + Q_Rerad)*area/Therm_mass;
    T2= TR(i,1) + 0.5*h*dT2;
    
    %Air Properties - Thermal mass already defined
    [rho_air, ~, ~, ~, dyn_visc_air, ~, ~] = Air_Properties((T2+((T_amb(i,1)+T_amb(i+1,1))/2))/2);
    
    %Third RK4 Fluxes
    [Q_Solar] = Direct_Solar_HT_ORP((GHI(i,1)+GHI(i+1,1))/2);
    [Q_Convection,~] = Convection_HT_ORPV3(area,perimeter,T2, (T_amb(i,1)+T_amb(i+1,1))/2, (WNDSPD(i,1)+WNDSPD(i+1,1))/2, dyn_visc_air, rho_air);
    [Q_Longwave_Atmo] = Longwave_Atmo_HT_ORPV4((T_amb(i,1)+T_amb(i+1,1))/2);
    [Q_Ground] = Ground_ConductionV2(T2,Tg,h);
    [M_Evap, Q_Evap,~] = EvaporationV2(area,perimeter,T2, ((T_amb(i,1)+T_amb(i+1))/2), (RH(i,1)+RH(i+1,1))/2, (WNDSPD(i,1)+WNDSPD(i+1,1))/2,length_orp, dyn_visc_air, rho_air);
    [Q_Inflow] = Inflow_HT_ORP(M_Evap, cp_algae, T2, (T_amb(i,1)+T_amb(i+1, 1))/2);
    [Q_Rerad] = Reradiation_HT_ORP(T2);
    
    %RK4_3 time derivative - K3 = dT3
    dT3 = (Q_Evap + Q_Convection + Q_Inflow + Q_Longwave_Atmo + Q_Ground + Q_Solar + Q_Rerad)*area/Therm_mass;
    T3 = TR(i,1) + h*dT3;
    
    %Air Properties - Thermal mass already defined
    [rho_air, ~, ~, ~, dyn_visc_air, ~, ~] = Air_Properties((T3+ T_amb(i+1,1))/2);
    
    %Fourth RK4 Fluxes
    [Q_Solar] = Direct_Solar_HT_ORP(GHI(i+1,1));
    [Q_Convection,~] = Convection_HT_ORPV3(area,perimeter,T3, T_amb(i+1,1), WNDSPD(i+1,1), dyn_visc_air, rho_air);
    [Q_Longwave_Atmo] = Longwave_Atmo_HT_ORPV4(T_amb(i+1,1));
    [Q_Ground] = Ground_ConductionV2(T3,Tg,h);
    [M_Evap, Q_Evap,~,~] = EvaporationV2(area,perimeter,T3, T_amb(i+1,1), RH(i+1,1), WNDSPD(i+1,1),length_orp, dyn_visc_air, rho_air);
    [Q_Inflow] = Inflow_HT_ORP(M_Evap, cp_algae, T3, T_amb(i+1,1));
    [Q_Rerad] = Reradiation_HT_ORP(T3);
    
    %RK4_1 time derivative
    dTdT = (Q_Evap + Q_Convection + Q_Inflow + Q_Longwave_Atmo + Q_Ground + Q_Solar + Q_Rerad)*area/Therm_mass;
    TR(i+1,1)= TR(i,1) + (1/6)*h*(dT1+2.0*(dT2+dT3)+dTdT);
end

    for i = 1:n2
    %  Air Properties
    [rho_air, ~, ~, ~, dyn_visc_air, ~, ~] = Air_Properties((TR(i,1)+T_amb(i,1))/2);
    [mass_flux(i,1),evaporation(i,1),~,~] = EvaporationV2(area,perimeter, TR(i,1), T_amb(i,1), RH(i,1), WNDSPD(i,1), length_orp, dyn_visc_air, rho_air);
    Water_consump(i,1) = mass_flux(i,1)*area*h; %[kg/m2*s]*[m2]*[3600 s] = [kg/hr]
    
    if Water_consump(i,1) < 0 || TR(i,1) < 273.15
        Water_consump(i,1) = 0;
    else
        Water_consump(i,1) = Water_consump(i,1);
    end
    
    % Heat Fluxes
    [solar(i,1)] = Direct_Solar_HT_ORP(GHI(i,1));
    [convection(i,1),~] = Convection_HT_ORPV3(area,perimeter, TR(i,1), T_amb(i,1), WNDSPD(i,1), dyn_visc_air, rho_air);
    [atmospheric(i,1),~] = Longwave_Atmo_HT_ORPV4(T_amb(i,1));
    [conduction(i,1)] = Ground_ConductionV2(TR(i,1),Tg,h);
    [inflow(i,1)] = Inflow_HT_ORP(M_Evap, cp_algae, TR(i,1), T_amb(i,1));
    [reradiation(i,1)] = Reradiation_HT_ORP(TR(i,1));   
    end

M_Evap = mass_flux; 
Water_con = Water_consump(n1:n2,1); %in kg/hr
%Water_consump_ave = sum(Water_con)/area; %L / m2/year
  
    % Seasonal Evaporation Rate (cm/day)
    evap_matrix = [sum(Water_con(1756:3963))/((3963-1756)/24),... 
    sum(Water_con(3964:6171))/((6171-3964)/24),...
    sum(Water_con(6175:8354))/((8354-6175)/24),...
    sum([Water_con(1:1755); Water_con(8355:8759)])/((1754+(8759-8355))/24),...
    sum(Water_con)/((n2-n1)/24)]/10/area;
  
%heat_fluxes = [solar convection atmospheric conduction evaporation inflow reradiation];
  
end

