function [Q_Longwave_Atmo,epsilon_air] = Longwave_Atmo_HT_ORPV4(T_amb)
%This function defines the heat flux from longwave atmospheric radiation
sigma = 5.67*(10^-8); %Stefan Boltzmann constant (W/m^2*k^4)
pa = 3385.5.*exp(-8.0929 + 0.97608.*(T_amb + 42.607 - 273.15).^0.5); % vapor pressure at air temperature Pa
pa = pa*.0075; %mmHg
% Epsilon air from Brunt equation
epsilon_air = 0.6+ 0.031*sqrt(pa); 
Q_Longwave_Atmo = epsilon_air*sigma.*((T_amb).^4);
end


