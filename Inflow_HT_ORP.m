function [Q_Inflow] = Inflow_HT_ORP(M_Evap, cp_algae, TR, T_amb)
%Calculates the Heat Flux due to the temperature difference between the
%algae pond and the makeup water. 

Q_Inflow = cp_algae*M_Evap*(T_amb - TR); % [J/kg*K] *[kg/m2*s] * [K] = [J/m2*s] = [w/m2]

% Yadala and Cremaschi 2016
end

