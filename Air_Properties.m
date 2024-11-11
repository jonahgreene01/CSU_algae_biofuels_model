function [rho_air, k_air, hfg_air, cp_air, dyn_visc_air, pr_air, spec_vol_air] = Air_Properties(T)
%This function returns the properties of air at a given temperature
%accepts temperature to provide the density, thermal conductivity, heat of water vaporization, specific heat, 
%dynamic viscosity, prandtl number, and specific volume

rho_air = 101.315*1000.0./(287.058.*T); %kg/m3
k_air = 0.02624*(T./300).^0.8646; %W/m*K
hfg_air = (-2E-05.*T.^3 + 0.0176.*T.^2 - 7.8474.*T + 3721.5); %kJ/kg
cp_air = 1002.5 + (275E-06).*(T - 200).^2; %J/kg*K
dyn_visc_air = (1.458E-06).*(T.^1.5)./(T + 110.4);
pr_air = cp_air.*dyn_visc_air./k_air; %dimless
spec_vol_air = 1./rho_air; %m3/kg

end
