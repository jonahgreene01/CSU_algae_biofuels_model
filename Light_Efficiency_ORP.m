function [Light_eff] = Light_Efficiency_ORP(GHI, Conc_eff, I_sat)
%This function determines light efficiency
% I_sat in W/m2 
I_o = GHI*.45; % W/m2

I_ave = I_o*Conc_eff; 
Light_eff = (I_ave/I_sat)*exp(1-(I_ave/I_sat)); 

end

