function [Temp_eff] =Temp_Efficiency_ORP (TR, T_opt, T_min, T_max)
% Cardinal Model based on Rosso et al. 1993
if TR < T_min
    Temp_eff = 0; 
end 

if TR >= T_min && TR <= T_max
    
    g_of_T = (T_opt - T_min)*(TR - T_opt); 
    f_of_T = (T_opt - T_max)*(T_opt + T_min - 2*TR); 
    Temp_eff = ((TR - T_max)*(TR - T_min)^2)/((T_opt - T_min)*(g_of_T - f_of_T));
end 

if TR > T_max
    Temp_eff = 0; 
end 


end

