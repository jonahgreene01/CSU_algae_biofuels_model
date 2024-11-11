function [Q_Ground] = Ground_ConductionV2(TR,Tg,dt)
%This function determines the heat flux between the ground and the ORP
% Inputs: diff_concrete = thermal diffusivity of concrete or soil (m2/s)
% k_concrete = thermal conductivity of concrete or soil (W/m-K)
% diffusivity values range from 9E-7 to 1E-6 , 
% Conductivity values range from 0.25 to 2.7. Ground temperature assumed to equal at the mean annual ambient temp. 



diff_concrete = 7.9E-6 ; % m2 per seconds
k_concrete = 1.7; % W per m*K
Q_Ground = -k_concrete*(TR-Tg)/sqrt(pi*diff_concrete*dt);

end


