function [Q_Convection,h_conv] = Convection_HT_ORPV3(area,perimeter, TR, T_amb, WNDSPD, dyn_visc_air, rho_air)
%This function calculates the convection heat flux from the ORP


kin_visc_a = dyn_visc_air/rho_air; 
L_c = area/perimeter; 
L_f = area/perimeter; 
alpha_a = 2.2*10^-5;
Pr_a = kin_visc_a/alpha_a;
lamda_a = 2.6*10^-2;
beta = 1/((TR+T_amb)/2);  % coefficient of expansion
n = 3; % Can be changed to 4 and 7/2 according to Incropera 


Re_L = WNDSPD*L_c/kin_visc_a; 

Gr = 9.81*beta*abs(TR-T_amb)*(L_f^(3))/(kin_visc_a^2);
 
if  order_mag(Gr) < order_mag(Re_L^2)  %(Re_L^2) > Gr % forced convection dominates
    if Re_L < (3*10^5) %then the flow is laminar and use:
        Nu_L = 0.75*(Re_L.^0.5)*(Pr_a^(1/3));
    elseif Re_L > (5*10^5) %then the flow is turbulent and use:
        Nu_L = 0.015*(Re_L.^0.8)*(Pr_a^(1/3));
    else   %Average the values sherwood values for laminar and turbulent
        Nu_L = ((0.65*(Re_L.^0.5).*(Pr_a^(1/3)))+(0.045*(Re_L.^0.8)*(Pr_a^(1/3))))/2;
    end
   h_conv = Nu_L*lamda_a/L_c; 
elseif order_mag(Gr) > order_mag(Re_L^2)  % natural convection dominates
    Ra_L = Gr*Pr_a;
    if Ra_L < 10^7 && Ra_L > 10^4 % laminar
        Nu_L = 0.54*Ra_L^(1/4);
    elseif Ra_L < 10^11 && Ra_L > 10^7 % turbulent
        Nu_L = 0.15*Ra_L^(1/3);
    else
        Nu_L = ((0.54*Ra_L^(1/4))+(0.15*Ra_L^(1/3)))/2;
    end
    h_conv = Nu_L*lamda_a/L_f; 
elseif order_mag(Gr) == order_mag(Re_L^2)  % mixed based on Incropera
    if Re_L < (3*10^5) %then the flow is laminar and use:
        Nu_f = 0.75*(Re_L^0.5)*(Pr_a^(1/3));
    elseif Re_L > (5*10^5) %then the flow is turbulent and use:
        Nu_f = 0.015*(Re_L^0.8)*(Pr_a^(1/3));
    else   %Average the values sherwood values for laminar and turbulent
        Nu_f = ((0.65*(Re_L^0.5)*(Pr_a^(1/3)))+(0.015*(Re_L^0.8)*(Pr_a^(1/3))))/2;
    end
    Ra_L = Gr*Pr_a;
    if Ra_L < 10^7 && Ra_L > 10^4 % laminar
        Nu_n = 0.54*Ra_L^(1/4);
    elseif Ra_L < 10^11 && Ra_L > 10^7 % turbulent
        Nu_n = 0.15*Ra_L^(1/3);
    else
        Nu_n = ((0.54*Ra_L^(1/4))+(0.15*Ra_L^(1/3)))/2;
    end
    %Calculate the convection coefficient given the Nusselt number
    Nu_L = ((Nu_f^n)+ (Nu_n^n))^(1/n);
    h_conv = Nu_L*lamda_a/L_c; 
end

%Calculate the convective heat transfer
Q_Convection = h_conv*(T_amb - TR);
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
