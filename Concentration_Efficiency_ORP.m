function [Conc_eff] = Concentration_Efficiency_ORP(CX, ODC, depth)
%This function determines the of impact of algal concentration on the
%average light intensity hitting the culture

gpl= CX/1000; %gpl = grams per liter
OD = gpl/ODC; %OD = optical density (750 nm)

Conc_eff = (1-exp(-depth*OD))./(OD*depth);
%Conc_eff = 100*exp(-1.574*OD)/100; 

end

