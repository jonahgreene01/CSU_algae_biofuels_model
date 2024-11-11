function [Q_Solar]= Direct_Solar_HT_ORP(GHI)
%This function defines the heat flux from direct solar radiation
absorp = 0.9;
f_a = 0.015; %fraction of sunlight converted to chemical energy during photosynthesis
Q_Solar = (1-f_a)* GHI*absorp; 

% Source: Yadala and Cremaschi, 2016

end

