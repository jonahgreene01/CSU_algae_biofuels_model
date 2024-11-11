function [decay] = Night_Respiration_ORP(GHI, CX, night_resp, volume, Temp_eff)


if (GHI < 5)
    
    decay_rate = (log(1 - night_resp))/10; %yields the decay rate 1/(dark period) as a negative decimal
    decay_specific = (CX*(decay_rate/3600)); 
    decay = (decay_specific*volume)*Temp_eff; 

else 
    
    decay = 0; 
    
end 

%based on Edmunson and Huessemann night respiration study. Study assumes 10
%hour dark period for the conversion from percentage to decay rate

end
