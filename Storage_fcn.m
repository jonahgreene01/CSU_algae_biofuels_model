function [storage_fcn_out,ann_ave_centr, biomass_balance_st, comp_ann_ave_centr] = Storage_fcn(prot, carb, lip, ash, centr_out, centr_out_dens, degradation, pump_eff)
% Anaerobic storage module estimates losses from anaerobic storage and
% resulting mass flow to HTL

%Determine annual average output from centrifuge
ann_ave_centr = mean(centr_out); %ave flow from centr in kg/hr

%Pre-allocate variables
comp_ann_ave_centr = zeros(1, 10); 
to_storage = zeros(1,8760);
from_storage = zeros(1, 8760); 
storage_vol = zeros(1, 8760); 

k = 1; 

while k < 10  %takes about 3-4 iterations to converge
for i = 1:8760
    
    if centr_out(i) > ann_ave_centr
        %divert to storage in covered settling ponds
        to_storage(i) = centr_out(i) - ann_ave_centr; 
        from_storage(i) = 0; 
        else 
        to_storage(i) = 0; 
        from_storage(i) = ann_ave_centr - centr_out(i); 
    end 
    
    %calculate volume in storage with each hour - in m3
    if i == 1
        storage_vol(i) = (to_storage(i) - from_storage(i))/centr_out_dens; %m3 in storage
        else 
        storage_vol(i) = (to_storage(i) - from_storage(i))/centr_out_dens + storage_vol(i-1); %m3 in storage
    end 
    
end

direct_HTL = sum(centr_out) - sum(to_storage); 

%caclulate required settling pond volume for storage
num_settling_ponds = ceil(max(storage_vol)/1000); %1000 m3 per settling pong; ceil fcn - round up to nearest integer

%Calculate biomass loss and resulting composition 
dw_biomass_stored = sum(to_storage)*0.20/(1-ash); %kg DW (includes ash)
dw_prot_stored = dw_biomass_stored*prot; 
dw_carb_stored = dw_biomass_stored*carb;
dw_lip_stored = dw_biomass_stored*lip; 
dw_ash_stored = dw_biomass_stored*ash; 
dw_biomass_deg = dw_biomass_stored*degradation; %losses from storage - set percentage on model interface

%Biomass output from storage (no degredation of lip, prot, ash - assumes only carbs degrade from anaerobic conditions)
dw_biomass_from_storage = dw_biomass_stored - dw_biomass_deg; 
dw_lip_from_storage = dw_lip_stored;
dw_prot_from_storage = dw_prot_stored; 
dw_ash_from_storage = dw_ash_stored; 
dw_carb_from_storage = dw_biomass_from_storage - dw_lip_from_storage - dw_prot_from_storage - dw_ash_from_storage; 

%Resulting composition of biomass pulled from anaerobic storage
stored_bm_prot = dw_prot_from_storage/dw_biomass_from_storage; 
stored_bm_lip = dw_lip_from_storage/dw_biomass_from_storage; 
stored_bm_carb = dw_carb_from_storage/dw_biomass_from_storage;
stored_bm_ash = dw_ash_from_storage/dw_biomass_from_storage; 

ann_ave_centr = (sum(centr_out) - dw_biomass_deg)/8760; 

%Track and compare annual average centrifuge flow to ensure convergence
comp_ann_ave_centr(k) = ann_ave_centr; 
k = k+1; 

end 

%resulting average composition 
total_bm_out = direct_HTL*0.20/(1-ash) + dw_biomass_from_storage; 

HTL_feed_carb = (direct_HTL*0.20/(1-ash))/total_bm_out*carb + dw_biomass_from_storage/total_bm_out*stored_bm_carb; 
HTL_feed_prot = (direct_HTL*0.20/(1-ash))/total_bm_out*prot + dw_biomass_from_storage/total_bm_out*stored_bm_prot; 
HTL_feed_lip = (direct_HTL*0.20/(1-ash))/total_bm_out*lip + dw_biomass_from_storage/total_bm_out*stored_bm_lip;
HTL_feed_ash = (direct_HTL*0.20/(1-ash))/total_bm_out*ash + dw_biomass_from_storage/total_bm_out*stored_bm_ash;

%Calculate slurry pumping energy to and from stroage
pump_pwr_to_storage = ((to_storage./centr_out_dens)./3600).*centr_out_dens.*9.81.*98.05./pump_eff./1000; %kW - 98.05 is pumping head from Li et al., 2020; 70% pump eff. 
pump_pwr_from_storage = ((from_storage./centr_out_dens)./3600).*centr_out_dens.*9.81.*98.05./pump_eff./1000; %kW - 98.05 is pumping head from Li et al., 2020; 70% pump eff. 

ann_pump_pwr_to_kWh = sum(pump_pwr_to_storage); 
ann_pump_pwr_from_kWh = sum(pump_pwr_from_storage); 

%MASS BALANCE!
%Biomass to storage (DW)
biomass_balance_st(1,1) = sum(centr_out)*0.20/(1-ash);
biomass_balance_st(2,1) = biomass_balance_st(1,1)*prot;
biomass_balance_st(3,1) = biomass_balance_st(1,1)*carb;
biomass_balance_st(4,1) = biomass_balance_st(1,1)*lip;
biomass_balance_st(5,1) = biomass_balance_st(1,1)*ash;

%Storage Losses (DW)
biomass_balance_st(1,2) = dw_biomass_deg;
biomass_balance_st(2,2) = (dw_prot_from_storage - dw_prot_stored)*-1; 
biomass_balance_st(3,2) = (dw_carb_from_storage - dw_carb_stored)*-1; 
biomass_balance_st(4,2) = (dw_lip_from_storage - dw_lip_stored)*-1; 
biomass_balance_st(5,2) = (dw_ash_from_storage - dw_ash_stored)*-1; 

%Biomass out of storage and to conversion
biomass_balance_st(1,3) = total_bm_out;
biomass_balance_st(2,3) = total_bm_out*HTL_feed_prot; 
biomass_balance_st(3,3) = total_bm_out*HTL_feed_carb; 
biomass_balance_st(4,3) = total_bm_out*HTL_feed_lip; 
biomass_balance_st(5,3) = total_bm_out*HTL_feed_ash; 

%New TSS of AFDW out of storage
water_in = sum(centr_out)*0.80; 
AFDW_out = sum(centr_out)*.20-dw_biomass_deg;
new_conc = AFDW_out/(AFDW_out + water_in); 

%Output matrix
storage_fcn_out(1,1) = "HTL Feed (kg/hr at 20% TSS)"; 
storage_fcn_out(1,2) = ann_ave_centr; 
storage_fcn_out(2,1) = "Composition to HTL"; 
storage_fcn_out(2,2) = "wt%"; 
storage_fcn_out(3,1) = "Proteins"; 
storage_fcn_out(3,2) = HTL_feed_prot; 
storage_fcn_out(4,1) = "Carbohydrates"; 
storage_fcn_out(4,2) = HTL_feed_carb; 
storage_fcn_out(5,1) = "Lipids"; 
storage_fcn_out(5,2) = HTL_feed_lip; 
storage_fcn_out(6,1) = "Ash"; 
storage_fcn_out(6,2) = HTL_feed_ash;
storage_fcn_out(7,1) = "Required Settling Pond Volume"; 
storage_fcn_out(7,2) = max(storage_vol); %m3
storage_fcn_out(8,1) = "Number of Settling Ponds"; 
storage_fcn_out(8,2) = num_settling_ponds; %number of 1000 m3 settling ponds
storage_fcn_out(9,1) = "Annual Pumping Energy to Storage (kWh)"; 
storage_fcn_out(9,2) = ann_pump_pwr_to_kWh; %number of 1000 m3 settling ponds
storage_fcn_out(10,1) = "Annual Pumping Energy from Storage (kWh)"; 
storage_fcn_out(10,2) = ann_pump_pwr_from_kWh; %number of 1000 m3 settling ponds
storage_fcn_out(11,1) = "Annual Biomass Loss from Storage kg/yr";
storage_fcn_out(11,2) = dw_biomass_deg; 
storage_fcn_out(12,1) = "New AFDW TSS%";
storage_fcn_out(12,2) = new_conc; 

end 
