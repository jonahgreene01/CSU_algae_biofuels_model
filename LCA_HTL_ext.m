function [net_process_lca, stage_breakdwn,stage_breakdwn_prc,consum_breakdwn, consum_breakdwn_prc, aware_outputs, traci_vec, cf_vec, nerc_out, bC_balance] = LCA_HTL_ext(co2_capt_lci,Cult_Out_Station_D,storage_fcn_out_D,HTL_process_mod_out_D,coord, direct_cap,file_id, aware_cult_1_D, aware_cult_2_D,transp_vec, annual_uptime, grid_red, bC_HTL_gas_D)

load background_data.mat aware_data 

% LCA module, take process outputs and load data for process consumables
cult_mean =  mean(Cult_Out_Station_D,2);
storage_mean =  mean(storage_fcn_out_D, 2);
HTL_process_mod_out =  mean(HTL_process_mod_out_D,2);
dewatering_energy = sum(cult_mean([11,14:16],:)); % cult_out(14:16) to include UV energy, cult_out(14:15) to exclude UV energy.
bC_HTL_gas = mean(bC_HTL_gas_D,2)*(12/44); % kg CO2 converted to kg C

co2_capt_lci(2:3) = co2_capt_lci(2:3)*cult_mean(10); % MJ natural gas per year and MJelec per year
co2_capt_lci(1) = co2_capt_lci(1)*cult_mean(10); % indirect kg CO2e per year
co2_capt_lci(3) = co2_capt_lci(3)/3.6; % kWh per year

fuel_output = HTL_process_mod_out(13); % Total MJ/yr
fuel_monthly =  HTL_process_mod_out(23); % Total MJ/month
monthly_conv_water = HTL_process_mod_out(24); %kg water/month
blue_wd_monthly = mean(aware_cult_1_D, 1); 
biomass_monthly = mean(aware_cult_2_D, 1); 

% Transportation and combustion mass balance (kg per MJ of fuel)
% [CO2; CH4; N2O; VOC; CO; NOx; PM10; PM2.5; SOx]
end_use = [0.000366653 0.071;
           4.54725E-07 0;
           5.4471E-09 1.78E-07;
           9.54914E-08 3.82E-05;
           5.71295E-07 9.42E-04;
           1.27208E-06 2.60E-05;
           3.84212E-08 7.75E-07;
           2.70104E-08 7.13E-07;
           2.137E-08 0.00E+00];

% Determine electricity impacts based on grid
[traci_vec, nerc_out] = grid_finder(coord);
traci_vec(1, 6) = traci_vec(1,6)*(1-grid_red); 

% Determine nutrient impacts based on transportation distances
[transp_impacts] = transp_calc(transp_vec);

% Carbon balance to determine additional biogenic loss - 
% all values in kg C/yr
bC_uptake = (-cult_mean(26))*annual_uptime*1000*(12/44); % kg C/yr 
bC_storage_loss = storage_mean(11)*(0.54); % kg C/yr - storage losses are subtracted from uptake credit
bC_HTL_CO2_vent = bC_HTL_gas; % kg C in CO2 emitted from HTL vent per year
bC_HTL_aqueous = cult_mean(26)*annual_uptime*1000*HTL_process_mod_out(25)*(12/44); %kg C/yr in aqueous phase (losses)
bC_fuel_comb = end_use(1,2)*fuel_output*(12/44); % kg C/yr in fuel combustion emissions - all fuels
bC_release = bC_storage_loss + bC_HTL_CO2_vent + bC_HTL_aqueous + bC_fuel_comb; 
net_bC = bC_uptake + bC_release;

    % Adjust carbon balance by adjusting C in char and aqueous
    if HTL_process_mod_out(20) == 0 % if no char yield from AFDW biomass (ash only)
        bC_char = 0; % carbon in char is zero
        bC_HTL_aqueous = bC_HTL_aqueous - net_bC; % adjust C in aqueous phase 
    else % if there is biochar
        bC_char = cult_mean(26)*annual_uptime*1000*(12/44)*HTL_process_mod_out(20); %bC in char = bC in biomass x char yield in wt% AFDW biomass
        net_bC = bC_uptake + bC_release + bC_char; 
        bC_HTL_aqueous = bC_HTL_aqueous - net_bC;
    end 

% Carbon Mass Balance Output
bC_balance = [bC_uptake bC_storage_loss bC_HTL_CO2_vent bC_HTL_aqueous bC_fuel_comb bC_char];  

% [1. Carbon capture, 2. Cultivation, 3. Dewatering, 4. Storage, 5. Conversion, 6. Transp, 7. Combustion, 8. Credits]
LCI_table = [0,0, 0, 0, 0, 0, 0,bC_uptake*(44/12)+bC_char*(44/12); % portion of CO2 uptake balanced by release of biogenic carbon, kg per yr
             0,cult_mean(8),0, 0, 0, 0, 0,0; % ammonia in kg per yr
             0,cult_mean(8),0, 0, 0, 0, 0,0; % ammonia transp in kg per yr
             0,cult_mean(9),0, 0, 0, 0, 0,0; % DAP in kg per yr
             0,cult_mean(9),0, 0, 0, 0, 0,0; % DAP transp in kg per yr
             co2_capt_lci(1),0,0,0,HTL_process_mod_out(22),end_use(1,:)*fuel_output,0; % direct CO2 in kg per yr
             0,cult_mean(18)*cult_mean(1),0,0, HTL_process_mod_out(4)/1000,0,0,0; % water consumption in m3 per yr
             co2_capt_lci(2)/35.396, 0,0,0, HTL_process_mod_out(3)/0.68,0,0,0; % natural gas consumption in m3 per yr
             0,0,0,0,HTL_process_mod_out(2),0,0,0;% hydrogen consumption in kg per yr
             co2_capt_lci(3),sum([cult_mean(12:13); cult_mean(17)]),dewatering_energy, sum(storage_mean(9:10)), HTL_process_mod_out(5),0,0,0; % electricity consumption in kWh per yr
             zeros(8,4), direct_cap(2:9)',end_use(2:9,:).*fuel_output,zeros(8,1);% direct emissions (not CO2) in kg per yr
             0,0,0,0,0,0,0,-HTL_process_mod_out(16); % ammonia credit in kg per yr
             0,0,0,0,0,0,0,-HTL_process_mod_out(16); % ammonia transp credit in kg per yr
             0,0,0,0,0,0,0,-HTL_process_mod_out(17);  % DAP credit in kg per yr
             0,0,0,0,0,0,0,-HTL_process_mod_out(17); % DAP transp credit in kg per yr
             0,0,0,0,bC_HTL_aqueous*(44/12),0,0,0;% kg CO2/yr from aqueous phase
             0,0,0,0,0,0,0,-bC_char*(44/12)];% biogenic CO2 stored in biochar kg CO2/yr 
        
% LCA data - Following TRACI 2.1 with GWP from the IPCC AR-6 Report: 
% 1 kg CO2 = 1 kg CO2eq,  1 kg fossil CH4 = 29.8 kg CO2eq, 
% 1 kg non fossil CH4 = 27.2 kg CO2eq, 1 kg N2O = 273 kg CO2eq
% Traci Impacts + water consumption: 
% [1. Respiratory effects (kg PM2.5), 2. Acidification (kg SO2eq), 3. Ecotoxicity (CTUe), 
%       4. Non carcinogenics (CTUh), 5. Carcinogenics (CTUh),
%       6. Global Warming Potential (kg CO2eq per MJ), 7. Smog formation (kg O3eq), 
%       8. ozone depletion (kg CFC-11), 9. Eutrophication (kg Neq),
%       10. Fossil fuel depletion (MJ surplus), 11. Water consumption (m3)] 
LCI_data = [0 0 0 0 0 1 0 0 0 0 0; % per kg of CO2
   0.00025	0.00206	13.54276 1.61E-07 8.13E-08 2.62741	0.04428	2.49E-08 0.00102 5.95224 1.73E-3; % per kg of ammonia
   transp_impacts(1,:),0; % per kg of ammonia transported
   1.96E-03 0.01456	101.35884 1.56E-06 1.76E-07 1.4729 0.08644 2.06E-08 0.01122 2.46E+00 3.70E-02; % per kg of DAP
   transp_impacts(2,:),0; % per kg of DAP transported
   0 0 0 0 0 1 0 0 0 0 0; % per kg of CO2
   0.000375291	0.000799737	2.297983096	5.23756E-08	1.97592E-08	0.209654243	0.008706602	1.30381E-08	0.001199824	0.219378355 1; % m3 of groundwater (no treatment)
   0.00011	0.00071	1.82E+00 2.99E-08 3.98E-08 0.50059	0.01681	8.32E-09 3.20E-04 6.37153 0.000; % per m3 of natural gas transported
   0.00095	0.00623	29.18753 3.53E-07 2.90E-07 11.39284 0.1122 2.67E-07 0.00241	29.30076 0.01098;% per kg of H2
   traci_vec; % per kWh of electricity
   0 0 0 0 0 29.8 0 0 0 0 0; % per kg of CH4
   0 0 0 0 0 273 0 0 0 0 0; % per kg of N2O
   0 0 0 0 0 0 3.595 0 0 0 0; % per kg of volatile organic compounds
   3.56E-04 0 0 0 0 0 0 0 0 0 0; % per kg of CO
   0.007222 0.7 0 0 0 0 0 0 0.04429 0 0; % per kg of NOx
   0 0 0 0 0 0 0 0 0 0 0; % per kg of PM10
   1 0 0 0 0 0 0 0 0 0 0; % per kg of PM2.5
   0.0611 1 0 0 0 0 0 0 0 0 0;% per kg of SOx
   0.00025	0.00206	13.54276 1.61E-07 8.13E-08 2.62741	0.04428	2.49E-08 0.00102 5.95224 1.73E-3; % per kg of ammonia
   transp_impacts(1,:),0; % per kg of ammonia transported
   1.96E-03 0.01456	101.35884 1.56E-06 1.76E-07 1.4729 0.08644 2.06E-08 0.01122 2.46E+00 3.70E-02; % per kg of DAP
   transp_impacts(2,:),0; % per kg of DAP transported
   0 0 0 0 0 1 0 0 0 0 0; % per kg of CO2 released 
   0 0 0 0 0 1 0 0 0 0 0]; % per kg of CO2
 
% Impacts by stage
co2_capt_lca = zeros(size(LCI_table,1), size(LCI_data,2)); 
cult_lca = zeros(size(LCI_table,1), size(LCI_data,2)); 
dew_lca = zeros(size(LCI_table,1), size(LCI_data,2)); 
storage_lca = zeros(size(LCI_table,1), size(LCI_data,2)); 
htl_lca = zeros(size(LCI_table,1), size(LCI_data,2)); 
transp_lca = zeros(size(LCI_table,1), size(LCI_data,2)); 
combustion_lca = zeros(size(LCI_table,1), size(LCI_data,2)); 
credits_lca =  zeros(size(LCI_table,1), size(LCI_data,2)); 

    for i = 1:size(LCI_data,2)
        co2_capt_lca(:,i) =  LCI_table(:,1).*LCI_data(:,i);
        cult_lca(:,i) = LCI_table(:,2).*LCI_data(:,i);
        dew_lca(:,i) = LCI_table(:,3).*LCI_data(:,i);
        storage_lca(:,i) =  LCI_table(:,4).*LCI_data(:,i);
        htl_lca(:,i) =  LCI_table(:,5).*LCI_data(:,i);
        transp_lca(:,i) =  LCI_table(:,6).*LCI_data(:,i);
        combustion_lca(:,i) =  LCI_table(:,7).*LCI_data(:,i);
        credits_lca(:,i) =  LCI_table(:,8).*LCI_data(:,i);
    end

co2_capt_lca = co2_capt_lca/fuel_output;  % normalize by functional unit
cult_lca = cult_lca/fuel_output; % normalize by functional unit
dew_lca = dew_lca/fuel_output; % normalize by functional unit
storage_lca = storage_lca/fuel_output; % normalize by functional unit
htl_lca = htl_lca/fuel_output; % normalize by functional unit
transp_lca = transp_lca/fuel_output; % normalize by functional unit
combustion_lca = combustion_lca/fuel_output; % normalize by functional unit
credits_lca = credits_lca/fuel_output; % normalize by functional unit
co2_credits = credits_lca(1,:);
ammonia_credits = sum(credits_lca(19:20,:)); 
dap_credits = sum(credits_lca(21:22,:)); 
char_credits = credits_lca(24,:); 

net_process_lca = sum([co2_capt_lca;cult_lca;dew_lca; storage_lca; htl_lca; transp_lca; combustion_lca; credits_lca], 1); % system traci impacts per MJ of fuel
gross_process_lca = sum([co2_capt_lca;cult_lca;dew_lca; storage_lca; htl_lca; transp_lca; combustion_lca], 1); 

% Break by processes
stage_breakdwn = [sum(co2_capt_lca,1); sum(cult_lca,1); sum(dew_lca,1); sum(storage_lca,1); sum(htl_lca,1); sum(transp_lca,1); sum(combustion_lca,1); co2_credits; ammonia_credits; dap_credits; char_credits]; 
stage_breakdwn_prc = [sum(co2_capt_lca,1); sum(cult_lca,1); sum(dew_lca,1); sum(storage_lca,1); sum(htl_lca,1); sum(transp_lca,1); sum(combustion_lca,1) ; sum(credits_lca,1)]*100./gross_process_lca; % in % total

% Break by consumables
for j = 1:size(LCI_data,1)
    consum_breakdwn(j,:) = sum([co2_capt_lca(j,:);cult_lca(j,:); dew_lca(j,:); storage_lca(j,:); htl_lca(j,:);transp_lca(j,:); combustion_lca(j,:) ; credits_lca(j,:)]); %#ok<AGROW> 
end
   consum_breakdwn_prc = consum_breakdwn*100./gross_process_lca; 

% Water scarcity
% find CF for corresponfing file
cf_vec = aware_data(aware_data(:,1) == str2double(file_id), 2:end);
% Cultivation
wsf_monthly_cult = [blue_wd_monthly(1:12).*cf_vec(2:end)./biomass_monthly(1:12), blue_wd_monthly(13)*cf_vec(1)/biomass_monthly(13)]; % [m3eq per tonne] 12 months and annual average
% Well to Wheels
wsf_monthly_fuels = [(blue_wd_monthly(1:12) + monthly_conv_water/1000).*cf_vec(2:end)/fuel_monthly, (blue_wd_monthly(13) + HTL_process_mod_out(4)/1000)*cf_vec(1)/fuel_output]; % [m3eq per MJ] 12 months and annual average

% Build Water Scarcity Outputs
aware_outputs = [wsf_monthly_cult; wsf_monthly_fuels]';

end %End of LCA Function
 

%% Grid finder function
function [traci_vec, nerc_out] = grid_finder(coord)

load background_data.mat S

% grid_shape matrix contains lat and long data for each of the NERC regions
% 1 = HICC; 2 = NPCC; 3 = RFC, 4 = WECC, 5 = TRE, 6 = SERC; 7 = MRO
% nerc_rgn = ["HICC", "NPCC", "RFC","WECC", "TRE", "SERC", "MRO"]; 

% Electricity Impacts matrix from ecoinvent_391_cutoff_regionalized, water from Lee et al. 2018 
% electricity high voltage production mix Cutoff, U, for 2018
% impacts at a NERC region level
% [Respiratory effects (kg PM2.5), acidification (kg SO2eq), ecotoxicity (CTUe), 
%       non carcinogenics (CTUh), carcinogenics (CTUh), GWP (kg CO2-eq), Smog formation (kg O3eq), 
%           ozone depletion (kg CFC-11), eutrophication (kg Neq),
%               Fossil fuel depletion (MJ surplus) Water consumption (m3)] - normalized by 1 kWh
elec_lci =  [0.001666443	0.005921149	4.079776851	1.21E-07	4.06E-08	0.9262526	0.063606275	9.88E-09	0.00342215	1.363982554	0; %HICC = hawaii island coordinating council  
            3.55E-05	0.000223711	0.647183731	2.21E-08	8.05E-09	0.23065419	0.005066792	1.14E-09	9.82E-05	0.572521208	0.00024;% NPCC = northeast power coordinating council
           0.000369183	0.001385715	2.213925871	7.61E-08	2.61E-08	0.503547502	0.013544062	6.69E-09	0.001568108	0.528636508	0.00143; % RFC = reliability first council
            0.00099995	0.000716388	2.893731424	9.13E-08	3.23E-08	0.388210782	0.011933808	2.40E-09	0.002549539	0.482137489	0.00354; % WECC = western electricity coordinating council
          0.001086464	0.000906571	3.318392754	9.59E-08	3.56E-08	0.47005713	0.009287823	2.65E-09	0.002731942	0.684279781	0.002953;% TRE = Texas reliability entity
           0.000426248	0.00100036	2.154228023	7.78E-08	2.60E-08	0.501282329	0.01244617	5.14E-09	0.001624789	0.709355056	0.004372;% SERC = southeastern reliability council
            0.002091226	0.0016689	5.620188375	1.71E-07	5.83E-08	0.512936408	0.016365639	2.59E-09	0.005092143	0.3559336	0.002877228]; %  MRO = midwest reliability orfanization

% Determine NERC region
for i = 1:7
     [idx_check, on] = inpolygon(round(coord(2),1), round(coord(1),1), S(i).X,S(i).Y);
    if idx_check == 1 || on == 1
        break
    end
end

if i == 10 && idx_check == 0 && on == 0 % Inderteminate region, take avg of SERC,MRO
    traci_vec = mean(elec_lci([2,5],:),1);
    nerc_out = "SERC/MRO"; 
else
    nerc_out = string(S(i).name);
    traci_vec = elec_lci(i,:); % get appropriate LCI data
end

end

%% Transportaion Impacts Function
function[transp_impacts] = transp_calc(transp_vec)
% transp_vec = [distance in freight train, distance in freight inland
% waterway, lorry distance, freight sea container ship distance]

% Environmental impacts of transportation in (km*tonne)
% [Respiratory effects (kg PM2.5), acidification (kg SO2eq), ecotoxicity (CTUe), 
%       non carcinogenics (CTUh), carcinogenics (CTUh), GWP (kg CO2-eq), Smog formation (kg O3eq), 
%           ozone depletion (kg CFC-11), eutrophication (kg Neq),
%               Fossil fuel depletion (MJ surplus) Water consumption (m3)]
transp_lci = [4.48E-05	0.00054	0.50501	6.82E-09	1.05E-08	0.06099	0.01636	9.04E-10	9.05E-05	0.10019; % transp freight train
2.84E-05	0.00042	0.24313	3.65E-09	4.44E-09	0.05142	0.01291	7.67E-10	6.38E-05	0.0835; % transp freight inland waterway
9.23E-05	0.00065	1.37755	3.56E-08	1.12E-08	0.1533	0.01698	2.61E-09	0.00015	0.30631; % transp lorry unspecified
1.55E-05	0.00026	0.03581	4.99E-10	6.44E-10	0.01016	0.00478	1.65E-10	1.24E-05	0.01834]; % transp freight sea container ship

transp_impacts(1,:) = sum(transp_lci.*transp_vec(1,:)')/1000; % per kg ammonia
transp_impacts(2,:) = sum(transp_lci.*transp_vec(2,:)')/1000; % per kg dap

end