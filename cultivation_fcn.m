function [CX, centr_out, centr_out_ann, cult_out, water_vol, CO2_demand_hourly, aquifer_hrly, area, total_recycled, num_settling_ponds, settlers_out, membrane_out, aware_cult, nutrients, recycled_membr_centr, dep_rep_scaling, harvest_water, Conc_at_Harvest, light_eff_track, aquifer, rec_algae, biomass_balance, energy_balance, dewatering_mat, check_biomass] = cultivation_fcn(weather_data,geom_vec,strain_vec, harvest_vec, comp_vector, rep_dep, prod_scale, carb, lip, prot)

% Geometry inputs
length_orp = geom_vec(1);
width = geom_vec(2);
depth = geom_vec(3);
max_depth = geom_vec(4);
min_depth = geom_vec(5);

% Weather inputs
TR = weather_data(:,1); % Pond temp in K
GHI =  weather_data(:,2); % W per m2
M_Evap = weather_data(:,3); % kg per m2 s
rainfall =  weather_data(:,4); % m per hr

%Strain Parameters
T_opt = strain_vec(1);
T_min = strain_vec(2);
T_max = strain_vec(3);
ODC = strain_vec(4);
I_sat = strain_vec(5)/4.56; % convert from micromol/m2/s to W/m2
night_resp= strain_vec(6);
phi = strain_vec(7); %photon efficiency is 12 g biomass/8 mol photons

%Define unchanging parameters
h = 3600; %Time step is 1 hour
n1 = 1;
n2 = length(GHI);
area = width*length_orp; % m2
volume = area*depth; %m3
Conversion = 4.56 * 10^-6 ; %mol/m2/s per W/m2

%Initialize Variables
marker = 0;

%Nutrients and CO2
Algae_C_content = comp_vector(1);
CO2_eff = comp_vector(2);
co2_deliv_power = 0.0391; % kWh per kg CO2
ash_content = comp_vector(9);

% Preallocate variables
water_vol = zeros(length(GHI),1);
depth_dyn_track = zeros(length(GHI),1);
salt_cx = zeros(length(GHI),1);
harvest_algae = zeros(length(GHI),1);
harvested_algae = zeros(length(GHI),1);
recycled_algae = zeros(length(GHI),1);
harvest_water = zeros(length(GHI),1);
CX = zeros(length(GHI),1);
Conc_at_Harvest = zeros(length(GHI),1);
salt_mass = zeros(length(GHI),1);
aquifer =zeros(length(GHI),1);
aquifer_hrly = zeros(length(GHI),1); 
blowdown = zeros(length(GHI),1);
blowdown_algae = zeros(length(GHI),1);
nutrients = zeros(length(GHI),3);
num_settling_ponds = zeros(length(GHI),1);
membrane_energy = zeros(length(GHI),1);
settlers_out = zeros(n2,1);
membrane_out  = zeros(n2,1);
centr_energy = zeros(length(GHI),1);
uv_energy = zeros(length(GHI),1);
makeup_energy  = zeros(length(GHI),1);
algae_to_con = zeros(length(GHI),1);
pump_energy = zeros(length(GHI),1);
turn_on = ones(n2,1);
algae_ponds = zeros(n2,1);
centr_out = zeros(n2,1); 
centr_out_algae = zeros(n2,1); 
centr_out_ann = zeros(n2,1);
rec_algae = zeros(n2,1); 
dep_rep_scaling = zeros(n2,1); 
light_eff_track = zeros(n2,1); 

for i = 1:length(GHI)-1
    
    %Assign IBC_o, harvest_conc, harv_days, i_salt, salt_max_cx 
    IBC_o = harvest_vec(1);
    harvest_conc = harvest_vec(2);
    harvest_days = harvest_vec(3);
    i_salt = harvest_vec(4); % g per m3
    salt_max_cx = harvest_vec(5); % g per m3
        
    if (i == n1)
        CX(i,1) = IBC_o;
        water_vol(i,1) = volume;
        salt_cx(i,1) = i_salt;
        salt_mass(i,1) = i_salt*volume;
    else
        CX(i,1) = CX(i,1);
        water_vol(i,1) = water_vol(i,1);
        salt_cx(i,1) = salt_cx(i,1);
        salt_mass(i,1) = salt_mass(i,1);
    end
    
    %Biological Growth Model
    [Temp_eff] = Temp_Efficiency_ORP (TR(i,1), T_opt, T_min, T_max);
    [Conc_eff] = Concentration_Efficiency_ORP(CX(i,1), ODC, depth);
    [Light_eff] = Light_Efficiency_ORP(GHI(i,1), Conc_eff, I_sat);
    [decay] = Night_Respiration_ORP(GHI(i,1), CX(i,1), night_resp, water_vol(i,1), Temp_eff);
    
    light_eff_track(i,1) = Light_eff; 

    % Nutrient deplete/replete productivity scaling function
    switch rep_dep 
        case 'Nutrient Replete'
        dep_rep_scaling(i,1) = 0;
        case 'Nutrient Deplete' 
        dep_rep_scaling(i,1) = -0.1616*cosd((i/n2)*360)-0.3116;
    end

    dCXdt = phi*Light_eff*Temp_eff*GHI(i,1)*0.458*0.95*Conversion*area*(1+dep_rep_scaling(i,1))*prod_scale + decay;
    CX(i+1,1) = ((dCXdt*h)/water_vol(i,1) + CX(i,1));
    
    %Marker and Operational Hours Counting
    marker = marker + 1;
    
    % Check pond depth
    depth_dyn = water_vol(i,1)*100/area; % dynamic depth in cm
    depth_dyn_track(i+1,1) = depth_dyn; 

    %Harvesting Sequence
    if (depth_dyn > max_depth  || (CX(i+1,1) >= harvest_conc) || (mod(marker, (harvest_days*24.0))== 0)  ||(i==(n2-1))) && turn_on(i+1,1) == 1
        
        % Calculate total precipitation in next harvest (p_future)
        [P]= p_future(i,area,rainfall,n2,harvest_days) ; % m3
        
        % Harvest mass and stream
        sc_percent = 1-(IBC_o/CX(i+1,1))*(volume/water_vol(i,1)); % fraction of pond volume that must be harvested
        harvest_algae(i+1,1) = CX(i+1,1)*(sc_percent*water_vol(i,1))/marker/1000; % Harvested algae kg/hr
        TSS = (CX(i+1,1)/1000/1000); % total suspended solids percent
        harvest_water(i+1,1) = (harvest_algae(i+1,1)/TSS); % Water in harvest stream kg per hour
        
        % Dewatering
        [pond_states,energy_matrix,dewatering_mat,~,~] = dewatering(area,CX(i+1,1),IBC_o,harvest_algae(i+1,1), harvest_water(i+1,1),volume,P,water_vol(i,1),salt_cx(i,1),i_salt,salt_max_cx,marker,max_depth,min_depth,ash_content);
        recycled_membr_centr(i+1-marker:i+1, 1) = dewatering_mat(5,5); % kg per hr total flow from membr and centr (including water, algae, ash, salt)

        % Parse dewatering streams
        %pond_states = [CX_np1 vol_refill salt_cx_np1 algae_added aquifer_makeup salt_mass];
        water_vol(i+1,1) = pond_states(2); % m3 - set to vol_refill determined in dewatering function
        aquifer(i+1,1) = pond_states(5);% m3 total aquifer makerup
        aquifer_hrly(i+1-marker:i+1,1) = aquifer(i+1,1)/marker; %m3/h flowrate from aquifer
        total_recycled(i+1-marker:i+1, 1) = sum(dewatering_mat(4:5,5)); % kg per hr total flow from settlers, mem, and centr
        settlers_out(i+1-marker:i+1, 1) = dewatering_mat(1,5); % kg per hr total flow from settlers to membranes
        membrane_out(i+1-marker:i+1, 1) = dewatering_mat(2,5); % kg per hr total flow from membranes to centrifuges
        blowdown(i+1,1) = dewatering_mat(6,1)*marker/1000; % blowdown water m3 - singular event
        algae_to_con(i+1,1) = dewatering_mat(3,2)*marker/1000;% tonnes AFDW algae out of centrifuge for harvesting event
        centr_out(i+1-marker:i+1, 1) = dewatering_mat(3,5); % kg per hr total flow from centrifuge

        %For mass balance check
        harvested_algae(i+1,1) = harvest_algae(i+1,1)*marker; % kg AFDW algae removed from ponds for harvesting event
        centr_out_algae(i+1,1) = dewatering_mat(3,2)*marker; % kg AFDW algae from centrifuge for harvesting event
        centr_out_ann(i+1,1) = dewatering_mat(3,5)*marker; 
        recycled_algae(i+1,1) = sum(dewatering_mat(4:5,2))*marker; %kg AFDW algae in recycle water from settlers, mem, and cent
        blowdown_algae(i+1,1) = dewatering_mat(6,2)*marker; % algae losses to blowdown in kg/hr - singular event

        % Settling pond quantitiy
        num_settling_ponds(i+1,1) = energy_matrix(1); 

        % Salts balance
        salt_cx(i+1,1) = pond_states(3); % g per m3
        salt_mass(i+1,1) = pond_states(6); % g
        
        % Parse energy vector
        %energy_matrix = [settling_pond_no membrane_energy centr_energy uv_energy pump_dew pump_makeup];
        membrane_energy(i+1,1) =  energy_matrix(2); % kWh/day
        centr_energy(i+1,1) =  energy_matrix(3); % kWh/day
        uv_energy(i+1,1) =  energy_matrix(4);% kWh/day
        pump_energy(i+1,1) = energy_matrix(5); % dewatering pumping kW
        makeup_energy(i+1,1) = energy_matrix(6); % evap makeup pumping in kW

        % Calculate net biomass harvested
        algae_added = pond_states(4); % grams - algae added to reach IBC_o
        rec_algae(i,1) = algae_added;
        algae_ponds(i+1,1) =  CX(i+1,1)*(sc_percent*water_vol(i,1)) - algae_added; % g AFDW harvested from ponds
        Conc_at_Harvest(i+1,1) = CX(i+1,1);
        
        % Restart Ponds
        CX(i+1,1) = pond_states(1);
        marker = 0;
        
    else
        
        Conc_at_Harvest(i+1,1) = 0;
        harvest_water(i+1,1) = 0;
        harvest_algae(i+1,1) = 0;
        algae_ponds(i+1,1) = 0; 
        water_vol(i+1,1) = water_vol(i,1)- (M_Evap(i,1)*area*h/1000) + rainfall(i,1)*area; % m3
        salt_cx(i+1,1) = salt_mass(i,1)/water_vol(i+1,1);
        salt_mass(i+1,1) = salt_mass(i,1);
        
    end
    
    % Nutrient and CO2 Mass Balance
    [dap_ammonia_CO2] = nutrient(CX(i+1,1),CX(i,1),water_vol(i,1),comp_vector);
    nutrients(i+1,:) = dap_ammonia_CO2;   
    
end

%Set final centr_out to zero to not add extra hour of flow to mass balance
centr_out(n2,1) = 0; 

%DEWATERING BIOMASS MASS BALANCE CHECK!
check_biomass = [sum(harvested_algae), sum(recycled_algae), sum(blowdown_algae), sum(algae_ponds)/1000, sum(centr_out_algae), sum(centr_out_ann)]; 

%net harvested DW biomass (including ash) and composition [prot; carb; lip; ash]
biomass_balance(1,1) = sum(algae_ponds)/1000/(1-ash_content); %net algae harvested
biomass_balance(2,1) = biomass_balance(1,1)*prot; 
biomass_balance(3,1) = biomass_balance(1,1)*carb; 
biomass_balance(4,1) = biomass_balance(1,1)*lip; 
biomass_balance(5,1) = biomass_balance(1,1)*ash_content; 

%centrifuge out DW biomass and composition
biomass_balance(1,2) = sum(centr_out)*0.20/(1-ash_content);
biomass_balance(2,2) = biomass_balance(1,2)*prot;
biomass_balance(3,2) = biomass_balance(1,2)*carb;
biomass_balance(4,2) = biomass_balance(1,2)*lip;
biomass_balance(5,2) = biomass_balance(1,2)*ash_content;

%ENERGY BALANCE CHECK!
energy_balance = [membrane_energy centr_energy uv_energy pump_energy makeup_energy];

% Water demand vector - net consumption (blue water)
wd_net =  aquifer - blowdown; % m3
wd_net(wd_net<0) = 0;

% % Calculate seasonal outputs
idx_seasonal = [1755 3963 6171 8355];
vec_idx(1:idx_seasonal(1),1) = 1; % winter
vec_idx((idx_seasonal(1)+1):idx_seasonal(2),1) = 2; % spring
vec_idx((idx_seasonal(2)+1):idx_seasonal(3),1) = 3; % summer
vec_idx((idx_seasonal(3)+1):idx_seasonal(4),1) = 4;% fall
vec_idx((idx_seasonal(4)+1):n2,1) = 1; % 2nd part of winter

% % monthly blue water demand for aware calculations
blue_wd_monthly = [sum(wd_net(1:744)),...
    sum(wd_net(745:1416)),...
    sum(wd_net(1417:2160)),...
    sum(wd_net(2161:2880)),...
    sum(wd_net(2881:3624)),...
    sum(wd_net(3625:4344)),...
    sum(wd_net(4345:5088)),...
    sum(wd_net(5089:5832)),...
    sum(wd_net(5833:6552)),...
    sum(wd_net(6553:7296)),...
    sum(wd_net(7297:8016)),...
    sum(wd_net(8017:8759)),...
    sum(wd_net)];

% %  monthly biomass output for aware calculations
algae_ponds_tonnes = algae_ponds/1000/1000; 

biomass_monthly = [sum(algae_ponds_tonnes(1:744)),...
    sum(algae_ponds_tonnes(745:1416)),...
    sum(algae_ponds_tonnes(1417:2160)),...
    sum(algae_ponds_tonnes(2161:2880)),...
    sum(algae_ponds_tonnes(2881:3624)),...
    sum(algae_ponds_tonnes(3625:4344)),...
    sum(algae_ponds_tonnes(4345:5088)),...
    sum(algae_ponds_tonnes(5089:5832)),...
    sum(algae_ponds_tonnes(5833:6552)),...
    sum(algae_ponds_tonnes(6553:7296)),...
    sum(algae_ponds_tonnes(7297:8016)),...
    sum(algae_ponds_tonnes(8017:8759)),...
    sum(algae_ponds_tonnes)];

aware_cult = [blue_wd_monthly; biomass_monthly]; 

% Operational days
operational_days = [accumarray(vec_idx,turn_on,[],@sum); sum(turn_on)]'/24; % days

% Biomass seasonal
biomass_seasonal = [accumarray(vec_idx,algae_ponds,[],@sum); sum(algae_ponds)]'; % grams

% Average annual and seasonal areal productivity [winter spring summer fall annual average];
Areal_prod_matrix = biomass_seasonal./operational_days/area;  % g m^-2 day^-1
Vol_prod_matrix = biomass_seasonal./operational_days/(volume*1000); % g L^-1 day^-1

% Clean consumable outputs, make sure to account only consumables
% used when ponds were operating
nutrients(nutrients < 0) = 0; % clear nightime consumption kg
nutrients(turn_on == 0 ,1) = 0;
nutrients(turn_on == 0 ,2) = 0;
nutrients(turn_on == 0 ,3) = 0;

CO2_demand_hrly = algae_ponds*Algae_C_content*(44/12)/CO2_eff/1000; % kg - hourly resolved but singular harvesting event values
CO2_demand_hourly = nutrients(:,3); 
consumables(1,1:2) = [sum(nutrients(:,1)) sum(nutrients(:,2))]/operational_days(5)/1000; % tonne per day
consumables(1,3) = sum(CO2_demand_hrly)/operational_days(5)/1000;  % tonne CO2 per day
CO2_uptake = sum(algae_ponds*Algae_C_content*(44/12)/1000)/operational_days(5)/1000; % tonne CO2 per day

% Calculate annual energy requirements
paddlewheel_energy = 55.050*(area/10000); % kWh per day - 55.050 kWh hectare^-1 day^-1
co2_deliv_energy = consumables(1,3)*1000*co2_deliv_power; % kWh per day
membrane_energy = mean(nonzeros(membrane_energy)); % kWh per day
centr_energy = mean(nonzeros(centr_energy)); % kWh per day
uv_energy = mean(nonzeros(uv_energy)); % kWh per day
makeup_energy = mean(nonzeros(makeup_energy))*24; % kWh per day
pump_energy = mean(nonzeros(pump_energy))*24; % kWh per day

energy_matrix = [membrane_energy centr_energy uv_energy pump_energy paddlewheel_energy co2_deliv_energy makeup_energy]*operational_days(5); % kWh per year

% Calculate evaporation rate
evap_rate = M_Evap*area*h; %[kg/m2*s]*[m2]*[3600 s] = [kg/hr]
for i = 1:n2
    if evap_rate(i,1) < 0 || TR(i,1) < 273.15
        evap_rate(i,1) = 0;
    else
        evap_rate(i,1) = evap_rate(i,1);
    end 
end

% Seasonal Evaporation Rate (cm/day)
evap_matrix = [accumarray(vec_idx,evap_rate,[],@mean); mean(evap_rate)]'*24/10/area;

% Blue and green water demand in m3/year
green_wd =  sum(rainfall(turn_on == 1))*area; % m3/year
blue_wd_seasonal = [accumarray(vec_idx, wd_net, [], @sum); sum(wd_net)]'; % m3 per season

% Prepare outputs
algae_out = sum(algae_to_con);% tonne AFDW per year
blue_wf_annual = (blue_wd_seasonal(5))/algae_out; % m3 per tonne AFDW
green_wf = green_wd/algae_out; % m3 per tonne AFDW


%Build final outputs matrix
cult_out(1,1) = "Annual Biomass Production (tonnes AFDW)"; 
cult_out(1,2) =  sum(algae_ponds)/1000/1000; %tonnes AFDW
cult_out(2,1) = "Annual Average Areal Productivity"; 
cult_out(2,2) = Areal_prod_matrix(1,5);
cult_out(3,1) = "Annual Average Volumetric Productivity";
cult_out(3,2) = Vol_prod_matrix(1,5);
cult_out(4,1) = "Spring Productivity"; 
cult_out(4,2) = Areal_prod_matrix(1,2); 
cult_out(5,1) = "Summer Productivity"; 
cult_out(5,2) = Areal_prod_matrix(1,3); 
cult_out(6,1) = "Fall Productivity"; 
cult_out(6,2) = Areal_prod_matrix(1,4); 
cult_out(7,1) = "Winter Productivity"; 
cult_out(7,2) = Areal_prod_matrix(1,1); 
cult_out(8,1) = "Ammonia Consumption (kg/yr)";
cult_out(8,2) = sum(nutrients(:,2)); 
cult_out(9,1) = "DAP Consumption (kg/yr)";
cult_out(9,2) = sum(nutrients(:,1)); 
cult_out(10,1) = "CO2 Consumption (tonnes CO2/yr)";
cult_out(10,2) = consumables(1,3)*operational_days(5);
cult_out(11,1) = "Pumping Energy (kWh/yr)";
cult_out(11,2) = energy_matrix(1,4); 
cult_out(12,1) = "Paddlewheel Energy (kWh/yr)";
cult_out(12,2) = energy_matrix(1,5); 
cult_out(13,1) = "Sparging Energy (kWh/yr)"; 
cult_out(13,2) = energy_matrix(1,6); 
cult_out(14,1) = "Membrane Energy (kWh/yr)"; 
cult_out(14,2) = energy_matrix(1,1); 
cult_out(15,1) = "Centrifuge Energy (kWh/yr)";
cult_out(15,2) = energy_matrix(1,2); 
cult_out(16,1) = "UV Sterilization Energy (kWh/yr)";
cult_out(16,2) = energy_matrix(1,3); 
cult_out(17,1) = "Makeup Water Pumping Energy (kWh/yr)";
cult_out(17,2) = energy_matrix(1,7); 
cult_out(18,1) = "Freshwater Consumption (m3 per tonne AFDW) at 20% solids";
cult_out(18,2) = blue_wf_annual;
cult_out(19,1) = "Rainwater Consumption (m3 per tonne AFDW) at 20% solids";
cult_out(19,2) = green_wf;
cult_out(20,1) = "Biomass Output at 20% solids (tonne AFDW per year)";
cult_out(20,2) = algae_out;
cult_out(21,1) = "Winter Evaporation (cm per day)";
cult_out(21,2) = evap_matrix(1);
cult_out(22,1) = "Spring Evaporation (cm per day)";
cult_out(22,2) = evap_matrix(2);
cult_out(23,1) = "Summer Evaporation (cm per day)";
cult_out(23,2) = evap_matrix(3);
cult_out(24,1) = "Fall Evaporation (cm per day)";
cult_out(24,2) = evap_matrix(4);
cult_out(25,1) = "Operational Days";
cult_out(25,2) = operational_days(5);
cult_out(26,1) = "Carbon Uptake (tonnes CO2 per day)";
cult_out(26,2) = CO2_uptake;
end

% Required functions
% Dewatering function mass and energy balance
function [pond_states,energy_matrix,summary_streams,conc,pump_electricity] = dewatering(area,CX,~,harvest_algae, harvest_water,volume_o,P_m3,water_vol,salt_cx,i_salt,salt_max_cx,marker,max_depth,min_depth,ash_content)
% Water volume in m3 and volume in m3
% TSS (%), harvest_algae (kg/hr), harvest_water (kg/hr),
% P (m3), volume_o (m3), water_vol (m3), salt_cx in (g per m3)
% salt_max_cx (g per mr)

% Concentration parameters for each dewatering stage
dewater_1 = 10; % g/L on an AFDW basis
dewater_2 = 130; % g/L on an AFDW basis
dewater_3 = 200; % g/L on an AFDW basis

% Efficiency parameters
sep_eff_1 = 0.90; 
sep_eff_2 = .995; 
sep_eff_3 = 0.97;
rec_eff = 0.97;

% Energy consumption parameters
biofloc_res_time = 4; % hrs
membr_en_in = 0.04; % kWh per m3
cent_en_in = 1.35; % kWh per m3 
UV_en_in = 2.71/1000; % kWh per m3

% Pumping power from Davis et al. NREL model
inlet_biofloc = 0.0189; % kWh per m3 stream 
rec_water_ponds = 0.01767; % kWh per m3 
sett_to_membranes = 0.1289; % kWh per m3 stream 
membr_to_centr = 0.0194; % kWh per m3 stream
pump_membr_centr_rec = 0.1838; % kWh per m3 water
aquifer_pump = 0.2571; % kWh per m3 water

% Current state and forecasted rain
depth = water_vol/area*100; % cm
p_depth = P_m3/area*100; %cm

%adjust algae flow rate for ash content as pond CX is calcuated as AFDW
%harvest_algae = harvest_algae/(1-ash_content); 

% Pre-Allocate variables
summary_streams = zeros(6,5); % [water biomass ash salt total] in kg per hr

% Inlet stream to dewatering in kg/hr (water  algae  ash  salt)
% Nutrients are assumed to be completed assimilated by the culture
inlet_dew_stream = [harvest_water harvest_algae harvest_algae*(ash_content/(1-ash_content)) harvest_water*salt_cx/1e6]; % in kg per hr [water AFDW_algae ash salt]

% % Settlers - settlers are co-located with ponds, not in central
% dewatering system (Davis et al.)
% Algae and ash balance
recycle_streams(1,2) = harvest_algae*(1-sep_eff_1); % AFDW algae in recycled water (kg per hr) (losses)
summary_streams(1,2) = inlet_dew_stream(1,2) - recycle_streams(1,2); % AFDW biomass out to membranes, including losses (kg per hr)
recycle_streams(1,3) = recycle_streams(1,2)/(1-ash_content)*ash_content; % Ash in recycled stream (kg per hr)
summary_streams(1,3) = summary_streams(1,2)/(1-ash_content)*ash_content;  % Ash in stream exiting bioflocculation (kg per hr)

% New concentration 
summary_streams(1,5) = summary_streams(1,2)*1000/dewater_1; % Mainstream out of settlers, kg per hr: algae_afdw[kg per hr]*1000[g/kg]/10 [g AFDW per L]

% Water and Salt Balance (kg per hr)
summary_streams(1,1) = summary_streams(1,5) - sum(summary_streams(1,2:3)); %  Water out stream from settlers to membranes = Total stream - algae - ash
recycle_streams(1,1) = inlet_dew_stream(1,1) - summary_streams(1,1); % Recycled water = harvest water - out from settlers
recycle_streams(1,4) = recycle_streams(1,1)*salt_cx/1e6; % Salt in recycled water
summary_streams(1,4) = inlet_dew_stream(1,4) - recycle_streams(1,4); % Salt exiting settlers to membranes

% Concentrations in biofloc streams [g AFDW per L,  %wt, g salt per L] 
conc(1,:) = [ inlet_dew_stream(1,2)*1000/inlet_dew_stream(1,1),  inlet_dew_stream(1,2)*100/sum(inlet_dew_stream),  inlet_dew_stream(1,4)*1000/inlet_dew_stream(1,1)]; % Concentrations in inlet stream
conc(2,:) =  [ recycle_streams(1,2)*1000/recycle_streams(1,1),  recycle_streams(1,2)*100/sum(recycle_streams(1,:)),  recycle_streams(1,4)*1000/recycle_streams(1,1)]; % Concentrations in settlers recycle stream
conc(3,:) = [ summary_streams(1,2)*1000/summary_streams(1,1),  summary_streams(1,2)*100/summary_streams(1,5),  summary_streams(1,4)*1000/summary_streams(1,1)]; % Concentration in settlers output stream

% Check for blowdown 
if salt_cx > salt_max_cx  || (depth + p_depth) > max_depth
    blowdown_stream = recycle_streams(1,:); % Do not recycle and dispose recycle stream
    recycle_streams(1,:) = [0 0 0 0]; 
else
    % Recycled Water
    blowdown_stream =zeros(1,4);
end

% Number of settlers needed
settling_pond_no = round(sum(inlet_dew_stream)*biofloc_res_time/1000/1000); % each has capacity of 1000 m3

% Electricity to pump biomass from ponds to settlers and
% clarified water back to ponds
pump_electricity(1,1) = inlet_biofloc*sum(inlet_dew_stream)/1000; % kW
pump_electricity(1,2) = rec_water_ponds*sum(recycle_streams(1,:))/1000; % kW

% % Membrane Filtration
% Algae and ash balance (kg per hr)
recycle_streams(2,2) = summary_streams(1,2)*(1-sep_eff_2); % Algae in recycled stream from membrane (losses)
summary_streams(2,2) = summary_streams(1,2) - recycle_streams(2,2); % Net biomass out to centrifuge, including losses (kg per hr)
recycle_streams(2,3) = recycle_streams(2,2)/(1-ash_content)*ash_content; % Ash in recycled stream from membrane (kg per hr)
summary_streams(2,3) = summary_streams(2,2)/(1-ash_content)*ash_content;  % Ash in stream exiting membrane to centrifuge (kg per hr)

% New concentration 
summary_streams(2,5) = summary_streams(2,2)*1000/dewater_2; % Mainstream out of membranes, kg per hr: algae_afdw[kg per hr]*1000[g/kg]/130 [g AFDW per L]

% Water and Salt Balance (kg per hr)
summary_streams(2,1) = summary_streams(2,5) - sum(summary_streams(2,2:3)); %  Water out from membranes = Total stream out - algae - ash
recycle_streams(2,1) = summary_streams(1,1) - summary_streams(2,1); % Recycled water from membranes to ponds
recycle_streams(2,4) = recycle_streams(2,1)*salt_cx/1e6; % Salt in recycled water
summary_streams(2,4) = summary_streams(1,4) - recycle_streams(2,4); % Salt exiting membrane

% Concentrations in membrane streams [g AFDW per L,  %wt, g salt per L] 
conc(4,:) =  [ recycle_streams(2,2)*1000/recycle_streams(2,1),  recycle_streams(2,2)*100/sum(recycle_streams(2,:)),  recycle_streams(2,4)*1000/recycle_streams(2,1)]; % Concentrations in membrane recycle stream
conc(5,:) = [ summary_streams(2,2)*1000/summary_streams(2,1),  summary_streams(2,2)*100/summary_streams(2,5),  summary_streams(2,4)*1000/summary_streams(2,1)]; % Concentration in membrane output stream

% Electricity consumed by membrane filtration and to pump from settlers to membrane
membrane_energy = (summary_streams(1,5)*24/1000)*membr_en_in;  % kWh per day
pump_electricity(1,3) = sett_to_membranes*summary_streams(1,5)/1000; % kW

% % Centrifugation Process 
% Algae and ash balance (kg per hr)
recycle_streams(3,2) = summary_streams(2,2)*(1-sep_eff_3); % Algae in recycled stream (losses)
summary_streams(3,2) = summary_streams(2,2) - recycle_streams(3,2); % Net biomass out to conversion, including losses (kg per hr)
recycle_streams(3,3) = recycle_streams(3,2)/(1-ash_content)*ash_content; % Ash in recycled stream (kg per hr)
summary_streams(3,3) = summary_streams(3,2)/(1-ash_content)*ash_content;  % Ash in stream exiting centrifuge (kg per hr)

% New concentration 
summary_streams(3,5) = summary_streams(3,2)*1000/dewater_3; % Mainstream out of centrifuge, kg per hr: algae_afdw + ash [kg per hr]*1000[g/kg]/200 [g DW per L]

% Water and Salt Balance (kg per hr)
summary_streams(3,1) = summary_streams(3,5) - sum(summary_streams(3,2:3)); %  Water out from cent. stream = Total stream - algae - ash
recycle_streams(3,1) = summary_streams(2,1) - summary_streams(3,1); % Recycled water
recycle_streams(3,4) = recycle_streams(3,1)*salt_cx/1e6; % Salt in recycled water
summary_streams(3,4) = summary_streams(2,4) - recycle_streams(3,4); % Salt exiting centrifuge

% Add recycled water, algae, and ash to summary streams
summary_streams(4,:) = [recycle_streams(1,:) sum(recycle_streams(1,:))]; % Recycled from settlers stage
summary_streams(5,:) = [sum(recycle_streams(2:3,:)) sum(sum(recycle_streams(2:3,:)))]; % Recycled from membrane and centrifuge stages

% Add blowdown to summary streams
summary_streams(6,:) = [blowdown_stream sum(blowdown_stream(1,:))];

% Concentrations in centrifuge streams [g AFDW per L,  %wt, g salt per L] 
conc(6,:) =  [recycle_streams(3,2)*1000/recycle_streams(3,1),  recycle_streams(3,2)*100/sum(recycle_streams(3,:)),  recycle_streams(3,4)*1000/recycle_streams(3,1)]; % Concentrations in centrifuge recycle stream
conc(7,:) = [summary_streams(3,2)*1000/summary_streams(3,1),  summary_streams(3,2)*100/summary_streams(3,5),  summary_streams(3,4)*1000/summary_streams(3,1)]; % Concentration in centrifuge output stream
conc(8,:) = [summary_streams(5,2)*1000/summary_streams(5,1),  summary_streams(5,2)*100/summary_streams(5,5),  summary_streams(5,4)*1000/summary_streams(5,1)]; % Concentration in total recycled water from membrane and centrifug

% Energy Consumption
centr_energy = (summary_streams(2,5)*24/1000)*cent_en_in; % kWh per day
pump_electricity(1,4) = membr_to_centr*(summary_streams(2,5)/1000); % kW

% % UV Sterilization - only membrane and centrifugation is passed through sterilizer, settler
% stream is assumed to be already clarified
water_rec = sum(recycle_streams(2:3,1))*rec_eff + recycle_streams(1,1); % kg per hour

% Electricity to sterilize and pump water from membrane and centrifuge back to ponds
uv_energy = (sum(recycle_streams(2:3,1))*24/1000)*UV_en_in; % kWh per day
pump_electricity(1,5) = pump_membr_centr_rec*sum(recycle_streams(2:3,1))/1000; % kW

% Global Water Balance
water_left = water_vol - harvest_water*marker/1000; % water left in ponds after harvest m3
water_req = (volume_o - water_left) - water_rec*marker/1000; % water required to reach inoculation volume m3

if water_req > P_m3 % if future rain (P) is not enough to reach inoc volume, add water, else no make up water is required
    aquifer_makeup = water_req - P_m3; % m3
else
    aquifer_makeup = 0;% m3
end

% Volume after refill
% aquifer + recycled + left in pond after harvest
vol_refill = water_left + (water_rec*marker/1000) + aquifer_makeup; % m3

% Make sure new volume is at least min depth, if not water demand from
% aquifer increases
if (vol_refill*100/area) < min_depth % cm
    aquifer_makeup = aquifer_makeup + ((min_depth*area/100)-vol_refill); % m3
    vol_refill = aquifer_makeup + (water_rec*marker/1000) + water_left; % m3
end

% Calculate new AFDW concentration
algae_left = CX*water_left; % g AFDW algae
algae_added = sum(summary_streams(4:5,2))*marker*1000; % g AFDW biomass that is returned to ponds and not cent to conversion - used in net productivity calc
IBC = (algae_left + algae_added)/vol_refill; % (g per m3) initial concentration
CX_np1 = IBC; % g per m3 - must always equal IBC

% Calculate energy to pump freshwater into ponds
pump_electricity(1,6) = aquifer_pump*aquifer_makeup/marker; % kW

% Total pump power
pump_dew = sum(pump_electricity(1:5)); % kW for dewatering and return water
pump_makeup = pump_electricity(6);% kW for evap makeup

% Global salt mass balance
% Salt mass balance in pond includes what is transported in the
% recycled stream (salt_dewater_in) minus blowdown + what is introduced by the
% aquifer stream, and what is left in th pond
salt_mass = salt_cx*water_left  +  sum(summary_streams(4:5,4))*marker*1000 + aquifer_makeup*i_salt; % g, net salt mass 
salt_cx_np1 =salt_mass/vol_refill; % g per m3

% Outputs 
pond_states = [CX_np1 vol_refill salt_cx_np1 algae_added aquifer_makeup salt_mass]; 
energy_matrix = [settling_pond_no membrane_energy centr_energy uv_energy pump_dew pump_makeup]; % kWh per day for mem, cent, and UV; kW for pumping
end

% Nutrients mass balance
function [consumable_vec] = nutrient(CX_np1,CX,volume,comp_vector)
%Calculate the CO2 demand kg/hr based on stoichiometric carbon balance 
  
   mass_frac_C = comp_vector(1); 
   CUE = comp_vector(2); 
   Algae_P_Content = comp_vector(3); 
   DAP_P_Content = comp_vector(4); 
   DAP_N_Content = comp_vector(5); 
   Algae_N_Content = comp_vector(6); 
   Ammonia_N_Content = comp_vector(7); 
   Net_surplus = comp_vector(8);
     
   CO2_demand_hourly = (((CX_np1 -CX)*volume/1000)*mass_frac_C*(44/12)/CUE); % kg CO2 per hour
   DAP_demand_hrly = (((CX_np1- CX)*volume/1000)*Algae_P_Content/DAP_P_Content); % kg per hour
   DAP_demand_hrly = DAP_demand_hrly*(1+Net_surplus) ;% kg per hour
   Ammonia_demand_hrly = (((CX_np1-CX)*volume/1000)*Algae_N_Content/Ammonia_N_Content)-(DAP_demand_hrly*DAP_N_Content); % kg per hour
   Ammonia_demand_hrly = Ammonia_demand_hrly*(1+Net_surplus);% kg per hour

   consumable_vec = [DAP_demand_hrly Ammonia_demand_hrly CO2_demand_hourly]; % kg per hour
   
end

% Precipitation forecast
function [P]= p_future(i,area,rainfall,n2,harvest_days)
    if i+harvest_days*24-1 < n2
        P = sum(rainfall(i:i+harvest_days*24-1))*area; % m3 in time period
    else
        P = sum(rainfall(i:end))*area; % m3 
    end
end
