%% CAPEX and OPEX Function for Cultivation and Dewatering
function[cult_TEA_out_CAPEX, cult_TEA_out_OPEX] = Cultivation_TEA_fcn(Cost_year, num_ponds, pond_size, module_size, CO2_demand_hourly, aquifer_hrly, total_recycled, storage_fcn_out, num_settling_ponds, settlers_out, membrane_out, cult_out, cost_matrix, recycled_membr_centr, nut_rem)

wetted_acres = pond_size*num_ponds; 

load background_data.mat %#ok<LOAD> 
CAPEX_INDEX = table2array(CAPEX_CEPI_INDEX);
CAP_IND_POS = Cost_year - 1989; 

% Land CAPEX
cost_per_acre = 3000*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2)); %$3000/acre (2014 dollars) from Davis et al., 2016
land_CAPEX = num_ponds*pond_size*1.52*cost_per_acre; %1.52 scaling value based on 7600 acres total for 5000 wetted acres - Davis
total_land = num_ponds*pond_size*1.52;


% Pond CAPEX
if pond_size == 2.2 
    
    civil_work_CAPX = num_ponds*pond_size*1.3*11462.42*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2));
    leakage_bar_CAPX = num_ponds*1630.17*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2)); 
    paddlewheel_CAPX = num_ponds*38800*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2)); 
    elec_wiring_CAPX = num_ponds*9513.18*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2));
    elec_inst_CAPX = num_ponds*7900*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2)); 
    piping_CAPX = num_ponds*30612.94*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2));
    CO2_diff_CAPX = num_ponds*500*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2)); 
    slide_gate_CAPX = num_ponds*5851*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2)); 
    concrete_CAPX = num_ponds*16770.32*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2)); 
    
elseif pond_size == 10
    
    civil_work_CAPX = num_ponds*pond_size*1.3*10907.17*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2));
    leakage_bar_CAPX = num_ponds*4368.65*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2)); 
    paddlewheel_CAPX = num_ponds*64000*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2)); 
    elec_wiring_CAPX = num_ponds*13701.95*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2));
    elec_inst_CAPX = num_ponds*7900*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2)); 
    piping_CAPX = num_ponds*82369.50*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2));
    CO2_diff_CAPX = num_ponds*1000*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2)); 
    slide_gate_CAPX = num_ponds*5851*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2)); 
    concrete_CAPX = num_ponds*32878.80*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2));    
    
end
total_pond_CAPEX = civil_work_CAPX + leakage_bar_CAPX + paddlewheel_CAPX + elec_wiring_CAPX + elec_inst_CAPX + piping_CAPX + CO2_diff_CAPX + slide_gate_CAPX + concrete_CAPX; 

% Inoculation System CAPEX
tub_PBR_CAPX = wetted_acres/5000*35242*18.22*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2));
cov_IP_CAPX = wetted_acres/5000*23/100*5427316*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2)); 
open_IP_CAPX = wetted_acres/5000*116/100*5427316*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2)); 
cov_IP_lining_CAPX = wetted_acres/5000*1106775/4765888*3097827*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2)); 
air_sup_GH_CAPX = wetted_acres/5000*1112756*3*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2)); 
open_IP_lining_CAPX = wetted_acres/5000*5533873/4765888*3097827*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2)); 

total_inoc_CAPEX = tub_PBR_CAPX + cov_IP_CAPX + open_IP_CAPX + cov_IP_lining_CAPX + air_sup_GH_CAPX + open_IP_lining_CAPX; 

% CO2 Delivery System CAPEX based on max summer demand (kg/day) delivered
% over 12 hours

start_day = 1; 
daily_CO2_kg = zeros(365,1);

for i = 1:365

    daily_CO2_kg(i,1) = sum(CO2_demand_hourly(start_day:start_day + 23, 1)); 
    start_day = start_day + 24;

end 

max_CO2_kg_hr = max(daily_CO2_kg)/12; 

% CO2 storage sphere
CO2_sphere_cost = 1400800; %2014 dollars
CO2_sphere_scaling_value = 68550; %kg CO2/hr
CO2_sphere_scaling_exp = 0.6; 
CO2_sphere_install = 1.25; 
CO2_sphere_new_value = max_CO2_kg_hr; 
CO2_sphere_size_ratio = CO2_sphere_new_value/CO2_sphere_scaling_value; 
CO2_sphere_scaled_purch_cost = CO2_sphere_size_ratio^CO2_sphere_scaling_exp*CO2_sphere_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2));
CO2_sphere_install_cost = CO2_sphere_scaled_purch_cost*CO2_sphere_install;

% CO2 storage tank immersion vaporizeres
CO2_vap_cost = 70500; %2014 dollars
CO2_vap_scaling_value = 68550; %kg CO2/hr
CO2_vap_scaling_exp = 1; 
CO2_vap_install = 1.76; 
CO2_vap_new_value = max_CO2_kg_hr; 
CO2_vap_size_ratio = CO2_vap_new_value/CO2_vap_scaling_value; 
CO2_vap_scaled_purch_cost = CO2_vap_size_ratio^CO2_vap_scaling_exp*CO2_vap_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2));
CO2_vap_install_cost = CO2_vap_scaled_purch_cost*CO2_vap_install;

% Other CO2 equipment
scaling_value = 68550; 
CO2_trunkline = max_CO2_kg_hr/scaling_value*1661900*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2)); 
CO2_branchline = max_CO2_kg_hr/scaling_value*912300*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2)); 
CO2_within_plot_piping = max_CO2_kg_hr/scaling_value*44200*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2)); 
CO2_supply_to_inoc = max_CO2_kg_hr/scaling_value*59300*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2)); 

CO2_delivery_CAPEX = CO2_sphere_install_cost + CO2_vap_install_cost + CO2_trunkline + CO2_branchline + CO2_within_plot_piping + CO2_supply_to_inoc;

% Makeup water delivery and onsite circulation equipment
channels = wetted_acres*100*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2));
pipe_34 = wetted_acres/5000*5883/5883*872158*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(23,2))*1.51; 
pipe_26 = wetted_acres/5000*9426/9426*852949*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(23,2))*1.48;
pipe_20 = wetted_acres/5000*9990/9990*469386*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(23,2))*1.72;
pipe_16 = wetted_acres/5000*9360/9360*416554*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(23,2))*1.58;
pipe_14 = wetted_acres/5000*9360/9360*321593*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(23,2))*1.63;
pipe_12 = wetted_acres/5000*9360/9360*267709*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(23,2))*1.67;

makeup_watr_flow_vel = 5; 
makeup_watr_pipeline_diam = 21; %inches

while makeup_watr_flow_vel > 4 %while loop to size makeup water pipeline

% Makeup water pipeline
makeup_watr_pipeline_diam = makeup_watr_pipeline_diam + 1; 
makeup_watr_pipeline_radius_m = (makeup_watr_pipeline_diam/2)/39.37; 
makeup_watr_pipeline_exp = 1.4; 
makeup_watr_pipeline_scaling_value = 36; %inch diameter
makeup_watr_pipeline_install_cost = (makeup_watr_pipeline_diam/makeup_watr_pipeline_scaling_value)^makeup_watr_pipeline_exp*2638100*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(20,2));

% Makeup water pump
makeup_watr_pump_scaling = 10626; %MMgal/yr - evaporative loss
makeup_watr_pump_exp = 0.6; 
makeup_watr_pump_new_value = sum(aquifer_hrly)*264.2/1000000; %MMgal/yr of makeup water from aquifer
makeup_watr_pump_size_ratio = makeup_watr_pump_new_value/makeup_watr_pump_scaling;
makeup_watr_pump_install_cost = makeup_watr_pump_size_ratio^makeup_watr_pump_exp*705500*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(20,2));

% Check pipe flow velocity and adjust if above 4 m/s
makeup_watr_pipe_area = 3.14*makeup_watr_pipeline_radius_m^2; 
makeup_watr_flow_vel = (max(aquifer_hrly)/3600)/makeup_watr_pipe_area; 

end 

water_del_circ_CAPEX = channels + pipe_34 + pipe_26 + pipe_20 + pipe_16 + pipe_14 + pipe_12 + makeup_watr_pipeline_install_cost + makeup_watr_pump_install_cost; 

%% Dewatering Equipment CAPEX

watr_return_pipe = wetted_acres/5000*168199/4300*284880*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2));

% Pumps for return water (from settlers, membrane, and centrifuge)
ret_watr_pump_cost = 121905; %2012 dollars
ret_watr_pump_scaling = 6653; %GPM per module (100 acre module)
ret_watr_pump_exp = 0.8; 
ret_watr_pump_install = 1.15; 
ret_watr_pump_new_value = max(total_recycled)/1000*264.2/60*(module_size/wetted_acres); %GPM per module
ret_watr_pump_size_ratio = ret_watr_pump_new_value/ret_watr_pump_scaling;
ret_watr_pump_install_cost = ret_watr_pump_size_ratio^ret_watr_pump_exp*ret_watr_pump_cost*(wetted_acres/module_size)*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(23,2))*ret_watr_pump_install; 

% Settlig ponds (including those for long-term storage) -2014 dollars
num_settling_ponds = max(num_settling_ponds) + str2double(storage_fcn_out(8,2));
settling_ponds_install_cost = num_settling_ponds*34300*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2));

% Liners for anaerobic storage settling ponds
num_settler_w_liners = str2double(storage_fcn_out(8,2));
settling_pond_area = (1000/4.572)*10.764; %area in sqft
settling_pond_liner_installed_cost = settling_pond_area*num_settler_w_liners*0.65; % $0.65/sqft - Davis

% Stage 2 and Stage 3
start_day = 1; 
daily_settlers_out = zeros(365,1);
daily_membrane_out = zeros(365,1); 

for i = 1:365

    daily_settlers_out(i,1) = sum(settlers_out(start_day:start_day + 23, 1)); 
    daily_membrane_out(i,1) = sum(membrane_out(start_day:start_day + 23, 1)); 
    start_day = start_day + 24;

end 

max_membrane_flow = max(daily_settlers_out)/1000*264.2; %gallons per day
max_cent_flow = max(membrane_out)/1000; %m3/hr

% Membranes
membrane_size_ratio = max_membrane_flow/20000000;
membrane_install_cost = membrane_size_ratio*4760000*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2));
membrane_other_install_cost = membrane_size_ratio^0.6*8104000*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2));
membrane_CAPEX = membrane_install_cost + membrane_other_install_cost; 

% Centrifuge
cent_exp = 0.6; 
cent_cost = 2242500; %2013 dollars
cent_install = 1.8; 
cent_scaling = 463; %m3/hr
cent_size_ratio = max_cent_flow/cent_scaling; 
cent_install_cost = cent_size_ratio^cent_exp*cent_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(24,2))*cent_install;

% UV Sterilization
UV_cost = 20000; %2014 dollars
UV_scaling = 500; %GPM
UV_install = 1.15; 
UV_new_value = max(recycled_membr_centr)/1000*264.2/60; %GPM (membrane and centrifuge only)
UV_size_ratio = UV_new_value/UV_scaling;
UV_install_cost = UV_size_ratio*UV_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(25,2))*UV_install; 

total_dewatering_CAPEX = watr_return_pipe + ret_watr_pump_install_cost + settling_ponds_install_cost + settling_pond_liner_installed_cost + membrane_CAPEX + cent_install_cost + UV_install_cost; 

%% Storage Equipment CAPX
% Product storage tank accounted for with additional settling ponds
% (anaerobic long-term storage)

% Firewater pump - scaled by facility size
fire_watr_pump_cost = 15000; %2009 dollars
fire_watr_pump_exp = 0.8; 
fire_watr_pump_install = 3.1; 
fire_watr_pump_scaling = 8343; 
fire_watr_pump_new_value = wetted_acres/5000*11051; 
fire_watr_pump_size_ratio = fire_watr_pump_new_value/fire_watr_pump_scaling;
fire_watr_pump_install_cost = fire_watr_pump_size_ratio^fire_watr_pump_exp*fire_watr_pump_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(20,2))*fire_watr_pump_install; 

% Firewater storage tank - scaled by facility size
fire_watr_ST_cost = 803000; %2009 dollars
fire_watr_ST_exp = 0.7; 
fire_watr_ST_install = 1.7; 
fire_watr_ST_scaling = 8343; 
fire_watr_ST_new_value = wetted_acres/5000*11051; 
fire_watr_ST_size_ratio = fire_watr_ST_new_value/fire_watr_ST_scaling;
fire_watr_ST_install_cost = fire_watr_ST_size_ratio^fire_watr_ST_exp*fire_watr_ST_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(20,2))*fire_watr_ST_install; 

% Make-up storage tank
makeup_watr_ST_cost = 418795; %2009 dollars
makeup_watr_ST_scaling = 750000; %gallons
makeup_watr_ST_exp = 0.6; 
makeup_watr_ST_install = 1.14; 
makeup_watr_ST_new_value = mean(nonzeros(aquifer_hrly))*1000/0.993*6/3.78541/4; 
makeup_watr_ST_size_ratio = makeup_watr_ST_new_value/makeup_watr_ST_scaling;
makeup_watr_ST_install_cost = makeup_watr_ST_size_ratio^makeup_watr_ST_exp*makeup_watr_ST_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(20,2))*makeup_watr_ST_install*4; 

% Tankage BOP
tankage_BOP = (fire_watr_pump_install_cost + fire_watr_ST_install_cost + makeup_watr_ST_install_cost)*0.20; 

storage_tank_CAPEX = fire_watr_pump_install_cost + fire_watr_ST_install_cost + makeup_watr_ST_install_cost + tankage_BOP; 

%% Additional Direct CAPEX
Cultivation_CAPEX = total_pond_CAPEX + total_inoc_CAPEX;
Dewatering_CAPEX = total_dewatering_CAPEX;
OSBL_CAPEX = CO2_delivery_CAPEX + water_del_circ_CAPEX + storage_tank_CAPEX;

warehouse_CAPEX = 0.04*Dewatering_CAPEX + 0.012*Cultivation_CAPEX;
site_devel_CAPEX = 0.09*Dewatering_CAPEX + 1534*total_land; %1534/acre for roads/fences/etc.
add_piping_CAPEX = 0.045*Dewatering_CAPEX; 

total_add_direct_CAPEX = warehouse_CAPEX + site_devel_CAPEX + add_piping_CAPEX;

%% Additional Indirect CAPEX
prorat_exp_CAPEX = 0.1*Dewatering_CAPEX + 0.04*Cultivation_CAPEX + 0.01*OSBL_CAPEX; 
fiel_exp_CAPEX = 0.1*Dewatering_CAPEX + 0.045*Cultivation_CAPEX + 0.01*OSBL_CAPEX; 
home_off_CAPEX = 0.2*Dewatering_CAPEX + 0.103*Cultivation_CAPEX + 0.01*OSBL_CAPEX;
proj_cont_CAPEX = 0.1*Dewatering_CAPEX + 0.1*Cultivation_CAPEX + 0.10*OSBL_CAPEX;
other_costs_CAPEX = 0.1*Dewatering_CAPEX + 0.026*Cultivation_CAPEX + 0.01*OSBL_CAPEX;

total_add_indirect_CAPEX = prorat_exp_CAPEX + fiel_exp_CAPEX + home_off_CAPEX + proj_cont_CAPEX + other_costs_CAPEX;

%% Total Capital Investment

FCI_cultivation = Cultivation_CAPEX + Dewatering_CAPEX + OSBL_CAPEX + total_add_direct_CAPEX + total_add_indirect_CAPEX;
working_CAPEX = 0.05*FCI_cultivation; 

TCI_cultivation = FCI_cultivation + working_CAPEX + land_CAPEX; 

%% Operational Expenses

CHEM_INDEX = table2array(CHEMICAL_INDEX); 
%load('ELEC_INDEX.mat'); electricity is indexed on the model interface page
%ELECT_INDEX = table2array(ELEC_INDEX);
LABOR_INDEX = table2array(LABOR_INDEX); %#ok<NODEF> 

CHEM_IND_POS = Cost_year - 1979; 
%ELEC_IND_POS = Cost_year - 2013; indexed on interface page
LABOR_IND_POS = Cost_year - 1989; 

%% Variable OPEX
% Energy
% Cultivation energy
pumping_en = str2double(cult_out(11,2)) + str2double(cult_out(17,2));
paddlewheel_en = str2double(cult_out(12,2));
sparging_en = str2double(cult_out(13,2));
membrane_en = str2double(cult_out(14,2));
centrifuge_en = str2double(cult_out(15,2));
UV_ster_en = str2double(cult_out(16,2));

total_cult_energy = pumping_en + paddlewheel_en + sparging_en + membrane_en + centrifuge_en + UV_ster_en; 
cost_of_energy = cost_matrix(1,1); %2019 dollars $/kWh
cult_energy_OPX = total_cult_energy*cost_of_energy; 

% Cultivation nutrients
ammonia_kg_yr = str2double(cult_out(8,2))*(1-nut_rem);
cost_of_ammonia = cost_matrix(1,2); % $/kg in 2019 dollars
ammonia_OPX = ammonia_kg_yr*cost_of_ammonia*(CHEM_INDEX(CHEM_IND_POS,2)/CHEM_INDEX(40,2));

dap_kg_yr = str2double(cult_out(9,2))*(1-nut_rem);
cost_of_dap = cost_matrix(1,3); % $/kg in 2019 dollars
dap_OPX = dap_kg_yr*cost_of_dap*(CHEM_INDEX(CHEM_IND_POS,2)/CHEM_INDEX(40,2)); 

% Cultivation CO2
annual_CO2 = str2double(cult_out(10,2)); %tonnes CO2 per year
cost_of_CO2 = cost_matrix(1,4); % 2019 dollars/tonne CO2 - Davis et al., 2016
CO2_OPX = annual_CO2*cost_of_CO2*(CHEM_INDEX(CHEM_IND_POS,2)/CHEM_INDEX(40,2)); 

%% Fixed OPEX (labor, maintenance, insurance, taxes)
%Labor - all salaries in 2018 USD
cult_plant_mngr = 1 * 155400 * (LABOR_INDEX(LABOR_IND_POS,2)/LABOR_INDEX(22,2));
cult_plant_civil_engr = 2 * 81935 * (LABOR_INDEX(LABOR_IND_POS,2)/LABOR_INDEX(22,2));
cult_plant_env_engr = 2 * 83244 * (LABOR_INDEX(LABOR_IND_POS,2)/LABOR_INDEX(22,2));
cult_maint_sup = 1 * 60257 * (LABOR_INDEX(LABOR_IND_POS,2)/LABOR_INDEX(22,2));
cult_maint_tech = 12 * 42286 * (LABOR_INDEX(LABOR_IND_POS,2)/LABOR_INDEX(22,2));
cult_lab_mngr = 1 * 59200 * (LABOR_INDEX(LABOR_IND_POS,2)/LABOR_INDEX(22,2));
cult_lab_tech = 1 * 42286 * (LABOR_INDEX(LABOR_IND_POS,2)/LABOR_INDEX(22,2));
cult_shft_sup = 4 * 50743 * (LABOR_INDEX(LABOR_IND_POS,2)/LABOR_INDEX(22,2));
cult_op_prod = 56 * 26872 * (LABOR_INDEX(LABOR_IND_POS,2)/LABOR_INDEX(22,2));
cult_op_inoc = 8 * 44038 * (LABOR_INDEX(LABOR_IND_POS,2)/LABOR_INDEX(22,2));
cult_op_dewater = 9 * 38536 * (LABOR_INDEX(LABOR_IND_POS,2)/LABOR_INDEX(22,2));
cult_clk_sec = 3 * 38057 * (LABOR_INDEX(LABOR_IND_POS,2)/LABOR_INDEX(22,2));

labor_total = (wetted_acres/5000)*(cult_plant_mngr + cult_plant_civil_engr + cult_plant_env_engr + cult_maint_sup + cult_maint_tech + cult_lab_mngr + cult_lab_tech + cult_shft_sup + cult_op_prod + cult_op_inoc + cult_op_dewater + cult_clk_sec);
labor_burden = labor_total * 0.90; 
maintenance = 0.03 * Dewatering_CAPEX + 0.005*(Cultivation_CAPEX + CO2_delivery_CAPEX + water_del_circ_CAPEX); 
prop_ins_tax = 0.007 * FCI_cultivation; 

%% Outputs Matrix

% CAPEX
cult_TEA_out_CAPEX(1,1) = "Civil Work CAPEX";
cult_TEA_out_CAPEX(1,2) = civil_work_CAPX;
cult_TEA_out_CAPEX(2,1) = "Leakage Barrier CAPEX";
cult_TEA_out_CAPEX(2,2) = leakage_bar_CAPX;
cult_TEA_out_CAPEX(3,1) = "Paddlewheel CAPEX";
cult_TEA_out_CAPEX(3,2) = paddlewheel_CAPX;
cult_TEA_out_CAPEX(4,1) = "Electrical Wiring CAPEX";
cult_TEA_out_CAPEX(4,2) = elec_wiring_CAPX;
cult_TEA_out_CAPEX(4,1) = "Electrical Instrumentation CAPEX";
cult_TEA_out_CAPEX(4,2) = elec_inst_CAPX;
cult_TEA_out_CAPEX(5,1) = "Pond Piping CAPEX";
cult_TEA_out_CAPEX(5,2) = piping_CAPX;
cult_TEA_out_CAPEX(6,1) = "CO2 Diffusers CAPEX";
cult_TEA_out_CAPEX(6,2) = CO2_diff_CAPX;
cult_TEA_out_CAPEX(7,1) = "Slide Gates CAPEX";
cult_TEA_out_CAPEX(7,2) = slide_gate_CAPX;
cult_TEA_out_CAPEX(8,1) = "Concrete CAPEX";
cult_TEA_out_CAPEX(8,2) = concrete_CAPX;
cult_TEA_out_CAPEX(9,1) = "Total ORP CAPEX";
cult_TEA_out_CAPEX(9,2) = total_pond_CAPEX;
cult_TEA_out_CAPEX(10,1) = "Tubular PBR (Inoculation System) CAPEX";
cult_TEA_out_CAPEX(10,2) = tub_PBR_CAPX;
cult_TEA_out_CAPEX(11,1) = "Covered Inoculum Ponds CAPEX";
cult_TEA_out_CAPEX(11,2) = cov_IP_CAPX;
cult_TEA_out_CAPEX(12,1) = "Open Inoculum Ponds CAPEX";
cult_TEA_out_CAPEX(12,2) = open_IP_CAPX;
cult_TEA_out_CAPEX(13,1) = "Covered Inoculum Pond Lining CAPEX";
cult_TEA_out_CAPEX(13,2) = cov_IP_lining_CAPX;
cult_TEA_out_CAPEX(14,1) = "Air Supported Greenhouse CAPEX";
cult_TEA_out_CAPEX(14,2) = air_sup_GH_CAPX;
cult_TEA_out_CAPEX(15,1) = "Open Inoculum Pond Lining CAPEX";
cult_TEA_out_CAPEX(15,2) = open_IP_lining_CAPX;
cult_TEA_out_CAPEX(16,1) = "Total Inoculation System CAPEX";
cult_TEA_out_CAPEX(16,2) = total_inoc_CAPEX;
cult_TEA_out_CAPEX(17,1) = "CO2 Sphere CAPEX";
cult_TEA_out_CAPEX(17,2) = CO2_sphere_install_cost;
cult_TEA_out_CAPEX(18,1) = "CO2 Vaporizer CAPEX";
cult_TEA_out_CAPEX(18,2) = CO2_vap_install_cost;
cult_TEA_out_CAPEX(19,1) = "CO2 Trunkline CAPEX";
cult_TEA_out_CAPEX(19,2) = CO2_trunkline;
cult_TEA_out_CAPEX(20,1) = "CO2 Branchline CAPEX";
cult_TEA_out_CAPEX(20,2) = CO2_branchline;
cult_TEA_out_CAPEX(21,1) = "CO2 within-plot Piping CAPEX";
cult_TEA_out_CAPEX(21,2) = CO2_within_plot_piping;
cult_TEA_out_CAPEX(22,1) = "CO2 Supply to Inoculum Ponds CAPEX";
cult_TEA_out_CAPEX(22,2) = CO2_supply_to_inoc;
cult_TEA_out_CAPEX(23,1) = "Total CO2 Delivery CAPEX";
cult_TEA_out_CAPEX(23,2) = CO2_delivery_CAPEX;
cult_TEA_out_CAPEX(24,1) = "Harvesting Chanels CAPEX";
cult_TEA_out_CAPEX(24,2) = channels;
cult_TEA_out_CAPEX(25,1) = "34 Inch Diameter Pipe CAPEX";
cult_TEA_out_CAPEX(25,2) = pipe_34;
cult_TEA_out_CAPEX(26,1) = "26 Inch Diameter Pipe CAPEX";
cult_TEA_out_CAPEX(26,2) = pipe_26;
cult_TEA_out_CAPEX(27,1) = "20 Inch Diameter Pipe CAPEX";
cult_TEA_out_CAPEX(27,2) = pipe_20;
cult_TEA_out_CAPEX(28,1) = "16 Inch Diameter Pipe CAPEX";
cult_TEA_out_CAPEX(28,2) = pipe_16;
cult_TEA_out_CAPEX(29,1) = "14 Inch Diameter Pipe CAPEX";
cult_TEA_out_CAPEX(29,2) = pipe_14;
cult_TEA_out_CAPEX(30,1) = "12 Inch Diameter Pipe CAPEX";
cult_TEA_out_CAPEX(30,2) = pipe_12;
cult_TEA_out_CAPEX(31,1) = "Makeup Water Pipeline CAPEX";
cult_TEA_out_CAPEX(31,2) = makeup_watr_pipeline_install_cost;
cult_TEA_out_CAPEX(32,1) = "Makeup Water Pipeline Diameter";
cult_TEA_out_CAPEX(32,2) = makeup_watr_pipeline_diam;
cult_TEA_out_CAPEX(33,1) = "Makeup Water Pump CAPEX";
cult_TEA_out_CAPEX(33,2) = makeup_watr_pump_install_cost;
cult_TEA_out_CAPEX(34,1) = "Water Delivery and Circulation CAPEX";
cult_TEA_out_CAPEX(34,2) = water_del_circ_CAPEX;
cult_TEA_out_CAPEX(35,1) = "Return Water Pipe CAPEX";
cult_TEA_out_CAPEX(35,2) = watr_return_pipe;
cult_TEA_out_CAPEX(36,1) = "Return Water Pump CAPEX";
cult_TEA_out_CAPEX(36,2) = ret_watr_pump_install_cost;
cult_TEA_out_CAPEX(37,1) = "UV Sterilizer CAPEX";
cult_TEA_out_CAPEX(37,2) = UV_install_cost;
cult_TEA_out_CAPEX(38,1) = "Settling Ponds CAPEX";
cult_TEA_out_CAPEX(38,2) = settling_ponds_install_cost;
cult_TEA_out_CAPEX(39,1) = "Settling Pond Liners (Anaerobic Storage) CAPEX";
cult_TEA_out_CAPEX(39,2) = settling_pond_liner_installed_cost;
cult_TEA_out_CAPEX(40,1) = "Membrane CAPEX";
cult_TEA_out_CAPEX(40,2) = membrane_CAPEX;
cult_TEA_out_CAPEX(41,1) = "Centrifuge CAPEX";
cult_TEA_out_CAPEX(41,2) = cent_install_cost;
cult_TEA_out_CAPEX(42,1) = "Total Dewatering CAPEX";
cult_TEA_out_CAPEX(42,2) = total_dewatering_CAPEX;
cult_TEA_out_CAPEX(43,1) = "Firewater Pump CAPEX";
cult_TEA_out_CAPEX(43,2) = fire_watr_pump_install_cost;
cult_TEA_out_CAPEX(44,1) = "Firewater Storage Tank CAPEX";
cult_TEA_out_CAPEX(44,2) = fire_watr_ST_install_cost;
cult_TEA_out_CAPEX(45,1) = "Makeup Water Storage Tank CAPEX";
cult_TEA_out_CAPEX(45,2) = makeup_watr_ST_install_cost;
cult_TEA_out_CAPEX(46,1) = "Tankage BOP CAPEX";
cult_TEA_out_CAPEX(46,2) = tankage_BOP;
cult_TEA_out_CAPEX(47,1) = "Total Storage CAPEX";
cult_TEA_out_CAPEX(47,2) = storage_tank_CAPEX;
cult_TEA_out_CAPEX(48,1) = "Warehouse CAPEX";
cult_TEA_out_CAPEX(48,2) = warehouse_CAPEX;
cult_TEA_out_CAPEX(49,1) = "Site Development CAPEX";
cult_TEA_out_CAPEX(49,2) = site_devel_CAPEX;
cult_TEA_out_CAPEX(50,1) = "Additional Piping CAPEX";
cult_TEA_out_CAPEX(50,2) = add_piping_CAPEX;
cult_TEA_out_CAPEX(51,1) = "Total Additional Direct CAPEX";
cult_TEA_out_CAPEX(51,2) = total_add_direct_CAPEX;
cult_TEA_out_CAPEX(52,1) = "Prorateable Expenses CAPEX";
cult_TEA_out_CAPEX(52,2) = prorat_exp_CAPEX;
cult_TEA_out_CAPEX(53,1) = "Field Expenses CAPEX";
cult_TEA_out_CAPEX(53,2) = fiel_exp_CAPEX;
cult_TEA_out_CAPEX(54,1) = "Home Office CAPEX";
cult_TEA_out_CAPEX(54,2) = home_off_CAPEX;
cult_TEA_out_CAPEX(55,1) = "Project Contingency CAPEX";
cult_TEA_out_CAPEX(55,2) = proj_cont_CAPEX;
cult_TEA_out_CAPEX(56,1) = "Other Costs CAPEX";
cult_TEA_out_CAPEX(56,2) = other_costs_CAPEX;
cult_TEA_out_CAPEX(57,1) = "Total Additional Indirect CAPEX";
cult_TEA_out_CAPEX(57,2) = total_add_indirect_CAPEX;
cult_TEA_out_CAPEX(58,1) = "Fixed Capital Investment";
cult_TEA_out_CAPEX(58,2) = FCI_cultivation;
cult_TEA_out_CAPEX(59,1) = "Land CAPEX";
cult_TEA_out_CAPEX(59,2) = land_CAPEX;
cult_TEA_out_CAPEX(60,1) = "Working Capital";
cult_TEA_out_CAPEX(60,2) = working_CAPEX;
cult_TEA_out_CAPEX(61,1) = "Total Capital Investment";
cult_TEA_out_CAPEX(61,2) = TCI_cultivation;

% OPEX
cult_TEA_out_OPEX(1,1) = "Cultivation and Dewatering Energy OPEX";
cult_TEA_out_OPEX(1,2) = cult_energy_OPX;
cult_TEA_out_OPEX(2,1) = "Ammonia OPEX";
cult_TEA_out_OPEX(2,2) = ammonia_OPX;
cult_TEA_out_OPEX(3,1) = "DAP OPEX";
cult_TEA_out_OPEX(3,2) = dap_OPX;
cult_TEA_out_OPEX(4,1) = "CO2 OPEX";
cult_TEA_out_OPEX(4,2) = CO2_OPX;
cult_TEA_out_OPEX(5,1) = "Total Labor OPEX";
cult_TEA_out_OPEX(5,2) = labor_total;
cult_TEA_out_OPEX(6,1) = "Total Burden OPEX";
cult_TEA_out_OPEX(6,2) = labor_burden;
cult_TEA_out_OPEX(7,1) = "Total Maintenance OPEX";
cult_TEA_out_OPEX(7,2) = maintenance;
cult_TEA_out_OPEX(8,1) = "Total Property Tax and Insurance OPEX";
cult_TEA_out_OPEX(8,2) = prop_ins_tax;
cult_TEA_out_OPEX(9,1) = "Total Cultivation and Dewatering OPEX";
cult_TEA_out_OPEX(9,2) = cult_energy_OPX + ammonia_OPX + dap_OPX + CO2_OPX + labor_total + labor_burden + maintenance + prop_ins_tax;

end
