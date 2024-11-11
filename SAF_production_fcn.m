function [LE_CAPEX_out, LE_OPEX_out, SAF_CAPEX_out, SAF_OPEX_out, SAF_proc_mod_out, LE_proc_mod_out, DCFROR_out, LE_coprod_mat, SAF_total_coprod_mat, SAF_LCA_consum] = SAF_production_fcn(comp_vector, cult_TEA_out_CAPEX, cult_TEA_out_OPEX, wetted_acres, cult_out, op_hours, ann_ave_centr, HPH_num_passes, HPH_flowrate,  SAF_consumables_costs, WHE_model_inputs, algae_comp, SAF_fuels_value_mat, Cost_year, IRR, SAF_dropdowns, ng_boiler_eff, nut_rem, protein_revenue)
%This function estimates energy consumption and costs associated with the
%pretreatment steps and HEFA process for the SAF pathway. The pretreatment steps include
%high pressure homogenization and hexane solvent extraction for maximum lipid
%recovery. Once lipids are extracted the algae-oil (triglycerides) are sent
%to the HEFA process to produce SAF, diesel, LPG, and naphtha. The lipid
%extracted algae (LEA) is sent to AD for nutrient and energy recovery. Net
%energy input is calculated based on the delta between total energy demand
%and total recovered energy from AD/CHP.

load background_data.mat %#ok<LOAD> 

% Cost Indices and establishing project year
CAPEX_INDEX = table2array(CAPEX_CEPI_INDEX);
CHEM_INDEX = table2array(CHEMICAL_INDEX); 
LABOR_INDEX = table2array(LABOR_INDEX); %#ok<NODEF> 

CAP_IND_POS = Cost_year - 1989; 
CHEM_IND_POS = Cost_year - 1979; 
LABOR_IND_POS = Cost_year - 1989; 

%Commodity prices
Elec_cost = SAF_consumables_costs(1,1); 
ammonia_cost = SAF_consumables_costs(2,1); 
dap_cost = SAF_consumables_costs(3,1); 
CO2_cost = SAF_consumables_costs(4,1); 
NG_cost = SAF_consumables_costs(5,1); 
% NG_transport = SAF_consumables_costs(6,1); 
H_cost = SAF_consumables_costs(7,1); 
Process_water_cost = SAF_consumables_costs(8,1); 
Hexane_cost = SAF_consumables_costs(9,1); 
carbon_offset_cost = SAF_consumables_costs(10,1); 
lip_purch_price = SAF_consumables_costs(11,1); 

value_propane_gal = SAF_fuels_value_mat(1,1);  
value_LNG_gal = SAF_fuels_value_mat(2,1);  
value_naphtha_gal = SAF_fuels_value_mat(3,1);  
value_jet_gal = SAF_fuels_value_mat(4,1);  
value_diesel_gal = SAF_fuels_value_mat(5,1);  

prot = algae_comp(1,1);
carb = algae_comp(2,1); 
lip = algae_comp(3,1);
ash = algae_comp(4,1); 
TSS = algae_comp(5,1); 

n_content = comp_vector(6); 
p_content = comp_vector(3);

LE_inc_AD_CHP = SAF_dropdowns(1,1);
HEFA_scenario = SAF_dropdowns(2,1);
HEFA_H2_scenario = SAF_dropdowns(3,1);
HEFA_CO2_rec = SAF_dropdowns(4,1);
fixed_opex_case = SAF_dropdowns(5,1);

%----------------------------------------------------------------------------------------------------------------------------------------------------

% SAF Production Model 
% HPH Energy consumption (Kang et al., 2019)
HPH_feed = ann_ave_centr; %annual average hourly flow from storage/ponds to HPH (kg/hr 20% solids)
HPH_spec_energy = 0.1784; % kWh_electric/kg dry biomass per unit
%HPH_disrup_eff = 0.90; % 90% of cells disrupted 
HPH_ann_energy_kWh = ann_ave_centr*TSS*op_hours*HPH_spec_energy*HPH_num_passes; % kWh_electric per year

% HPH capital costs (Kang et al., 2019)
HPH_scaling = HPH_flowrate; 
HPH_new_value = HPH_feed/1000; %m3/hour (average annual feed to HPH)
HPH_cost = 475209; % 2016 dollars (Kang et al., 2019)
HPH_install_factor = 1.0; % estimated
HPH_scaling_exp = 0.7; % estimated
HPH_scaled_install_cost = HPH_cost*(HPH_new_value/HPH_scaling)^HPH_scaling_exp*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(27,2))*HPH_install_factor*HPH_num_passes; 

% HPH operational costs
HPH_energy_opex = HPH_ann_energy_kWh*Elec_cost; %#ok<NASGU> 

% WHE Inputs
WHE_thermal_en = WHE_model_inputs(1,1);
WHE_electric_en = WHE_model_inputs(2,1);
WHE_solvent_to_biomass = WHE_model_inputs(3,1);
WHE_hexane_loss = WHE_model_inputs(4,1);
WHE_extraction_eff = WHE_model_inputs(5,1);

% Wet Hexane Extraction (WHE) energy consumption (Argonne National Lab, 2011) - Currently most robust (transparent) source (2/20/23)
WHE_oil_out = ann_ave_centr*TSS*lip*WHE_extraction_eff; % kg oil/hr (lipid) exiting WHE
WHE_ann_thermal_kWh = WHE_oil_out*op_hours*WHE_thermal_en; % kWh_thermal/yr
WHE_ann_elec_kWh = WHE_oil_out*op_hours*WHE_electric_en; % kWh_elec/yr
WHE_hexane_input = ann_ave_centr*WHE_solvent_to_biomass; % kg/hr hexane required to process hourly flowrate of wet biomass
WHE_ann_hexane_loss = WHE_oil_out*op_hours*WHE_hexane_loss/1000; % kg hexane loss per year

% Wet Hexane Extraction (WHE) energy consumption (Shi et al., 2019)
% WHE_thermal_en = 3.08; % kWh_thermal/kg oil
% WHE_electric_en = 0.61; % kWh_elec/kg oil
% WHE_hexane_loss = 0.055; % kg hexane/kg oil
% WHE_extraction_eff = 0.95; % percent of lipids extracted
% WHE_solvent_to_biomass = 1.8; % 1.8 g hexane/g wet biomass
% WHE_oil_out = ann_ave_centr*0.20*lip*WHE_extraction_eff; % kg oil/hr (lipid) exiting WHE
% WHE_ann_thermal_kWh = WHE_oil_out*op_hours*WHE_thermal_en; % kWh_thermal/yr
% WHE_ann_elec_kWh = WHE_oil_out*op_hours*WHE_electric_en; % kWh_elec/yr
% WHE_hexane_input = ann_ave_centr*WHE_solvent_to_biomass; % kg hexane required to process hourly flowrate of wet biomass
% WHE_ann_hexane_loss = WHE_oil_out*op_hours*WHE_hexane_loss/1000; % kg hexane loss per year

% Wet Hexane Extraction (Kang et al., 2019)
% Inputs
% WHE_thermal_en = 1.3; % kWh_thermal/kg dry biomass
% WHE_electric_en = 0.00035; % kWh_elec/kg dry biomass
% WHE_hexane_loss = 0.01; % 1% loss of hexane input 99% recovery efficiency
% WHE_extraction_eff = 0.80; % percent of lipids extracted
% WHE_solvent_to_biomass = 10; % 10 kg hexane/kg dry biomass
% WHE_heat_eff = 0.8; % 80% heat efficiency
% % Outputs
% WHE_oil_out = ann_ave_centr*0.20*lip*WHE_extraction_eff; % kg oil/hr (lipid) exiting WHE
% WHE_ann_thermal_kWh = ann_ave_centr*op_hours*0.20*WHE_thermal_en; % kWh_thermal/yr
% WHE_ann_elec_kWh = ann_ave_centr*op_hours*0.20*WHE_electric_en; % kWh_elec/yr
% WHE_hexane_input = ann_ave_centr*op_hours*0.20*WHE_solvent_to_biomass; 
% WHE_ann_hexane_loss = WHE_hexane_input*WHE_hexane_loss; % kg hexane loss per year

% WHE capital costs (Delrue et al., 2012)
WHE_cost = 84.51; % 2011 $/tonne dry weight biomass per year
WHE_install_cost = WHE_cost*(ann_ave_centr*TSS*op_hours/1000)*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(22,2));
WHE_initial_hexane_charge = WHE_hexane_input*Hexane_cost; % kg/hr hexane * $/kg hexane

% WHE operational costs
WHE_heating_eff = 0.80; % assumption
WHE_NG_input = ((WHE_ann_thermal_kWh/WHE_heating_eff)*3412/983)*0.022; % kg NG - 3412 BTU/kWh; LHV = 983 BTU/SCF; 0.022 kg/SCF
WHE_elec_opex = WHE_ann_elec_kWh*Elec_cost; %#ok<NASGU> 
WHE_NG_opex = WHE_NG_input*NG_cost;
WHE_hexane_makeup = WHE_ann_hexane_loss*Hexane_cost; 

% WHE Outputs
WHE_outputs = [WHE_oil_out; WHE_ann_thermal_kWh; WHE_ann_elec_kWh; WHE_hexane_input; WHE_ann_hexane_loss]; %#ok<NASGU> 

%----------------------------------------------------------------------------------------------------------------------------------------------------

% HEFA Process
% HEFA product yields and hydrogen input (Pearlson et al., 2013)
HEFA_oil_in = WHE_oil_out; % kg oil per hour from hexane extraction
HEFA_oil_BPD = WHE_oil_out/0.916/3.785*24/42; % BPD of algae oil (approximated density as soybean oil of 0.916 g/mL); and 42 gallons per barrel
HEFA_oil_in_lbs = HEFA_oil_in*2.205; % lbs/hr

    switch HEFA_scenario
        case 'Maximum Distillate'
            HEFA_hydrogen_in = HEFA_oil_in/100*2.7; 
            HEFA_water_out = HEFA_oil_in/100*8.7;
            HEFA_CO2_out = HEFA_oil_in/100*5.5; 
            HEFA_propane_out = HEFA_oil_in/100*4.2; 
            HEFA_LNG_out = HEFA_oil_in/100*1.6; 
            HEFA_naphtha_out = HEFA_oil_in/100*1.8; 
            HEFA_jet_out = HEFA_oil_in/100*12.8; 
            HEFA_diesel_out = HEFA_oil_in/100*68.1;
        case 'Maximum Jet'
            HEFA_hydrogen_in = HEFA_oil_in/100*4; 
            HEFA_water_out = HEFA_oil_in/100*8.7; 
            HEFA_CO2_out = HEFA_oil_in/100*5.4; 
            HEFA_propane_out = HEFA_oil_in/100*4.2; 
            HEFA_LNG_out = HEFA_oil_in/100*6; 
            HEFA_naphtha_out = HEFA_oil_in/100*7; 
            HEFA_jet_out = HEFA_oil_in/100*49.4; 
            HEFA_diesel_out = HEFA_oil_in/100*23.3;
    end 

% HEFA total liquid fuel output in gallons per day
liq_propane_density = 0.493; % kg/L - density of liquid propane
LNG_density = 0.55; % kg/L - ranges from 0.525 to 0.580 kg/L (for pressurized liquid LPG)
diesel_density = 0.8366; % kg/L - Chen & Quinn 2021
jet_density = 0.8000; % 0.78 - 0.84 kg/L 
naphtha_density = 0.7447; % kg/L - Chen & Quinn 2021

HEFA_propane_out_gal_day = HEFA_propane_out/liq_propane_density/3.785*24; % gallons/day
HEFA_LNG_gal_day = HEFA_LNG_out/LNG_density/3.785*24; % gallons/day
HEFA_diesel_out_gal_day = HEFA_diesel_out/diesel_density/3.785*24; % gallons/day
HEFA_jet_out_gal_day = HEFA_jet_out/jet_density/3.785*24; % gallons/day
HEFA_naphtha_out_gal_day = HEFA_naphtha_out/naphtha_density/3.785*24; % gallons/day

HEFA_propane_out_gal_yr = HEFA_propane_out/liq_propane_density/3.785*op_hours; % gallons/yr
HEFA_LNG_out_gal_yr = HEFA_LNG_out/LNG_density/3.785*op_hours; % gallons/yr
HEFA_diesel_out_gal_yr = HEFA_diesel_out/diesel_density/3.785*op_hours; % gallons/yr
HEFA_jet_out_gal_yr = HEFA_jet_out/jet_density/3.785*op_hours; % gallons/yr
HEFA_naphtha_out_gal_yr = HEFA_naphtha_out/naphtha_density/3.785*op_hours; % gallons/yr

total_liquid_fuels = HEFA_diesel_out_gal_yr + HEFA_jet_out_gal_yr + HEFA_naphtha_out_gal_yr + HEFA_LNG_out_gal_yr; %#ok<NASGU> % total gallons fuel per year

%Total Energy in GGE
% Fuel energy density mmBTU per gal
% [renewable diesel, jet, naphtha, LNG, propane]
% Retrieved from GREET
fuel_lhv = [0.1295, 0.1198, 0.11692, 0.07472,  0.084250]'*1055; % MJ per gal
propane_GGE_yr = HEFA_propane_out_gal_yr*(fuel_lhv(5)/fuel_lhv(3)); 
LNG_GGE_yr = HEFA_LNG_out_gal_yr*(fuel_lhv(4)/fuel_lhv(3)); 
diesel_GGE_yr = HEFA_diesel_out_gal_yr*(fuel_lhv(1)/fuel_lhv(3)); 
jet_GGE_yr = HEFA_jet_out_gal_yr*(fuel_lhv(2)/fuel_lhv(3)); 
naphtha_GGE_yr = HEFA_naphtha_out_gal_yr; 
total_GGE_yr = propane_GGE_yr + LNG_GGE_yr + diesel_GGE_yr + jet_GGE_yr + naphtha_GGE_yr; 

% HEFA CAPEX ISBL (Pearlson et al., 2013; Gary et al., 2007; Chen & Quinn 2021)
% HEFA Hydrotreater capital costs
HEFA_HT_scaling_value = 6524; %bpd
HEFA_HT_price = 27000000; %2007 dollars 
HEFA_HT_new_scaling_value = HEFA_oil_BPD; 
HEFA_HT_install_factor = 1.51; 
HEFA_HT_scaling_exp = 0.75; 
HEFA_HT_size_ratio = HEFA_HT_new_scaling_value/HEFA_HT_scaling_value;
HEFA_HT_scaled_cost = HEFA_HT_price*HEFA_HT_size_ratio^HEFA_HT_scaling_exp; 
HEFA_HT_scaled_cost_proj_yr = HEFA_HT_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(18,2)); 
HEFA_HT_instal_cost_proj_yr = HEFA_HT_scaled_cost_proj_yr*HEFA_HT_install_factor;

% HEFA Isomerizer capital costs
ISO_feed_rate_MBSPD = HEFA_oil_BPD/1000; % thousand barrels per stream day
ISO_purchase_price = 10.404*(ISO_feed_rate_MBSPD^0.5065)*10^6; % purchase cost in 2005 USD - see excel HEFA for equation derivision
ISO_install_factor = 1.0; % assumption based on investment price from Gary et al., 2007
ISO_scaled_cost_proj_yr = ISO_purchase_price*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(16,2));
ISO_instal_cost_proj_yr = ISO_scaled_cost_proj_yr*ISO_install_factor; % intall cost in project year

% On-site SMR unit (hydrogen island) CAPEX
switch HEFA_H2_scenario
    case 'Purchased'
    SMR_instal_cost_proj_yr = 0; % no capex for SMR is hydrogen is purchased offsite
    case 'Onsite SMR'
    SMR_H2_out = HEFA_hydrogen_in; % kg H2/hr
    SMR_H2_out_MMSCFD = SMR_H2_out*24*423/1000000; % H2 out of SMR in MMSCF/day gas produced
    SMR_purchase_price = 5.2293*(SMR_H2_out_MMSCFD^0.6474)*10^6; % purchase cost in 2005 USD - see excel HEFA for equation derivision
    SMR_install_factor = 1.0; 
    SMR_scaled_cost_proj_yr = SMR_purchase_price*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(16,2));
    SMR_instal_cost_proj_yr = SMR_scaled_cost_proj_yr*SMR_install_factor; % intall cost in project year
end

% Gas processing plant (to recover valuable smaller hydrocarbons) CAPEX
HEFA_C3_C4_out = (HEFA_propane_out_gal_day + HEFA_LNG_gal_day); % gallons propane and LPG per day
GPP_gas_feed_MMSCF_day = HEFA_C3_C4_out/16.2*1000/1000000; % MMSCF of gas feed/day; Gary et al., 2007 - table 13.6 relates total liq products from GPP to gas feed
GPP_C3_C4_per_MSCF = HEFA_C3_C4_out/(GPP_gas_feed_MMSCF_day*1000); %#ok<NASGU> % Relates amount of propane and heavier to gas feed (closest to 15 line in Figure 13.4 in Gary)
GPP_purchase_price = 5.1531*(GPP_gas_feed_MMSCF_day^0.6038)*10^6; % purchase cost in 2005 USD - see excel for derivision
GPP_install_factor = 1.0; 
GPP_scaled_cost_proj_yr = GPP_purchase_price*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(16,2));
GPP_instal_cost_proj_yr = GPP_scaled_cost_proj_yr*GPP_install_factor; 

% % HEFA CAPEX OSBL
% %Cooling tower system (Sized with Gary et al., 2007; costed with Chen & Quinn 2021)
HEFA_CW_HT = 500*(HEFA_oil_BPD/24)*3.785*2.205; % lbs water circulated per hour for hydrotreating (Gary)
HEFA_CW_ISO = 800*(HEFA_oil_BPD/24)*3.785*2.205; % lbs water per hour for isomerization (Gary)
HEFA_CW_SMR = 65*(HEFA_hydrogen_in*2.205)*3.785*2.205; % lbs water per hour for SMR (Gary)
HEFA_CW_GPP = 100*(HEFA_C3_C4_out)*3.785*2.205; % lbs water per hour for GPP (Gary)
HEFA_total_cooling_water = HEFA_CW_HT + HEFA_CW_ISO + HEFA_CW_GPP + HEFA_CW_SMR; 
HEFA_CTS_scaling_value = 35631668; %lb/hr
HEFA_CTS_price = 2000000; %2009 dollars 
HEFA_CTS_new_scaling_value = HEFA_total_cooling_water; 
HEFA_CTS_install_factor = 2.95; 
HEFA_CTS_scaling_exp = 0.6; 
HEFA_CTS_size_ratio = HEFA_CTS_new_scaling_value/HEFA_CTS_scaling_value;
HEFA_CTS_scaled_cost = HEFA_CTS_price*HEFA_CTS_size_ratio^HEFA_CTS_scaling_exp; 
HEFA_CTS_scaled_cost_proj_yr = HEFA_CTS_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(20,2)); 
HEFA_CTS_instal_cost_proj_yr = HEFA_CTS_scaled_cost_proj_yr*HEFA_CTS_install_factor;

% %Cooling water pump (Sized with Gary et al., 2007; costed with Chen & Quinn 2021)
HEFA_CWP_scaling_value = 35631668; %lb/hr
HEFA_CWP_price = 445700; %2009 dollars 
HEFA_CWP_new_scaling_value = HEFA_total_cooling_water; 
HEFA_CWP_install_factor = 2.95; 
HEFA_CWP_scaling_exp = 0.6; 
HEFA_CWP_size_ratio = HEFA_CWP_new_scaling_value/HEFA_CWP_scaling_value;
HEFA_CWP_scaled_cost = HEFA_CWP_price*HEFA_CWP_size_ratio^HEFA_CWP_scaling_exp; 
HEFA_CWP_scaled_cost_proj_yr = HEFA_CWP_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(20,2)); 
HEFA_CWP_instal_cost_proj_yr = HEFA_CWP_scaled_cost_proj_yr*HEFA_CWP_install_factor;

% Crude Storage - 13 days (Gary et al., 2007)
total_crude_storage = HEFA_oil_BPD*13; % barrels storage
crude_storage_CAPEX = total_crude_storage*100; % $100 per barrel storage capacity ($50 - $130)

% Liquid products storage - 25 days (Gary et al., 2007)
total_liquid_storage = HEFA_naphtha_out_gal_day + HEFA_jet_out_gal_day + HEFA_diesel_out_gal_day + HEFA_LNG_gal_day + HEFA_propane_out_gal_day;
liquid_storage_CAPEX = 25*(total_liquid_storage/42)*100; % $100 per barrel storage capacity ($50 - $130) - 42 gallons per barrel

% Firewater pump (Chen & Quinn 2021)
HEFA_FWP_scaling_value = 2000; %tonnes per day
HEFA_FWP_price = 184000; %1997 dollars 
HEFA_FWP_new_scaling_value = (ann_ave_centr*TSS)*24/1000; % AFDW Flowrate to refinery
HEFA_FWP_install_factor = 2.95; 
HEFA_FWP_scaling_exp = 0.79; 
HEFA_FWP_size_ratio = HEFA_FWP_new_scaling_value/HEFA_FWP_scaling_value;
HEFA_FWP_scaled_cost = HEFA_FWP_price*HEFA_FWP_size_ratio^HEFA_FWP_scaling_exp; 
HEFA_FWP_scaled_cost_proj_yr = HEFA_FWP_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(8,2)); 
HEFA_FWP_instal_cost_proj_yr = HEFA_FWP_scaled_cost_proj_yr*HEFA_FWP_install_factor;

% Firewater storage tank (Chen & Quinn 2021)
HEFA_FWST_scaling_value = 2000; %tonnes per day
HEFA_FWST_price = 166100; %1997 dollars 
HEFA_FWST_new_scaling_value = (ann_ave_centr*TSS)*24/1000; % AFDW Flowrate to refinery
HEFA_FWST_install_factor = 2.95; 
HEFA_FWST_scaling_exp = 0.51; 
HEFA_FWST_size_ratio = HEFA_FWST_new_scaling_value/HEFA_FWST_scaling_value;
HEFA_FWST_scaled_cost = HEFA_FWST_price*HEFA_FWST_size_ratio^HEFA_FWST_scaling_exp; 
HEFA_FWST_scaled_cost_proj_yr = HEFA_FWST_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(8,2)); 
HEFA_FWST_instal_cost_proj_yr = HEFA_FWST_scaled_cost_proj_yr*HEFA_FWST_install_factor;

% HEFA variable operational costs
% HEFA process utilites (annual values)
switch HEFA_H2_scenario
    case 'Purchased'
    HEFA_BFW = (0.25)*HEFA_oil_in_lbs*op_hours/2.205; % kg/yr 
    HEFA_cooling_water = (2.55 + 5.26)*HEFA_oil_in_lbs*op_hours/2.205; % kg/yr 
    HEFA_elec_en = (0.01 + 0.01)*HEFA_oil_in_lbs*op_hours; % kWh/yr
    HEFA_H2_opex = HEFA_hydrogen_in*op_hours*H_cost; 
    HEFA_SMR_cat_opex = 0; 

        switch HEFA_scenario
            case 'Maximum Distillate'
            HEFA_NG = (0.02 + 0.03)*HEFA_oil_in_lbs*op_hours/2.205; % kg/yr
            case 'Maximum Jet'
            HEFA_NG = (0.02 + 0.03 + 0.06)*HEFA_oil_in_lbs*op_hours/2.205; % kg/yr
        end
       
    case 'Onsite SMR'
    HEFA_BFW = (0.25 + 0.5)*HEFA_oil_in_lbs*op_hours/2.205; % kg/yr 
    HEFA_cooling_water = (2.55 + 5.26 + 5.33)*HEFA_oil_in_lbs*op_hours/2.205; % kg/yr 
    HEFA_elec_en = (0.01 + 0.01 + 0.02)*HEFA_oil_in_lbs*op_hours; % kWh/yr
    HEFA_H2_opex = 0; % for onsite SMR, opex for H2 is covered by additional capex and utilities for SMR
    HEFA_SMR_cat_opex = 0.8/100*(HEFA_hydrogen_in*op_hours*2.205); 

        switch HEFA_scenario
            case 'Maximum Distillate'
            HEFA_NG = (0.02 + 0.03 + 0.06)*HEFA_oil_in_lbs*op_hours/2.205; % kg/yr
            case 'Maximum Jet'
            HEFA_NG = (0.02 + 0.03 + 0.06 + 0.06)*HEFA_oil_in_lbs*op_hours/2.205; % kg/yr
        end
end 

% HEFA opex
HEFA_elec_en_opex = HEFA_elec_en*Elec_cost; %#ok<NASGU> 
HEFA_NG_opex = HEFA_NG*NG_cost;
HEFA_ann_therm_kWh = HEFA_NG*ng_boiler_eff*(1/0.022)*983/3412; % kWh_therm per year including boiler eff.
HEFA_water_opex = (HEFA_BFW + HEFA_cooling_water)*Process_water_cost;

% HEFA catalyst CAPEX
HEFA_ISO_icat_CAPEX = 150*HEFA_oil_BPD*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(16,2)); 

% HEFA catalyst OPEX - (catalysts and chemicals for SMR is cases above)
HEFA_HT_cat_opex = 0.06*(HEFA_oil_BPD/24)*op_hours*(CHEM_INDEX(CHEM_IND_POS,2)/CHEM_INDEX(26,2));  
HEFA_ISO_cat_opex = 0.08*(HEFA_oil_BPD/24)*op_hours*(CHEM_INDEX(CHEM_IND_POS,2)/CHEM_INDEX(26,2));  
HEFA_SMR_cat_opex = HEFA_SMR_cat_opex*1;

%----------------------------------------------------------------------------------------------------------------------------------------------------

% Anaerobic digestion of LEA (ANL 2011)
AD_inflow_water = ann_ave_centr*(1-TSS); %#ok<NASGU> 
AD_inflow_carbs = ann_ave_centr*TSS*carb; 
AD_inflow_proteins = ann_ave_centr*TSS*prot; 
AD_inflow_ash = ann_ave_centr*TSS*ash; 
AD_inflow_lip = ann_ave_centr*TSS*lip - HEFA_oil_in; 
AD_inflow_TS = AD_inflow_carbs + AD_inflow_proteins + AD_inflow_ash + AD_inflow_lip; % kg TS/hr
AD_inflow_VS = 0.9*AD_inflow_TS; % kg VS/hr

% AD methane yields
AD_CH4_yield = 0.3*(AD_inflow_TS*1000); % L CH4 per hour; ANL 2011 - 0.3 L CH4/g-TS entering AD
AD_biogas_yield = AD_CH4_yield/0.67; % L of biogas per hour; 67% CH4 in biogas; ANL 2011
AD_ann_CH4_yield = AD_CH4_yield*op_hours; %#ok<NASGU> 
AD_ann_biogas_yield = AD_biogas_yield*op_hours; %#ok<NASGU> 

% AD outputs (nutrient recovery - C, N, and P mass balances)
% Carbon balance
C_content_LEA = 0.44; % ANL 2011 - carbon into AD
AD_C_in = (AD_inflow_TS - AD_inflow_ash)*C_content_LEA; % kg C/hr
AD_C_to_methane = AD_CH4_yield/1000*0.657*(12/16.04); % kg C per hour from methane
AD_C_to_CO2_biogas = (AD_biogas_yield*(1-0.67))/1000*1.87*(12/44); % kg C per hour from CO2
AD_C_out = AD_C_in - (AD_C_to_methane + AD_C_to_CO2_biogas); % kg C per hour
AD_C_to_digestate = 0.5*AD_C_out; % kg C per hour 
AD_C_to_supernatant = 0.5*AD_C_out; % kg C per hour
AD_seq_credit = 0.08*AD_C_to_digestate*(44/12); % kg CO2-eq/hr
AD_ann_seq_credit = AD_seq_credit*op_hours/1000; % tCO2-eq/yr

% Nitrogen balance
AD_N_in = (AD_inflow_TS - AD_inflow_ash)*n_content; % kg N/hr
AD_N_out_NH3eq_aqueous = AD_N_in*0.76*(17.031/14.0067); %kg NH3eq/hr - 80% N retained as NH3, 5% volatized - 76% total retention of N in aqueous
AD_N_out_digestate = (AD_N_in*0.20); % kg N/hr
AD_N_out_BA = (AD_N_in*0.20)*0.40; % 20% of N to organic N in digestate with 40% bioavailability (BA)
AD_ann_NH3_credit = AD_N_out_NH3eq_aqueous*op_hours*(1-nut_rem); 
WWT_N_removal_credit = AD_ann_NH3_credit*nut_rem*0.82*15; % Units are $/yr 82% N in NH3; $15/kg N
AD_ann_oN_credit = AD_N_out_BA*op_hours; %#ok<NASGU> 

% N2O emissions from soil application of AD digestate
AD_N2O_soil_em = 0.01*AD_N_out_digestate; % kg N2O-N; 0.01 kg N2O-N/kg N applied
AD_N2O_ann_soil_em = AD_N2O_soil_em*op_hours; % kg N2O-N per year

% Phosphours balance
AD_P_in = (AD_inflow_TS - AD_inflow_ash)*p_content; % kg P/hr
AD_P_out_DAPeq_SN = 0.50*AD_P_in*(132.06/30.973762); % kg DAP-eq/hr
AD_P_out_oP_BA = 0.5*AD_P_in; % kg organic N per hr
AD_ann_DAP_credit = AD_P_out_DAPeq_SN*op_hours*(1-nut_rem); % kg DAP/yr 
WWT_P_removal_credit = AD_ann_DAP_credit*nut_rem*0.20*100; % units are $/yr; 20% N in NH3; $100/kg P
AD_ann_oP_credit = AD_P_out_oP_BA*op_hours; %#ok<NASGU> % kg organic P/yr

% AD utilities
AD_ann_heat_input = 0.68*(AD_inflow_TS)*op_hours; % kWh/yr - 0.68 kWh_thermal per kg TS
AD_ann_elec_input = 0.122*(AD_inflow_VS)*op_hours; % kWh/yr - 0.122 kWh_elec per kg VS
AD_ann_cent_elec_input = 0.014*(AD_inflow_VS)*op_hours; % kWh/yr - 0.014 kWh_elec per kg VS
AD_ann_elec_kWh_yr = AD_ann_elec_input+AD_ann_cent_elec_input;

% AD CAPEX
total_cult_area_ha = wetted_acres/2.471; 
areal_prod = str2double(cult_out(2,2)); 
AD_capex = 242.57*(areal_prod^1.1153)*total_cult_area_ha*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(21,2)); %Zamalloa et al., 2011 - Cost curve in excel based on areal prod.

% AD OPEX
AD_heating_eff = 0.80; % assumption
AD_NG_input = ((AD_ann_heat_input/AD_heating_eff)*3412/983)*0.022; % kg NG - 3412 BTU/kWh; LHV = 983 BTU/SCF; 0.022 kg/SCF
AD_elec_opex = AD_ann_elec_kWh_yr*Elec_cost; %#ok<NASGU> 
AD_NG_opex = AD_NG_input*NG_cost;

%----------------------------------------------------------------------------------------------------------------------------------------------------

% Biogas cleanup with pressure swing absorption
PSA_CH4_in_kg = AD_CH4_yield/1000*0.657; % kg CH4/hr
PSA_biogas_in_kg = AD_biogas_yield/1000*1.158; %#ok<NASGU> % kg biogas/hr
PSA_biogas_in_scf = (AD_biogas_yield/1000)*35.31; % scf biogas/hr
PSA_fug_CH4 = 0.03*PSA_CH4_in_kg; % fugitive methane loss from AD and cleaning - Borjesson and Berglund 2006
PSA_ann_fug_CH4 = PSA_fug_CH4*op_hours; % kg CH4 leaked per year - Borjesson and Berglund 2006

% PSA outputs
PSA_CH4_out_kg_hr = PSA_CH4_in_kg - PSA_fug_CH4; % kg CH4 coming from PSA including fugitive losses
PSA_CO2_out_kg_hr = AD_biogas_yield*(1-0.67)/1000*1.87; % kg CO2 per hour out of PSA; CO2 in = CO2 out

% PSA CAPEX - scaled on total volume of biogas cleaned per day
PSA_scaling_value = 10; %MMscf/d
PSA_price = 1750000; %2004 dollars 
PSA_new_scaling_value = PSA_biogas_in_scf*24/10^6; % MMSCF/d
PSA_install_factor = 2.47; 
PSA_scaling_exp = 0.8; 
PSA_size_ratio = PSA_new_scaling_value/PSA_scaling_value;
PSA_scaled_cost = PSA_price*PSA_size_ratio^PSA_scaling_exp; 
PSA_scaled_cost_proj_yr = PSA_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(15,2)); 
PSA_instal_cost_proj_yr = PSA_scaled_cost_proj_yr*PSA_install_factor;

% PSA - energy OPEX
HEFA_PSA_ann_elec_en = 0.25*(AD_biogas_yield/1000*1.158)*op_hours; % biogas density of 1.158 kg/Nm3 - 1000L/m3; 0.25 kWh/Nm3 biogas for PSA
HEFA_PSA_ann_elec_en_opex = HEFA_PSA_ann_elec_en*Elec_cost; %#ok<NASGU> 

%----------------------------------------------------------------------------------------------------------------------------------------------------

% Combined Heat and Power
CHP_CH4_in_kg = PSA_CH4_out_kg_hr;
CHP_CH4_leaked = 0.02*CHP_CH4_in_kg; 
CHP_CH4_utilized = CHP_CH4_in_kg - CHP_CH4_leaked; 
CHP_energy_in = CHP_CH4_utilized*13.89; % kWh per hour - ANL 2011 - 13.89 kWh/kg CH4

% CHP outputs (electric and thermal energy)
CHP_ann_elec_out = 0.33*CHP_energy_in*op_hours; % 33% electrical efficiency kWh_elec out per year
CHP_ann_therm_out = 0.76*(CHP_energy_in*op_hours) - CHP_ann_elec_out; % kWh_therm out per year

% CHP emissions matrix [VOC; CO; NOx; PM10; PM2.5; SOx; CH4; N2O; CO2] all
% in g/mmBTU of fuel input
CHP_em_mat = [1; 24; 113; 3.1; 3.1; 0.269; 4.26; 1.5; 59360]; 
CHP_fuel_input_mmBTU = CHP_energy_in*op_hours*3412/1000000; % MMBTU/yr of fuel input
CHP_emissions = CHP_em_mat.*CHP_fuel_input_mmBTU./1000; %kg of each 

% CHP CAPEX
CHP_kW_rating = CHP_ann_elec_out/op_hours; 
CHP_MW_rating = CHP_kW_rating/1000; 
CHP_installed_CAPEX = 1384.2*(CHP_MW_rating^-0.249)*CHP_kW_rating; % installed CAPEX in USD/kW x kW rating

%----------------------------------------------------------------------------------------------------------------------------------------------------

% Total CAPEX
% LE CAPEX: HPH and WHE total CAPEX - for minimum lipid selling price - option to include AD and CHP for energy and nutrient recovery from LEA

    switch LE_inc_AD_CHP
        case 'Yes'
        LE_CAPEX_MAT = [HPH_scaled_install_cost; WHE_install_cost; WHE_initial_hexane_charge; AD_capex; PSA_instal_cost_proj_yr; CHP_installed_CAPEX];
        case 'No'
        LE_CAPEX_MAT = [HPH_scaled_install_cost; WHE_install_cost; WHE_initial_hexane_charge; 0; 0; 0];
    end 

LE_CAPEX = sum(LE_CAPEX_MAT);

% Entire SAF Process CAPEX: HPH/WHE/HEFA/AD/CHP total CAPEX
SAF_CAPEX_MAT = [HPH_scaled_install_cost; WHE_install_cost; WHE_initial_hexane_charge; HEFA_HT_instal_cost_proj_yr; ISO_instal_cost_proj_yr; HEFA_ISO_icat_CAPEX; SMR_instal_cost_proj_yr; GPP_instal_cost_proj_yr; HEFA_CTS_instal_cost_proj_yr; HEFA_CWP_instal_cost_proj_yr; crude_storage_CAPEX; liquid_storage_CAPEX; HEFA_FWP_instal_cost_proj_yr; HEFA_FWST_instal_cost_proj_yr; AD_capex; PSA_instal_cost_proj_yr; CHP_installed_CAPEX];
SAF_CAPEX = sum(SAF_CAPEX_MAT); 

SAF_ISBL_LE = LE_CAPEX; 
SAF_ISBL_total = SAF_CAPEX; 

% CAPEX for Lipid Extraction system boundary
    %Direct costs
    LE_warehouse = 0.04*SAF_ISBL_LE; 
    LE_site_devel = 0.09*SAF_ISBL_LE; 
    LE_add_piping = 0.045*SAF_ISBL_LE; 
    LE_TDC = LE_warehouse + LE_site_devel + LE_add_piping + SAF_ISBL_LE; 
    LE_TDC_MAT = [SAF_ISBL_LE; LE_warehouse; LE_site_devel; LE_add_piping; LE_TDC];
   
    %Indirect costs
    %LE_prorate_exp = 0.10*LE_TDC; 
    %LE_field_exp = 0.10*LE_TDC;
    %LE_home_off_const = 0.20*LE_TDC;
    LE_contingency = 0.10*LE_TDC;
    %LE_other_costs = 0.10*LE_TDC;
    
    %LE_TINC = LE_prorate_exp + LE_field_exp + LE_home_off_const + LE_contingency + LE_other_costs;
    %LE_TINC_MAT = [LE_prorate_exp; LE_field_exp; LE_home_off_const; LE_contingency; LE_other_costs; LE_TINC]; 
    LE_TINC = LE_contingency;

    %Fixed capital investment
    LE_FCI = LE_TDC + LE_TINC; 
    LE_working_cap = 0.05*LE_FCI; 
    
    %Total capital investment (TCI)
    LE_TCI = LE_FCI + LE_working_cap; 
    
    % Grouped outputs
    LE_CAPEX_SUM(1:5,1) = LE_TDC_MAT; 
    LE_CAPEX_SUM(6,1) = LE_contingency; 
    LE_CAPEX_SUM(7,1) = LE_FCI;
    LE_CAPEX_SUM(8,1) = LE_working_cap; 
    LE_CAPEX_SUM(9,1) = LE_TCI; 

% CAPEX for complete SAF system boundary
    SAF_warehouse = 0.04*SAF_ISBL_total; 
    SAF_site_devel = 0.09*SAF_ISBL_total; 
    SAF_add_piping = 0.045*SAF_ISBL_total; 
    SAF_TDC = SAF_warehouse + SAF_site_devel + SAF_add_piping + SAF_ISBL_total; 
    SAF_TDC_MAT = [SAF_ISBL_total; SAF_warehouse; SAF_site_devel; SAF_add_piping; SAF_TDC]; 

    %Indirect costs
    %SAF_prorate_exp = 0.10*SAF_TDC; 
    %SAF_field_exp = 0.10*SAF_TDC;
    %SAF_home_off_const = 0.20*SAF_TDC;
    SAF_contingency = 0.10*SAF_TDC;
    %SAF_other_costs = 0.10*SAF_TDC;
    
    %SAF_TINC = SAF_prorate_exp + SAF_field_exp + SAF_home_off_const + SAF_contingency + SAF_other_costs;
    %SAF_TINC_MAT = [SAF_prorate_exp; SAF_field_exp; SAF_home_off_const; SAF_contingency; SAF_other_costs; SAF_TINC];
    SAF_TINC = SAF_contingency; 

    %Fixed capital investment
    SAF_FCI = SAF_TDC + SAF_TINC; 
    SAF_working_cap = 0.05*SAF_FCI; 
    
    %Total capital investment (TCI)
    SAF_TCI = SAF_FCI + SAF_working_cap; 

    % Grouped outputs
    SAF_CAPEX_SUM(1:5,1) = SAF_TDC_MAT; 
    SAF_CAPEX_SUM(6,1) = SAF_contingency; 
    SAF_CAPEX_SUM(7,1) = SAF_FCI;
    SAF_CAPEX_SUM(8,1) = SAF_working_cap; 
    SAF_CAPEX_SUM(9,1) = SAF_TCI; 

% Formulate LE CAPEX matrix
LE_CAPEX_out(1,1) = "High-Pressure Homogenizer"; 
LE_CAPEX_out(2,1) = "Wet Hexane Extraction Equipment";
LE_CAPEX_out(3,1) = "Initial Hexane Charge"; 
LE_CAPEX_out(4,1) = "Anaerobic Digester";
LE_CAPEX_out(5,1) = "Pressure Swing Absorber (PSA)";
LE_CAPEX_out(6,1) = "Combined Heat and Power Unit";
LE_CAPEX_out(1:6,2) = LE_CAPEX_MAT; 
LE_CAPEX_out(7,1) = "Total ISBL"; 
LE_CAPEX_out(8,1) = "Warehouse"; 
LE_CAPEX_out(9,1) = "Additional Site Development";
LE_CAPEX_out(10,1) = "Additional Piping";
LE_CAPEX_out(11,1) = "Total Direct Costs";
LE_CAPEX_out(12,1) = "Project Contingency";
LE_CAPEX_out(13,1) = "Fixed Capital Investment";
LE_CAPEX_out(14,1) = "Working Capital"; 
LE_CAPEX_out(15,1) = "Total Capital Investment"; 
LE_CAPEX_out(7:15,2) = LE_CAPEX_SUM; 

% Formulate SAF CAPEX matrix
SAF_CAPEX_out(1,1) = "High-Pressure Homogenizer"; 
SAF_CAPEX_out(2,1) = "Wet Hexane Extraction Equipment";
SAF_CAPEX_out(3,1) = "Initial Hexane Charge";
SAF_CAPEX_out(4,1) = "Hydrotreater";
SAF_CAPEX_out(5,1) = "Isomerizer"; 
SAF_CAPEX_out(6,1) = "Isomerizer Initial Catalyst Charge"; 
SAF_CAPEX_out(7,1) = "Steam Methane Reformer (SMR)";
SAF_CAPEX_out(8,1) = "Gas Processing Plant";
SAF_CAPEX_out(9,1) = "Cooling Tower System";
SAF_CAPEX_out(10,1) = "Cooling Water Pump"; 
SAF_CAPEX_out(11,1) = "Biocrude Storage (13 days)";
SAF_CAPEX_out(12,1) = "Liquid Products Storage (25 days)";
SAF_CAPEX_out(13,1) = "Fire Water Pump";
SAF_CAPEX_out(14,1) = "Fire Water Storage Tank";
SAF_CAPEX_out(15,1) = "Anaerobic Digester";
SAF_CAPEX_out(16,1) = "Pressure Swing Absorber (PSA)";
SAF_CAPEX_out(17,1) = "Combined Heat and Power Unit";
SAF_CAPEX_out(1:17,2) = SAF_CAPEX_MAT; 
SAF_CAPEX_out(18,1) = "Total ISBL"; 
SAF_CAPEX_out(19,1) = "Warehouse"; 
SAF_CAPEX_out(20,1) = "Additional Site Development";
SAF_CAPEX_out(21,1) = "Additional Piping";
SAF_CAPEX_out(22,1) = "Total Direct Costs";
SAF_CAPEX_out(23,1) = "Project Contingency";
SAF_CAPEX_out(24,1) = "Fixed Capital Investment";
SAF_CAPEX_out(25,1) = "Working Capital"; 
SAF_CAPEX_out(26,1) = "Total Capital Investment"; 
SAF_CAPEX_out(18:26,2) = SAF_CAPEX_SUM;

%----------------------------------------------------------------------------------------------------------------------------------------------------

% Energy Balance and Total OPEX
% Plant Energy Balance (electricity and NG OPEX)
% Lipid Extraction System Boundary
    switch LE_inc_AD_CHP
        case 'Yes'
        LE_elec_energy_matrix = [HPH_ann_energy_kWh; WHE_ann_elec_kWh; AD_ann_elec_kWh_yr; HEFA_PSA_ann_elec_en; -1*CHP_ann_elec_out]; 
        LE_therm_energy_matrix = [WHE_ann_thermal_kWh; AD_ann_heat_input; -1*CHP_ann_therm_out];
        case 'No'
        LE_elec_energy_matrix = [HPH_ann_energy_kWh; WHE_ann_elec_kWh]; 
        LE_therm_energy_matrix = WHE_ann_thermal_kWh;
    end 
LE_net_electric = sum(LE_elec_energy_matrix); 
LE_net_thermal = sum(LE_therm_energy_matrix); 

% Electricity and NG OPEX (LE)
    % Determine net elec. consumption and determine OPEX
    if LE_net_electric < 0 
        LE_ex_elec = -1*LE_net_electric; %excess elec. becomes a revenue stream
        LE_elec_opex = 0; 
    else 
        LE_ex_elec = 0; 
        LE_elec_opex = LE_net_electric*Elec_cost; % USD/yr
    end 
    
    % Determine net therm. energy consumption and determine OPEX
    if LE_net_thermal < 0 
        LE_ex_therm = -1*LE_net_thermal; %#ok<NASGU> %excess therm. is a zero-value waste product
        LE_NG_opex = 0; 
    else 
        LE_ex_therm = 0; %#ok<NASGU> 
        LE_NG_opex = (LE_net_thermal/ng_boiler_eff*3412/983*0.022)*NG_cost; % USD/yr
    end 

% Full SAF Pathway System Boundary
SAF_elec_energy_matrix = [HPH_ann_energy_kWh; WHE_ann_elec_kWh; HEFA_elec_en; AD_ann_elec_kWh_yr; HEFA_PSA_ann_elec_en; -1*CHP_ann_elec_out]; 
SAF_therm_energy_matrix = [WHE_ann_thermal_kWh; HEFA_ann_therm_kWh; AD_ann_heat_input; -1*CHP_ann_therm_out];
SAF_net_electric = sum(SAF_elec_energy_matrix); % negative means more energy produced than consumed; 
SAF_net_thermal = sum(SAF_therm_energy_matrix); 

% Electricity and NG OPEX (SAF)
    % Determine net elec. consumption and determine OPEX
    if SAF_net_electric < 0 
        SAF_ex_elec = -1*SAF_net_electric; %excess elec. becomes a revenue stream
        SAF_elec_opex = 0; 
    else 
        SAF_ex_elec = 0; 
        SAF_elec_opex = SAF_net_electric*Elec_cost; % USD/yr
    end 
    
    % Determine net therm. energy consumption and determine OPEX
    if SAF_net_thermal < 0 
        SAF_ex_therm = -1*SAF_net_thermal; %#ok<NASGU> %excess therm. is a zero-value waste product
        SAF_NG_opex = 0; 
    else 
        SAF_ex_therm = 0; %#ok<NASGU> 
        SAF_NG_opex = (SAF_net_thermal/ng_boiler_eff*3412/983*0.022)*NG_cost; % USD/yr - apply ng_boiler_eff because supplemental heat comes from NG bioler
    end 

SAF_NG_gross_opex = WHE_NG_opex + HEFA_NG_opex + AD_NG_opex; %#ok<NASGU> % Gross NG opex to look at savings from heat int. with AD

% CO2 Recycle from HEFA and AD/PSA
total_HEFA_CO2_out = (HEFA_CO2_out + PSA_CO2_out_kg_hr)*op_hours/1000; % tCO2 per year
total_LE_CO2_out = PSA_CO2_out_kg_hr*op_hours/1000; 
total_cult_CO2 = str2double(cult_out(10,2)); % tCO2/yr
net_HEFA_cult_CO2 = total_cult_CO2 - total_HEFA_CO2_out; 
net_LE_cult_CO2 = total_cult_CO2 - total_LE_CO2_out; 

if net_HEFA_cult_CO2 < 0 
    HEFA_CO2_credit = total_cult_CO2*CO2_cost; 
else
    HEFA_CO2_credit = total_HEFA_CO2_out*CO2_cost; 
end

if net_LE_cult_CO2 < 0 
    LE_CO2_credit = total_cult_CO2*CO2_cost; 
else
    LE_CO2_credit = total_LE_CO2_out*CO2_cost; 
end


% Consumables Matrix (hexane, H2, catalysts, water) 
SAF_consum = [WHE_hexane_makeup; HEFA_HT_cat_opex; HEFA_ISO_cat_opex; HEFA_water_opex; HEFA_SMR_cat_opex; HEFA_H2_opex; SAF_elec_opex; SAF_NG_opex]; 
LE_consum = [WHE_hexane_makeup; LE_elec_opex; LE_NG_opex]; 

% Total variable OPEX
SAF_total_var_opex = sum(SAF_consum); 
LE_total_var_opex = sum(LE_consum); 

SAF_consum_mat = [WHE_hexane_makeup; HEFA_HT_cat_opex; HEFA_ISO_cat_opex; HEFA_water_opex; HEFA_SMR_cat_opex; HEFA_H2_opex; SAF_elec_opex; SAF_NG_opex; SAF_total_var_opex]; 
SAF_OPEX_out(1,1) = "WHE Annual Hexane OPEX";
SAF_OPEX_out(2,1) = "HEFA HT Catalyst OPEX";
SAF_OPEX_out(3,1) = "HEFA ISO Catalyst OPEX";
SAF_OPEX_out(4,1) = "HEFA Water OPEX"; 
SAF_OPEX_out(5,1) = "HEFA SMR Catalyst OPEX"; 
SAF_OPEX_out(6,1) = "HEFA Purchased H2 OPEX"; 
SAF_OPEX_out(7,1) = "Supplemental Electricity OPEX";
SAF_OPEX_out(8,1) = "Supplemental Natural Gas OPEX"; 
SAF_OPEX_out(9,1) = "Total Variable OPEX"; 
SAF_OPEX_out(1:9,2) = SAF_consum_mat; 

LE_consum_mat = [WHE_hexane_makeup; LE_elec_opex; LE_NG_opex; LE_total_var_opex]; 
LE_OPEX_out(1,1) = "WHE Annual Hexane OPEX";
LE_OPEX_out(2,1) = "Supplemental Electricity OPEX";
LE_OPEX_out(3,1) = "Supplemental Natural Gas OPEX"; 
LE_OPEX_out(4,1) = "Total Variable OPEX"; 
LE_OPEX_out(1:4,2) = LE_consum_mat; 

% Plant fixed OPEX (labor, maintenance, insurance, taxes, etc.)
%Labor - all salaries quoted in 2009 USD
switch fixed_opex_case
    case 'Chen and Quinn'
        plant_mngr = 1 * 147000 * (LABOR_INDEX(LABOR_IND_POS,2)/LABOR_INDEX(20,2));
        plant_engr = 2 * 70000 * (LABOR_INDEX(LABOR_IND_POS,2)/LABOR_INDEX(20,2));
        maint_sup = 1 * 57000 * (LABOR_INDEX(LABOR_IND_POS,2)/LABOR_INDEX(20,2));
        maint_tech = 5 * 40000 * (LABOR_INDEX(LABOR_IND_POS,2)/LABOR_INDEX(20,2));
        lab_mngr = 1 * 56000 * (LABOR_INDEX(LABOR_IND_POS,2)/LABOR_INDEX(20,2));
        lab_tech = 2 * 40000 * (LABOR_INDEX(LABOR_IND_POS,2)/LABOR_INDEX(20,2));
        shft_sup = 2 * 48000 * (LABOR_INDEX(LABOR_IND_POS,2)/LABOR_INDEX(20,2));
        shft_op = 9 * 40000 * (LABOR_INDEX(LABOR_IND_POS,2)/LABOR_INDEX(20,2));
        yrd_emp = 2 * 28000 * (LABOR_INDEX(LABOR_IND_POS,2)/LABOR_INDEX(20,2));
        clk_sec = 2 * 36000 * (LABOR_INDEX(LABOR_IND_POS,2)/LABOR_INDEX(20,2));

        SAF_labor_total = (wetted_acres/5000)*(plant_mngr + plant_engr + maint_sup + maint_tech + lab_mngr + lab_tech + shft_sup + shft_op + yrd_emp + clk_sec);
        SAF_labor_burden = SAF_labor_total * 0.90; 
        
        LE_labor_total = SAF_labor_total/2; 
        LE_labor_burden = SAF_labor_burden/2;

        % Maintenance and Inssurance

        SAF_maintenance = 0.03 * SAF_ISBL_total; 
        SAF_prop_ins_tax = 0.01 * SAF_FCI;
        SAF_misc_supplies = 0; 
        SAF_fixed_opex_cont = 0;
        SAF_Total_fixed_OPEX = SAF_labor_total + SAF_labor_burden + SAF_maintenance + SAF_prop_ins_tax + SAF_misc_supplies + SAF_fixed_opex_cont;
        SAF_fixed_opex_mat = [SAF_labor_total; SAF_labor_burden; SAF_maintenance; SAF_prop_ins_tax; SAF_misc_supplies; SAF_fixed_opex_cont; SAF_Total_fixed_OPEX]; 

        LE_maintenance = 0.03 * SAF_ISBL_LE; 
        LE_prop_ins_tax = 0.01 * LE_FCI;
        LE_misc_supplies = 0;
        LE_fixed_opex_cont = 0;
        LE_Total_fixed_OPEX = LE_labor_total + LE_labor_burden + LE_maintenance + LE_prop_ins_tax + LE_misc_supplies + LE_fixed_opex_cont;
        LE_fixed_opex_mat = [LE_labor_total; LE_labor_burden; LE_maintenance; LE_prop_ins_tax; LE_misc_supplies; LE_fixed_opex_cont; LE_Total_fixed_OPEX]; 

    case 'Pearlson' 
        SAF_labor_total = 12*72000*(LABOR_INDEX(LABOR_IND_POS,2)/LABOR_INDEX(24,2)); % 12 employees at 72k (2013 USD)
        SAF_labor_burden = 0;
        
        LE_labor_total = SAF_labor_total/2; 
        LE_labor_burden = 0; 

        SAF_maintenance = 0.055 * SAF_FCI; % 5.5% for maintenance
        SAF_prop_ins_tax = 0.015 * SAF_FCI; % 0.5% for insu. and 1% for taxes
        SAF_misc_supplies = 0.002 * SAF_FCI; 
        SAF_fixed_opex_cont = 0.10*(SAF_labor_total+SAF_maintenance+SAF_prop_ins_tax+ SAF_misc_supplies);
        SAF_Total_fixed_OPEX = SAF_labor_total + SAF_labor_burden + SAF_maintenance + SAF_prop_ins_tax + SAF_misc_supplies + SAF_fixed_opex_cont;
        SAF_fixed_opex_mat = [SAF_labor_total; SAF_labor_burden; SAF_maintenance; SAF_prop_ins_tax; SAF_misc_supplies; SAF_fixed_opex_cont; SAF_Total_fixed_OPEX]; 

        LE_maintenance = 0.055 * LE_FCI; 
        LE_prop_ins_tax = 0.015 * LE_FCI;
        LE_misc_supplies = 0.002 * LE_FCI; 
        LE_fixed_opex_cont = 0.10*(LE_labor_total + LE_maintenance + LE_prop_ins_tax + LE_misc_supplies);
        LE_Total_fixed_OPEX = LE_labor_total + LE_labor_burden + LE_maintenance + LE_prop_ins_tax + LE_misc_supplies + LE_fixed_opex_cont;
        LE_fixed_opex_mat = [LE_labor_total; LE_labor_burden; LE_maintenance; LE_prop_ins_tax; LE_misc_supplies; LE_fixed_opex_cont; LE_Total_fixed_OPEX]; 

end 

% add fixed opex to output matrix:
SAF_OPEX_out(10,1) = "Labor OPEX";
SAF_OPEX_out(11,1) = "Labor Burden OPEX";
SAF_OPEX_out(12,1) = "Equipment Maintenance";
SAF_OPEX_out(13,1) = "Property Tax and Insurance"; 
SAF_OPEX_out(14,1) = "Misc. Supplies"; 
SAF_OPEX_out(15,1) = "Fixed OPEX Contingency"; 
SAF_OPEX_out(16,1) = "Total Fixed OPEX"; 
SAF_OPEX_out(10:16,2) = SAF_fixed_opex_mat; 
        
LE_OPEX_out(5,1) = "Labor OPEX";
LE_OPEX_out(6,1) = "Labor Burden OPEX";
LE_OPEX_out(7,1) = "Equipment Maintenance";
LE_OPEX_out(8,1) = "Property Tax and Insurance"; 
LE_OPEX_out(9,1) = "Misc. Supplies"; 
LE_OPEX_out(10,1) = "Fixed OPEX Contingency"; 
LE_OPEX_out(11,1) = "Total Fixed OPEX"; 
LE_OPEX_out(5:11,2) = LE_fixed_opex_mat;

%----------------------------------------------------------------------------------------------------------------------------------------------------

% Revnue Streams and Carbon Sequestration from AD
% CO-Products:NH3-eq; DAP-eq; CO2 sequestration; Excess energy (elec/therm)
SAF_coprod_NH3eq_rev = AD_ann_NH3_credit*ammonia_cost + WWT_N_removal_credit; % kg NH3_eq/yr * $/kg NH3 + $/yr
SAF_coprod_DAPeq_rev = AD_ann_DAP_credit*dap_cost + WWT_P_removal_credit; % kg DAP_eq/yr * $/kg DAP + $/yr
SAF_coprod_carb_offset_rev = AD_ann_seq_credit*carbon_offset_cost; %#ok<NASGU> % tCO2-eq seq/yr * $/tCO2 seq
SAF_excess_en_rev = SAF_ex_elec*Elec_cost; % kWh/yr * $/kWh - excess thermal energy is a burden free waste product

    switch LE_inc_AD_CHP
        case 'Yes'
        LE_excess_en_rev = LE_ex_elec*Elec_cost; 
        LE_coprod_NH3eq_rev = SAF_coprod_NH3eq_rev; 
        LE_coprod_DAPeq_rev = SAF_coprod_DAPeq_rev; 
        case 'No'
        LE_excess_en_rev = 0;
        LE_coprod_NH3eq_rev = 0;
        LE_coprod_DAPeq_rev = 0;
    end 

    switch HEFA_CO2_rec
        case 'Yes' 
        LE_coprod_CO2_rev = LE_CO2_credit; 
        SAF_coprod_CO2_rev = HEFA_CO2_credit; 
        case 'No'
        LE_coprod_CO2_rev = 0; 
        SAF_coprod_CO2_rev = 0;
    end 

% Total coproduct revenue for lipid extraction system boundary
LE_total_coprod_rev = LE_excess_en_rev + LE_coprod_NH3eq_rev + LE_coprod_DAPeq_rev + LE_coprod_CO2_rev; 
LE_coprod_mat = [LE_excess_en_rev; LE_coprod_NH3eq_rev; LE_coprod_DAPeq_rev; LE_coprod_CO2_rev];

% Fuel Products: Propane; LPG; Naphtha; Diesel; Jet
SAF_rev_propane = HEFA_propane_out_gal_yr*value_propane_gal;
SAF_rev_LNG = HEFA_LNG_out*op_hours*2.205/6.059*value_LNG_gal; 
SAF_rev_diesel =  HEFA_diesel_out_gal_yr*value_diesel_gal; 
SAF_rev_jet = HEFA_jet_out_gal_yr*value_jet_gal;
SAF_rev_naphtha = HEFA_naphtha_out_gal_yr*value_naphtha_gal;

switch HEFA_scenario
    case 'Maximum Distillate'
    SAF_total_coprod_rev = SAF_coprod_NH3eq_rev + SAF_coprod_DAPeq_rev + SAF_excess_en_rev + SAF_rev_propane + SAF_rev_LNG + SAF_rev_naphtha + SAF_rev_jet + SAF_coprod_CO2_rev + protein_revenue; %#ok<NASGU> 
    SAF_total_coprod_mat = [SAF_coprod_NH3eq_rev; SAF_coprod_DAPeq_rev; SAF_excess_en_rev; SAF_rev_propane; SAF_rev_LNG; SAF_rev_naphtha; SAF_rev_jet; SAF_coprod_CO2_rev];
    case 'Maximum Jet'
    SAF_total_coprod_rev = SAF_coprod_NH3eq_rev + SAF_coprod_DAPeq_rev + SAF_excess_en_rev + SAF_rev_propane + SAF_rev_LNG + SAF_rev_naphtha + SAF_rev_diesel + SAF_coprod_CO2_rev + protein_revenue; %#ok<NASGU> 
    SAF_total_coprod_mat = [SAF_coprod_NH3eq_rev; SAF_coprod_DAPeq_rev; SAF_excess_en_rev; SAF_rev_propane; SAF_rev_LNG; SAF_rev_naphtha; SAF_rev_diesel; SAF_coprod_CO2_rev];
end

%----------------------------------------------------------------------------------------------------------------------------------------------------

% Process Model Outputs Matrix
SAF_liquid_fuels = [HEFA_diesel_out_gal_yr; HEFA_jet_out_gal_yr; HEFA_naphtha_out_gal_yr; HEFA_LNG_out_gal_yr; HEFA_propane_out_gal_yr];
SAF_proc_mod_mat_1  = [WHE_ann_hexane_loss; HEFA_BFW + HEFA_cooling_water; HEFA_hydrogen_in*op_hours; HPH_ann_energy_kWh; WHE_ann_elec_kWh; HEFA_elec_en; AD_ann_elec_kWh_yr; HEFA_PSA_ann_elec_en; CHP_ann_elec_out; SAF_elec_opex/Elec_cost; WHE_ann_thermal_kWh; HEFA_ann_therm_kWh; AD_ann_heat_input; CHP_ann_therm_out; SAF_NG_opex/NG_cost]; 
SAF_proc_mod_mat_2 = [AD_ann_NH3_credit; AD_ann_DAP_credit; AD_N2O_ann_soil_em; HEFA_CO2_out*op_hours/1000; PSA_CO2_out_kg_hr*op_hours/1000; AD_ann_seq_credit; CHP_kW_rating; CHP_MW_rating; CHP_CH4_leaked*op_hours]; 

switch LE_inc_AD_CHP
    case 'Yes'
    LE_proc_mod_mat_1  = [WHE_ann_hexane_loss; HPH_ann_energy_kWh; WHE_ann_elec_kWh; AD_ann_elec_kWh_yr; HEFA_PSA_ann_elec_en; CHP_ann_elec_out; LE_elec_opex/Elec_cost; WHE_ann_thermal_kWh; AD_ann_heat_input; CHP_ann_therm_out; LE_NG_opex/NG_cost];
    LE_proc_mod_mat_2 = [AD_ann_NH3_credit; AD_ann_DAP_credit; AD_N2O_ann_soil_em; PSA_CO2_out_kg_hr*op_hours/1000; AD_ann_seq_credit; CHP_kW_rating; CHP_MW_rating; CHP_CH4_leaked*op_hours; HEFA_oil_in*op_hours/1000]; 
    case 'No'
    LE_proc_mod_mat_1  = [WHE_ann_hexane_loss; HPH_ann_energy_kWh; WHE_ann_elec_kWh; 0; 0; 0; LE_elec_opex/Elec_cost; WHE_ann_thermal_kWh; 0; 0; LE_NG_opex/NG_cost];
    LE_proc_mod_mat_2 = [0; 0; 0; 0; 0; 0; 0; 0; HEFA_oil_in*op_hours/1000];
end

% SAF Process Model Outputs Matrix
% Consumables
SAF_proc_mod_out(1,1) = "WHE Annual Hexane kg/yr";
SAF_proc_mod_out(2,1) = "HEFA Water Consumption kg/yr"; 
SAF_proc_mod_out(3,1) = "HEFA H2 Input"; 
SAF_proc_mod_out(4,1) = "HPH Electric Energy kWh/yr";
SAF_proc_mod_out(5,1) = "WHE Electric Energy kWh/yr";
SAF_proc_mod_out(6,1) = "HEFA Electric Energy kWh/yr";
SAF_proc_mod_out(7,1) = "AD Electric Energy kWh/yr";
SAF_proc_mod_out(8,1) = "Biogas Cleanup Electric Energy kWh/yr";
SAF_proc_mod_out(9,1) = "CHP Electricity Out kWh/yr";
SAF_proc_mod_out(10,1) = "Supplemental Electricity kWh/yr";
SAF_proc_mod_out(11,1) = "WHE Thermal Energy kWh/yr";
SAF_proc_mod_out(12,1) = "HEFA Thermal Energy kWh/yr";
SAF_proc_mod_out(13,1) = "AD Thermal Energy kWh/yr"; 
SAF_proc_mod_out(14,1) = "CHP Thermal Energy Out kWh/yr";
SAF_proc_mod_out(15,1) = "Supplemental Natural Gas kg/yr"; 
% Fuel output
SAF_proc_mod_out(16,1) = "Diesel gal/yr";
SAF_proc_mod_out(17,1) = "Jet gal/yr"; 
SAF_proc_mod_out(18,1) = "Naphtha gal/yr";
SAF_proc_mod_out(19,1) = "LNG gal/yr";
SAF_proc_mod_out(20,1) = "Propane gal/yr";
% co products
SAF_proc_mod_out(21,1) = "NH3-eq Recovery kg/yr";
SAF_proc_mod_out(22,1) = "DAP-eq Recovery kg/yr";
SAF_proc_mod_out(23,1) = "N2O Emissions to Soil kg N2O-N/yr";
% Carbon
SAF_proc_mod_out(24,1) = "HEFA CO2 Out tCO2/yr";
SAF_proc_mod_out(25,1) = "PSA CO2 Out tCO2/yr"; 
SAF_proc_mod_out(26,1) = "Digestate CO2 Sequestration tCO2/yr";
% CHP ratings
SAF_proc_mod_out(27,1) = "CHP Power Rating kW";
SAF_proc_mod_out(28,1) = "CHP Power Rating MW"; 
SAF_proc_mod_out(29,1) = "CHP Leaked kg CH4/yr";

SAF_proc_mod_out(1:15,2) = SAF_proc_mod_mat_1; 
SAF_proc_mod_out(16:20,2) = SAF_liquid_fuels;
SAF_proc_mod_out(21:29,2) = SAF_proc_mod_mat_2;

% LE Process Model Outputs Matrix
% Consumables
LE_proc_mod_out(1,1) = "WHE Annual Hexane kg/yr";
LE_proc_mod_out(2,1) = "HPH Electric Energy kWh/yr";
LE_proc_mod_out(3,1) = "WHE Electric Energy kWh/yr";
LE_proc_mod_out(4,1) = "AD Electric Energy kWh/yr";
LE_proc_mod_out(5,1) = "Biogas Cleanup Electric Energy kWh/yr";
LE_proc_mod_out(6,1) = "CHP Electricity Out kWh/yr";
LE_proc_mod_out(7,1) = "Supplemental Electricity kWh/yr";
LE_proc_mod_out(8,1) = "WHE Thermal Energy kWh/yr";
LE_proc_mod_out(9,1) = "AD Thermal Energy kWh/yr"; 
LE_proc_mod_out(10,1) = "CHP Thermal Energy Out kWh/yr";
LE_proc_mod_out(11,1) = "Supplemental Natural Gas kg/yr"; 
% co products
LE_proc_mod_out(12,1) = "NH3-eq Recovery kg/yr";
LE_proc_mod_out(13,1) = "DAP-eq Recovery kg/yr";
LE_proc_mod_out(14,1) = "N2O Emissions to Soil kg N2O-N/yr";
% Carbon
LE_proc_mod_out(15,1) = "PSA CO2 Out tCO2/yr"; 
LE_proc_mod_out(16,1) = "Digestate CO2 Sequestration tCO2/yr";
% CHP ratings
LE_proc_mod_out(17,1) = "CHP Power Rating kW";
LE_proc_mod_out(18,1) = "CHP Power Rating MW"; 
LE_proc_mod_out(19,1) = "CHP Leaked kg CH4/yr";
% Lipid Output
LE_proc_mod_out(20,1) = "Annual Lipid Production (tonnes algae-oil/yr)"; 
    
LE_proc_mod_out(1:11,2) = LE_proc_mod_mat_1; 
LE_proc_mod_out(12:20,2) = LE_proc_mod_mat_2;

%----------------------------------------------------------------------------------------------------------------------------------------------------

% Discounted Cash Flow Rate of Return
% Minimum Fuel Selling Prices - Conversion and Cultivation
switch HEFA_scenario
    case 'Maximum Distillate'
    SAF_total_coprod_rev = SAF_coprod_NH3eq_rev + SAF_coprod_DAPeq_rev + SAF_excess_en_rev; %+ SAF_rev_propane + SAF_rev_LNG + SAF_rev_naphtha + SAF_rev_jet; 
    
    % MDSP (diesel)
    % DCFROR for Minimum Diesel Selling Price
    Gal_per_yr = total_GGE_yr; %GGE/yr
    land_CAPEX = str2double(cult_TEA_out_CAPEX(59,2)); %from cultivation (includes area for downstream processing)
    total_equip_fac_CAPEX = str2double(cult_TEA_out_CAPEX(58,2)) + SAF_FCI;
    total_working_CAPEX = 0.05*total_equip_fac_CAPEX;
    total_combined_OPEX = str2double(cult_TEA_out_OPEX(9,2)) + SAF_total_var_opex + SAF_Total_fixed_OPEX;
    total_coprod_rev = SAF_total_coprod_rev;

    options = optimset('Display','off');
    
    fun = @NPV_MFSP;
    x0 = 0.1;
    MDSP = fsolve(fun, x0, options, Gal_per_yr, land_CAPEX, total_equip_fac_CAPEX, total_working_CAPEX, total_combined_OPEX, IRR, total_coprod_rev);
    [NPV_mfsp, cash_flow_mfsp] = NPV_MFSP(MDSP, Gal_per_yr, land_CAPEX, total_equip_fac_CAPEX, total_working_CAPEX, total_combined_OPEX, IRR, total_coprod_rev); %#ok<ASGLU> 
  
    DCFROR_out(1,1) = "Minimum Fuel Selling Price w/ Cultivation"; 
    DCFROR_out(1,2) = MDSP; 
    DCFROR_out(2,1) = "Minimum Jet Fuel Selling Price w/ Cultivation"; 
    DCFROR_out(2,2) = "NA"; 

    case 'Maximum Jet'
    SAF_total_coprod_rev = SAF_coprod_NH3eq_rev + SAF_coprod_DAPeq_rev + SAF_excess_en_rev + protein_revenue; % + SAF_rev_propane + SAF_rev_LNG + SAF_rev_naphtha + SAF_rev_diesel   --> uncomment for MFSP in $/gal SAF
    
    % MFSP (per GGE basis)
    % DCFROR for Minimum Fuel Selling Price
    Gal_per_yr = total_GGE_yr; %GGE/yr
    land_CAPEX = str2double(cult_TEA_out_CAPEX(59,2)); %from cultivation (includes area for downstream processing)
    total_equip_fac_CAPEX = str2double(cult_TEA_out_CAPEX(58,2)) + SAF_FCI;
    total_working_CAPEX = 0.05*total_equip_fac_CAPEX;
    total_combined_OPEX = str2double(cult_TEA_out_OPEX(9,2)) + SAF_total_var_opex + SAF_Total_fixed_OPEX;
    total_coprod_rev = SAF_total_coprod_rev;
    
    options = optimset('Display','off');
    
    fun = @NPV_MFSP;
    x0 = 0.1;
    MJSP = fsolve(fun, x0, options, Gal_per_yr, land_CAPEX, total_equip_fac_CAPEX, total_working_CAPEX, total_combined_OPEX, IRR, total_coprod_rev);
    [NPV_mfsp, cash_flow_mfsp] = NPV_MFSP(MJSP, Gal_per_yr, land_CAPEX, total_equip_fac_CAPEX, total_working_CAPEX, total_combined_OPEX, IRR, total_coprod_rev); %#ok<ASGLU> 
    
    DCFROR_out(1,1) = "Minimum Diesel Selling Price w/ Cultivation"; 
    DCFROR_out(1,2) = "NA"; 
    DCFROR_out(2,1) = "Minimum Fuel Selling Price w/ Cultivation"; 
    DCFROR_out(2,2) = MJSP; 
end 
               
% MLSP (lipids)
% DCFROR for Minimum Lipid Selling Price
kg_lip_per_yr = WHE_oil_out*op_hours; %kg lip/yr
land_CAPEX = str2double(cult_TEA_out_CAPEX(59,2)); %from cultivation (includes area for downstream processing)
total_equip_fac_CAPEX = str2double(cult_TEA_out_CAPEX(58,2)) + LE_FCI;
total_working_CAPEX = 0.05*total_equip_fac_CAPEX;
total_combined_OPEX = str2double(cult_TEA_out_OPEX(9,2)) + LE_total_var_opex + LE_Total_fixed_OPEX; 
total_coprod_rev = LE_total_coprod_rev;

options = optimset('Display','off');

fun = @NPV_MFSP;
x0 = 0.1;
MLSP = fsolve(fun, x0, options, kg_lip_per_yr, land_CAPEX, total_equip_fac_CAPEX, total_working_CAPEX, total_combined_OPEX, IRR, total_coprod_rev);
[NPV_mfsp, cash_flow_mfsp] = NPV_MFSP(MLSP, kg_lip_per_yr, land_CAPEX, total_equip_fac_CAPEX, total_working_CAPEX, total_combined_OPEX, IRR, total_coprod_rev); %#ok<ASGLU> 

DCFROR_out(3,1) = "Minimum Lipid Selling Price"; 
DCFROR_out(3,2) = MLSP;

    % Minimum Fuel Selling Prices - Conversion Only
    switch HEFA_scenario
        case 'Maximum Distillate'
        SAF_total_coprod_rev = SAF_coprod_NH3eq_rev + SAF_coprod_DAPeq_rev + SAF_excess_en_rev; %+ SAF_rev_propane + SAF_rev_LNG + SAF_rev_naphtha + SAF_rev_jet; 
        
        % MDSP (diesel)
        % DCFROR for Minimum Diesel Selling Price
        Gal_per_yr = total_GGE_yr; %GGE/yr
        land_CAPEX = 0; %land assumed to be available near/on cultivation site
        total_equip_fac_CAPEX = SAF_FCI;
        total_working_CAPEX = 0.05*total_equip_fac_CAPEX;
        total_combined_OPEX = SAF_total_var_opex + SAF_Total_fixed_OPEX + HEFA_oil_in*op_hours*lip_purch_price; 
        total_coprod_rev = SAF_total_coprod_rev;
        
        options = optimset('Display','off');
        
        fun = @NPV_MFSP;
        x0 = 0.1;
        MDSP_CO = fsolve(fun, x0, options, Gal_per_yr, land_CAPEX, total_equip_fac_CAPEX, total_working_CAPEX, total_combined_OPEX, IRR, total_coprod_rev);
        [NPV_mfsp, cash_flow_mfsp] = NPV_MFSP(MDSP_CO, Gal_per_yr, land_CAPEX, total_equip_fac_CAPEX, total_working_CAPEX, total_combined_OPEX, IRR, total_coprod_rev); %#ok<ASGLU> 
        
        DCFROR_out(4,1) = "Minimum Fuel Selling Price (Conversion Only)"; 
        DCFROR_out(4,2) = MDSP_CO;
        DCFROR_out(5,1) = "Minimum Jet Fuel Selling Price (Conversion Only)"; 
        DCFROR_out(5,2) = "NA";
    
        case 'Maximum Jet'
        SAF_total_coprod_rev = SAF_coprod_NH3eq_rev + SAF_coprod_DAPeq_rev + SAF_excess_en_rev; % + SAF_rev_propane + SAF_rev_LNG + SAF_rev_naphtha + SAF_rev_diesel; 
        
        % MFSP (per GGE)
        % DCFROR for Minimum Fuel Selling Price (Conversion Only)
        Gal_per_yr = total_GGE_yr; %Gal/yr
        land_CAPEX = 0; %land assumed to be available near/on cultivation site
        total_equip_fac_CAPEX = SAF_FCI;
        total_working_CAPEX = 0.05*total_equip_fac_CAPEX;
        total_combined_OPEX = SAF_total_var_opex + SAF_Total_fixed_OPEX + HEFA_oil_in*op_hours*lip_purch_price;
        total_coprod_rev = SAF_total_coprod_rev;
        
        options = optimset('Display','off');
        
        fun = @NPV_MFSP;
        x0 = 0.1;
        MJSP_CO = fsolve(fun, x0, options, Gal_per_yr, land_CAPEX, total_equip_fac_CAPEX, total_working_CAPEX, total_combined_OPEX, IRR, total_coprod_rev);
        [NPV_mfsp, cash_flow_mfsp] = NPV_MFSP(MJSP_CO, Gal_per_yr, land_CAPEX, total_equip_fac_CAPEX, total_working_CAPEX, total_combined_OPEX, IRR, total_coprod_rev); %#ok<ASGLU> 
      
        DCFROR_out(4,1) = "Minimum Diesel Selling Price (Conversion Only)"; 
        DCFROR_out(4,2) = "NA";
        DCFROR_out(5,1) = "Minimum Fuel Selling Price (Conversion Only)"; 
        DCFROR_out(5,2) = MJSP_CO;
    end

%----------------------------------------------------------------------------------------------------------------------------------------------------

% LCA Consumables Matrix

% Electric Energy 
SAF_LCA_consum(1,1) = "HPH Electricity kWh/yr"; 
SAF_LCA_consum(2,1) = "WHE Electricity kWh/yr"; 
SAF_LCA_consum(3,1) = "HEFA Electricity kWh/yr"; 
SAF_LCA_consum(4,1) = "AD Electricity kWh/yr"; 
SAF_LCA_consum(5,1) = "PSA Electricity kWh/yr"; 
SAF_LCA_consum(6,1) = "CHP Electricity Out (produced = negative) kWh/yr"; 
SAF_LCA_consum(7,1) = "Net Supplemental Electricity (from GRID) kWh/yr"; 
SAF_LCA_consum(1:6,2) = SAF_elec_energy_matrix; %[HPH_ann_energy_kWh; WHE_ann_elec_kWh; HEFA_elec_en; AD_ann_elec_kWh_yr; HEFA_PSA_ann_elec_en; -1*CHP_ann_elec_out]; 
SAF_LCA_consum(7,2) = sum(SAF_elec_energy_matrix);

% Thermal Energy and supplemental NG - if excess thermal energy from CHP it
% is a zero-value and zero-burden waste stream
SAF_LCA_consum(8,1) = "WHE Thermal Energy kWh/yr";
SAF_LCA_consum(9,1) = "HEFA Thermal Energy kWh/yr";
SAF_LCA_consum(10,1) = "AD Thermal Energy kWh/yr"; 
SAF_LCA_consum(11,1) = "CHP Thermal Energy Out (produced = negative) kWh/yr";
SAF_LCA_consum(8:11,2) = SAF_therm_energy_matrix;
SAF_LCA_consum(12,1) = "Net Supplemental Thermal Energy (from NG) kWh_thermal/yr";
SAF_LCA_consum(12,2) = sum(SAF_therm_energy_matrix); 
SAF_LCA_consum(13,1) = "Supplemental Natural Gas kg/yr"; 
SAF_LCA_consum(13,2) = SAF_NG_opex/NG_cost;

% Consumables (hexane, water consumption, hydrogen)
SAF_LCA_consum(14,1) = "WHE Initial Hexane (one-time emission in year 1) kg"; 
SAF_LCA_consum(15,1) = "WHE Annual Hexane Makeup kg/yr";
SAF_LCA_consum(16,1) = "HEFA Water Consumption kg/yr"; 
SAF_LCA_consum(17,1) = "HEFA H2 Input (kg/yr) - if onsite SMR selected then H2 input is zero"; 
SAF_materials = [WHE_hexane_input; WHE_ann_hexane_loss; HEFA_BFW + HEFA_cooling_water; HEFA_hydrogen_in*op_hours]; 
SAF_LCA_consum(14:17,2) = SAF_materials; 

% Direct process CO2 emissions - option to recycle to ponds
SAF_LCA_consum(18,1) = "HEFA CO2 Out tCO2/yr";
SAF_LCA_consum(18,2) = HEFA_CO2_out*op_hours/1000; 
SAF_LCA_consum(19,1) = "PSA CO2 Out tCO2/yr"; 
SAF_LCA_consum(19,2) = PSA_CO2_out_kg_hr*op_hours/1000;

% Combustion Emissions
% CHP emissions matrix [VOC; CO; NOx; PM10; PM2.5; SOx; CH4; N2O; CO2]
% kg/yr of each 
SAF_LCA_consum(20,1) = "CHP - kg VOC/yr";
SAF_LCA_consum(21,1) = "CHP - kg CO/yr";
SAF_LCA_consum(22,1) = "CHP - kg NOx/yr";
SAF_LCA_consum(23,1) = "CHP - kg PM10/yr";
SAF_LCA_consum(24,1) = "CHP - kg PM2.5/yr";
SAF_LCA_consum(25,1) = "CHP - kg SOx/yr";
SAF_LCA_consum(26,1) = "CHP - kg CH4/yr";
SAF_LCA_consum(27,1) = "CHP - kg N2O/yr";
SAF_LCA_consum(28,1) = "CHP - kg CO2/yr";
SAF_LCA_consum(20:28,2) = CHP_emissions; 
SAF_LCA_consum(29,1) = "CHP - leaked (fugitive) methane - kg CH4/yr";
SAF_LCA_consum(29,2) = CHP_CH4_leaked*op_hours; 

% Fuel Products (combusted)
SAF_LCA_consum(30,1) = "Renewable Diesel gal/yr";
SAF_LCA_consum(31,1) = "Jet (kerosene) gal/yr"; 
SAF_LCA_consum(32,1) = "Naphtha gal/yr";
SAF_LCA_consum(33,1) = "LNG gal/yr";
SAF_LCA_consum(34,1) = "Propane gal/yr";
SAF_LCA_consum(30:34,2) = SAF_liquid_fuels;

% Nutrient Recovery from AD
SAF_LCA_consum(35,1) = "NH3-eq Recovery kg/yr";
SAF_LCA_consum(35,2) = AD_ann_NH3_credit; 
SAF_LCA_consum(36,1) = "DAP-eq Recovery kg/yr";
SAF_LCA_consum(36,2) = AD_ann_DAP_credit;

% Other relevant flows
SAF_LCA_consum(37,1) = "N2O Emissions to Soil kg N2O-N/yr";
SAF_LCA_consum(37,2) = AD_N2O_ann_soil_em;
SAF_LCA_consum(38,1) = "Digestate CO2 Sequestration tCO2/yr";
SAF_LCA_consum(38,2) = AD_ann_seq_credit;
SAF_LCA_consum(39,1) = "Waste Water to Wastewater Treatment (m3/yr)"; 
SAF_LCA_consum(39,2) = HEFA_water_out*op_hours/1000;    

SAF_LCA_consum(40,1) = "Monthly water consumption (kg per month)";
SAF_LCA_consum(40,2) = (HEFA_BFW + HEFA_cooling_water)/12; %kg/yr spread over 12 months

SAF_LCA_consum(41,1) = "Renewable Diesel (gal/month)";
SAF_LCA_consum(42,1) = "Jet (kerosene) (gal/month)"; 
SAF_LCA_consum(43,1) = "Naphtha (gal/month)";
SAF_LCA_consum(44,1) = "LNG (gal/month)";
SAF_LCA_consum(45,1) = "Propane (gal/month)";
SAF_LCA_consum(41:45,2) = SAF_liquid_fuels/12; %gallons per year spread over 12 months

SAF_LCA_consum(46,1) = "C to AD Digestate kg C/yr";
SAF_LCA_consum(46,2) = AD_C_to_digestate*op_hours; 
SAF_LCA_consum(47,1) = "C to AD Supernatant kg C/yr";
SAF_LCA_consum(47,2) = AD_C_to_supernatant*op_hours; 

SAF_LCA_consum(48,1) = "PSA Fugitive CH4 kg/yr";
SAF_LCA_consum(48,2) = PSA_ann_fug_CH4;



end