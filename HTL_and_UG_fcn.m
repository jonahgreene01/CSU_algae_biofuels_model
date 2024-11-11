function[HTL_process_mod_out, HTL_UG_CAPX_OPX, direct_cap, case_number, biogenic_c_in_gas_phase] = HTL_and_UG_fcn(prot, carb, lip, ash, areal_prod, fac_size, res_time, reactor_diam, annual_uptime, Char_recovery, NH3EQ_recovery, DAPEQ_recovery, Cost_year, NG_cost, H_cost, Elec_cost, Process_water_cost, NH3EQ_value, DAPEQ_value, Char_value, cult_out, nut_rem) 
%HTL and upgrading module

load background_data; %#ok<*LOAD>
comp_yield = table2array(comp_and_yield);
comp_diff = zeros(231,4); 
cult_outs = str2double(cult_out(:,2)); 

% Calculate AFDW composition
prot_AFDW = prot/(prot + carb + lip); 
carb_AFDW = carb/(prot + carb + lip);
lip_AFDW = lip/(prot + carb + lip);

%% Determine case number and aspen model outputs based on composition of algae stream
for i = 1:231

    comp_diff(i,1) = abs(comp_yield(i, 2) - prot_AFDW); 
    comp_diff(i,2) = abs(comp_yield(i, 3) - carb_AFDW);
    comp_diff(i,3) = abs(comp_yield(i, 4) - lip_AFDW);
    comp_diff(i,4) =  comp_diff(i,1) + comp_diff(i,2) + comp_diff(i,3);
    
    if i == 1 
        min_diff = comp_diff(i,4); 
        case_number = 1; 
    else
        if comp_diff(i,4) < min_diff
            min_diff = comp_diff(i,4); 
            case_number = i; 
        end
    end
    
end

%Determine HTL phase yields based on AFDW composition and ash - Based on
%additive phase model from Leow et al., 2018 
biocrude_yld = comp_yield(case_number, 5)*(1-ash); 
aqueous_yld = comp_yield(case_number, 6)*(1-ash); 
solids_yld = comp_yield(case_number, 7)*(1-ash); 
gas_yld = comp_yield(case_number, 8)*(1-ash); 

%Determine ASPEN model outputs based on case number and interpolate based
%on ash content
    if ash > 0 && ash <= 0.10
        zero_ash = table2array(zero_percent_ash);
        zero_ash_BCN = zero_ash(case_number,:);
        ten_ash = table2array(ten_percent_ash);
        ten_ash_BCN = ten_ash(case_number,:);

        aspen_outputs = zero_ash_BCN + (ash - 0).*((ten_ash_BCN - zero_ash_BCN)./(0.1-0)); 

    elseif ash > 0.10 && ash <= 0.20
        ten_ash = table2array(ten_percent_ash);
        ten_ash_BCN = ten_ash(case_number,:);
        twenty_ash = table2array(twenty_percent_ash);
        twenty_ash_BCN = twenty_ash(case_number,:);

        aspen_outputs = ten_ash_BCN + (ash - 0.1).*((twenty_ash_BCN - ten_ash_BCN)./(0.2-0.1)); 

    elseif ash > 0.20 && ash <= 0.30
        twenty_ash = table2array(twenty_percent_ash);
        twenty_ash_BCN = twenty_ash(case_number,:);
        thirty_ash = table2array(thirty_percent_ash);
        thirty_ash_BCN = thirty_ash(case_number,:);

        aspen_outputs = twenty_ash_BCN + (ash - 0.2).*((thirty_ash_BCN - twenty_ash_BCN)./(0.3-0.2)); 

    elseif ash > 0.30 && ash <= 0.40
        thirty_ash = table2array(thirty_percent_ash);
        thirty_ash_BCN = thirty_ash(case_number,:);
        fourty_ash = table2array(fourty_percent_ash);
        fourty_ash_BCN = fourty_ash(case_number,:);

        aspen_outputs = thirty_ash_BCN + (ash - 0.3).*((fourty_ash_BCN - thirty_ash_BCN)./(0.4-0.3)); 

    end

%Scale aspen outputs by facility size and productivity (some variables do
%not scale with size or productivity)
aspen_outputs_scaled(1,1) = aspen_outputs(1,1);
aspen_outputs_scaled(1, 2:7) = aspen_outputs(1,2:7).*(areal_prod/25).*(fac_size/5000);
aspen_outputs_scaled(1,8:9) = aspen_outputs(1,8:9);
aspen_outputs_scaled(1, 10:20) = aspen_outputs(1,10:20).*(areal_prod/25).*(fac_size/5000);
aspen_outputs_scaled(1,21) = aspen_outputs(1,21);
aspen_outputs_scaled(1, 22) = aspen_outputs(1,22).*(areal_prod/25).*(fac_size/5000);
aspen_outputs_scaled(1,23:25) = aspen_outputs(1,23:25);
aspen_outputs_scaled(1, 26:28) = aspen_outputs(1,26:28).*(areal_prod/25).*(fac_size/5000);
aspen_outputs_scaled(1,29) = aspen_outputs(1,29);
aspen_outputs_scaled(1, 30) = aspen_outputs(1,30).*(areal_prod/25).*(fac_size/5000);
aspen_outputs_scaled(1,31:39) = aspen_outputs(1,31:39);
aspen_outputs_scaled(1, 40:45) = aspen_outputs(1,40:45).*(areal_prod/25).*(fac_size/5000);
aspen_outputs_scaled(1,46:47) = aspen_outputs(1,46:47);
aspen_outputs_scaled(1, 48:51) = aspen_outputs(1,48:51).*(areal_prod/25).*(fac_size/5000);

%get rid of case number so position one is first variable
asp_out_fin = aspen_outputs_scaled(1, 2:51); 

%% Calculate heat exchanger areas and reactor lengths

%Heat transfer coefficients in kW/sqft/K
U_E301 = 0.07595; 
U_E304 = 0.08122; 
U_E305 = 0.05010; 
U_E306 = 0.05274; 

%Cold stream temps in degrees C
CS_in_E301 = asp_out_fin(1, 30);
CS_in_E304 = asp_out_fin(1, 31);
CS_in_E305 = 90.0; %Knorr et al., 2013
CS_in_E306 = asp_out_fin(1, 46);

CS_out_E301 = asp_out_fin(1, 31);
CS_out_E304 = asp_out_fin(1, 32);
CS_out_E305 = 100.0; %Knorr et al., 2013
CS_out_E306 = asp_out_fin(1, 45);

%Hot stream temps in dgrees C
HS_in_E301 = asp_out_fin(1, 35);
HS_in_E304 = 380; %hot oil assumption
HS_in_E305 = asp_out_fin(1, 33);
HS_in_E306 = asp_out_fin(1, 37);

HS_out_E301 = asp_out_fin(1, 36);
HS_out_E304 = 350; %hot oil assumption
HS_out_E305 = asp_out_fin(1, 34);
HS_out_E306 = asp_out_fin(1, 38);

%Calculate log mean temp diff
DTLM_E301 = ((CS_out_E301 - HS_in_E301)-(CS_in_E301 - HS_out_E301))/log((CS_out_E301 - HS_in_E301)/(CS_in_E301-HS_out_E301));
DTLM_E304 = ((CS_out_E304 - HS_in_E304)-(CS_in_E304 - HS_out_E304))/log((CS_out_E304 - HS_in_E304)/(CS_in_E304-HS_out_E304));
DTLM_E305 = ((CS_out_E305 - HS_in_E305)-(CS_in_E305 - HS_out_E305))/log((CS_out_E305 - HS_in_E305)/(CS_in_E305-HS_out_E305));
DTLM_E306 = ((CS_out_E306 - HS_in_E306)-(CS_in_E306 - HS_out_E306))/log((CS_out_E306 - HS_in_E306)/(CS_in_E306-HS_out_E306));

%Determine heat duty in kW
HD_E301 = asp_out_fin(1, 10) - asp_out_fin(1, 9);
HD_E304 = asp_out_fin(1, 25)*1163*-1;
HD_E305 = asp_out_fin(1, 26)*1163;
HD_E306 = asp_out_fin(1, 27)*1163;

%HX area (sqft)
HX_A_sqft_E301 = -1*HD_E301/(U_E301*DTLM_E301);
HX_A_sqft_E304 = -1*HD_E304/(U_E304*DTLM_E304);
HX_A_sqft_E305 = -1*HD_E305/(U_E305*DTLM_E305);
HX_A_sqft_E306 = -1*HD_E306/(U_E306*DTLM_E306);

%Calculate reactor length
kg_hr = asp_out_fin(1, 21);
density = asp_out_fin(1, 22)*1000; %kg/m3
reactor_length_ft = 4*(kg_hr/density)*0.589*res_time/(3.14*reactor_diam);

%% Determine HTL and upgrading consumables and total fuel output
op_hours = annual_uptime*24; 

%Fuel properties
%Diesel
diesel_LHV = 128450; %BTU/gal
diesel_dens = 0.8366; %kg/L

%Gasoline(naptha)
naphtha_LHV = 116090; %BTU/gal
naphtha_dens = 0.7447; %kg/L

%Biocrude
biocrude_dens = asp_out_fin(1, 22); %kg/L

%Hourly flowrates (kg/hr) or (kW for electricity)
hydrogen = asp_out_fin(1, 41) + asp_out_fin(1, 42);

%HTL solids (either for waste disposal or char credits) 
HTL_solids = asp_out_fin(1, 43);

%Utilities
NG_supp = asp_out_fin(1, 2);
process_water = asp_out_fin(1, 48);
grid_elec = asp_out_fin(1, 49);

%Fuel products
diesel_kg_hr = asp_out_fin(1, 1);
naphtha_kg_hr = asp_out_fin(1, 44);
biocrude_kg_hr = asp_out_fin(1, 39);

%Co-Products
NH3EQ_rec = cult_outs(8,1)*NH3EQ_recovery; 
DAPEQ_rec = cult_outs(9,1)*DAPEQ_recovery; 

%Annual outputs (XX/yr)
annual_hydrogen = hydrogen * op_hours; 
annual_NG = NG_supp * op_hours; 
annual_water = process_water * op_hours; 
annual_elec = grid_elec * op_hours;

%Annual products (MMgal/yr) and co-products (kg/yr)
diesel_MMgal_yr = diesel_kg_hr/diesel_dens*op_hours/3.78541/1000000; 
naphtha_MMgal_yr = naphtha_kg_hr/naphtha_dens*op_hours/3.78541/1000000;
biocrude_MMgal_yr = biocrude_kg_hr/biocrude_dens*op_hours/3.78541/1000000;

%Convert to GGE
diesel_MMGGE_yr = diesel_MMgal_yr*(diesel_LHV/naphtha_LHV);
naphtha_MMGGE_yr = naphtha_MMgal_yr;
total_MMGGE_yr = diesel_MMGGE_yr + naphtha_MMGGE_yr;
total_MJ_yr = total_MMGGE_yr*1000000*naphtha_LHV*0.001055;

%Monthly products and water
total_MJ_monthly = total_MJ_yr/12; 
process_water_monthly = annual_water/12; 

annual_HTL_solids = HTL_solids * op_hours;
annual_NH3EQ = NH3EQ_rec; %previous code used ammonia output from Peters model --> NH3EQ * op_hours;
annual_DAPEQ = DAPEQ_rec; %previous code used DAP output from Peters model --> DAPEQ * op_hours;

%% Determine additional emissions from boiler combustion
direct_co2 = asp_out_fin(1,4)*op_hours; % kg CO2 per year
gas_in_boiler = direct_co2*.20; % kg per yr, total recycled gas combusted in the boiler 0.2 factor calculated from ratio of CO2 OUT to mass in from Aspen
htl_gas = asp_out_fin(1, 21)*(1-asp_out_fin(1, 20)*(1-ash))*gas_yld*op_hours; 

heat_input = (gas_in_boiler + htl_gas + annual_NG)*47.1/1055; % heat input to boiler in mmBTU per yr; gas_in_boiler and htl_gas modeled as CH4

% Natural gas boiler emission factors from GREET 2021 - g per mmBTU
% [CO2, CH4, N2O, VOC, CO, NOx, PM10 PM2.5, SOx]
if (heat_input/op_hours) > 100  % utility boiler 
    ef = [59367, 1.060, 0.75, 2.54, 22.21, 36.4, 3.507, 3.507, 0.26857]; % 
elseif (heat_input/op_hours) <= 100 % small industrial boiler
    ef = [59363, 1.060, 0.35, 2.54, 24.97, 41.05, 3.507, 3.507, 0.26857]; % 
end

direct_cap = annual_NG*ef/1000; % kg per yr
fossil_co2 = direct_cap(1); 
biogenic_c_in_gas_phase = (direct_co2 - fossil_co2); % biogenic CO2 in kg CO2/yr

% carbon from aqueous phase
c_aq_dist = .336 - lip*.422; %from Leow et al., carbon in aqueous phase based on linear function

%% Process model output matrix 
HTL_process_mod_out(1,1) = "Process Consumables";
HTL_process_mod_out(1,2) = "Annual Consumption";
HTL_process_mod_out(2,1) = "Hydrogen (kg/yr)"; 
HTL_process_mod_out(2,2) = annual_hydrogen;
HTL_process_mod_out(3,1) = "Natural Gas (kg/yr)";
HTL_process_mod_out(3,2) = annual_NG;
HTL_process_mod_out(4,1) = "Process Water (kg/yr)";
HTL_process_mod_out(4,2) = annual_water;
HTL_process_mod_out(5,1) = "Grid Electricity (kWh/yr)";
HTL_process_mod_out(5,2) = annual_elec;
HTL_process_mod_out(6,1) = "Process Outputs";
HTL_process_mod_out(6,2) = "Annual Production";
HTL_process_mod_out(7,1) = "Diesel (gal/yr)";
HTL_process_mod_out(7,2) = diesel_MMgal_yr*1000000;
HTL_process_mod_out(8,1) = "Naphtha (gal/yr)";
HTL_process_mod_out(8,2) = naphtha_MMgal_yr*1000000;
HTL_process_mod_out(9,1) = "Biocrude (gal/yr)";
HTL_process_mod_out(9,2) = biocrude_MMgal_yr*1000000;
HTL_process_mod_out(10,1) = "Diesel (GGE/yr)";
HTL_process_mod_out(10,2) = diesel_MMGGE_yr*1000000;
HTL_process_mod_out(11,1) = "Naphtha (GGE/yr)";
HTL_process_mod_out(11,2) = naphtha_MMGGE_yr*1000000;
HTL_process_mod_out(12,1) = "Total GGE/yr";
HTL_process_mod_out(12,2) = total_MMGGE_yr*1000000;
HTL_process_mod_out(13,1) = "Total MJ/yr";
HTL_process_mod_out(13,2) = total_MJ_yr;
HTL_process_mod_out(14,1) = "Co-Product Output";
HTL_process_mod_out(14,2) = "Annual Production";
HTL_process_mod_out(15,1) = "HTL Solids Recovered as Char kg/yr";
HTL_process_mod_out(15,2) = annual_HTL_solids*Char_recovery;
HTL_process_mod_out(16,1) = "NH3 Equivalent";
HTL_process_mod_out(16,2) = annual_NH3EQ; 
HTL_process_mod_out(17,1) = "DAP Equivalent";
HTL_process_mod_out(17,2) = annual_DAPEQ; 
HTL_process_mod_out(18,1) = "Biocrude Yield";
HTL_process_mod_out(18,2) = biocrude_yld;
HTL_process_mod_out(19,1) = "Aqueous Phase Yield";
HTL_process_mod_out(19,2) = aqueous_yld;
HTL_process_mod_out(20,1) = "Solids Yield";
HTL_process_mod_out(20,2) = solids_yld;
HTL_process_mod_out(21,1) = "Gas Yield";
HTL_process_mod_out(21,2) = gas_yld;
HTL_process_mod_out(22,1) = "Direct CO2 emissions (kg per year)";
HTL_process_mod_out(22,2) = direct_co2; %kg CO2 vented per year
HTL_process_mod_out(23,1) = "Monthly Energy Output (MJ per month)"; 
HTL_process_mod_out(23,2) = total_MJ_monthly; 
HTL_process_mod_out(24,1) = "Monthly Water Consumption (kg per month)"; 
HTL_process_mod_out(24,2) = process_water_monthly; 
HTL_process_mod_out(25,1) = "Carbon in Aqueous Phase (percent of C in)"; 
HTL_process_mod_out(25,2) = c_aq_dist;

%% HTL Capital expenses
CAPEX_INDEX = table2array(CAPEX_CEPI_INDEX);
CHEM_INDEX = table2array(CHEMICAL_INDEX); 
ELECT_INDEX = table2array(ELEC_INDEX);
LABOR_INDEX = table2array(LABOR_INDEX); %#ok<NODEF>

CAP_IND_POS = Cost_year - 1989; 
CHEM_IND_POS = Cost_year - 1979; 
ELEC_IND_POS = Cost_year - 2013; 
LABOR_IND_POS = Cost_year - 1989; 

%Purge/Reactor exchanger (Knorr et al., 2013)
PR_scaling_value = 49756; %sqft
PR_price = 44604000; %2012 dollars 
PR_new_scaling_value = HX_A_sqft_E301;
PR_install_factor = 2.2; 
PR_scaling_exp = 0.7; 
PR_size_ratio = PR_new_scaling_value/PR_scaling_value;
PR_scaled_cost = PR_price*PR_size_ratio^PR_scaling_exp; 
PR_scaled_cost_proj_yr = PR_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(23,2)); 
PR_instal_cost_proj_yr = PR_scaled_cost_proj_yr*PR_install_factor;

%HTL reactor
HTL_scaling_value = 480; %ft
HTL_price = 272788; %2013 dollars 
HTL_new_scaling_value = reactor_length_ft;
HTL_install_factor = 2; 
HTL_scaling_exp = 1; 
HTL_size_ratio = HTL_new_scaling_value/HTL_scaling_value;
HTL_scaled_cost = HTL_price*HTL_size_ratio^HTL_scaling_exp; 
HTL_scaled_cost_proj_yr = HTL_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(24,2)); 
HTL_instal_cost_proj_yr = HTL_scaled_cost_proj_yr*HTL_install_factor;

%Reactor gas KO Drum
KO_scaling_value = 1; %cubic foot
KO_price = 5600000; %2012 dollars 
KO_new_scaling_value = 1;
KO_install_factor = 2; 
KO_scaling_exp = 0.7; 
KO_size_ratio = KO_new_scaling_value/KO_scaling_value;
KO_scaled_cost = KO_price*KO_size_ratio^KO_scaling_exp; 
KO_scaled_cost_proj_yr = KO_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(23,2)); 
KO_instal_cost_proj_yr = KO_scaled_cost_proj_yr*KO_install_factor;

%Solids Filter
SF_scaling_value = 3689; %gpm
SF_price = 1311000; %2011 dollars 
SF_new_scaling_value = (0.264/60)*asp_out_fin(1,29)/asp_out_fin(1,28);
SF_install_factor = 1.7; 
SF_scaling_exp = 0.6; 
SF_size_ratio = SF_new_scaling_value/SF_scaling_value;
SF_scaled_cost = SF_price*SF_size_ratio^SF_scaling_exp; 
SF_scaled_cost_proj_yr = SF_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(22,2)); 
SF_instal_cost_proj_yr = SF_scaled_cost_proj_yr*SF_install_factor;

%Separator
SEP_scaling_value = 3689; %gpm
SEP_price = 3561000; %2011 dollars 
SEP_new_scaling_value = (0.264/60)*asp_out_fin(1,29)/asp_out_fin(1,28);
SEP_install_factor = 2; 
SEP_scaling_exp = 0.7; 
SEP_size_ratio = SEP_new_scaling_value/SEP_scaling_value;
SEP_scaled_cost = SEP_price*SEP_size_ratio^SEP_scaling_exp; 
SEP_scaled_cost_proj_yr = SEP_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(22,2)); 
SEP_instal_cost_proj_yr = SEP_scaled_cost_proj_yr*SEP_install_factor;

%Purge water cooler
PWC_scaling_value = 13020; %sqft
PWC_price = 255600; %2013 dollars 
PWC_new_scaling_value = HX_A_sqft_E306;
PWC_install_factor = 2.3; 
PWC_scaling_exp = 0.8; 
PWC_size_ratio = PWC_new_scaling_value/PWC_scaling_value;
PWC_scaled_cost = PWC_price*PWC_size_ratio^PWC_scaling_exp; 
PWC_scaled_cost_proj_yr = PWC_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(24,2)); 
PWC_instal_cost_proj_yr = PWC_scaled_cost_proj_yr*PWC_install_factor;

%HTL system total
HTL_INST_EQUIP_CAPEX = PR_instal_cost_proj_yr + HTL_instal_cost_proj_yr + KO_instal_cost_proj_yr + SF_instal_cost_proj_yr + SEP_instal_cost_proj_yr + PWC_instal_cost_proj_yr;

%% Upgrading Capital Expenses

%HT reactor
HT_scaling_value = 6524; %bpd
HT_price = 27000000; %2007 dollars 
HT_new_scaling_value = 0.00629*24*asp_out_fin(1,39)/asp_out_fin(1,22);
HT_install_factor = 1.51; 
HT_scaling_exp = 0.75; 
HT_size_ratio = HT_new_scaling_value/HT_scaling_value;
HT_scaled_cost = HT_price*HT_size_ratio^HT_scaling_exp; 
HT_scaled_cost_proj_yr = HT_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(18,2)); 
HT_instal_cost_proj_yr = HT_scaled_cost_proj_yr*HT_install_factor;

%Hydrogen compressor
HYD_COMP_scaling_value = 17.1; %MMscf/d
HYD_COMP_price = 1385600; %2011 dollars 
HYD_COMP_new_scaling_value = 0.0000000353*24*asp_out_fin(1,41)/0.00009;
HYD_COMP_install_factor = 1.1; 
HYD_COMP_scaling_exp = 0.8; 
HYD_COMP_size_ratio = HYD_COMP_new_scaling_value/HYD_COMP_scaling_value;
HYD_COMP_scaled_cost = HYD_COMP_price*HYD_COMP_size_ratio^HYD_COMP_scaling_exp; 
HYD_COMP_scaled_cost_proj_yr = HYD_COMP_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(22,2)); 
HYD_COMP_instal_cost_proj_yr = HYD_COMP_scaled_cost_proj_yr*HYD_COMP_install_factor;

%Pressure swing absorber
PSA_scaling_value = 10; %MMscf/d
PSA_price = 1750000; %2004 dollars 
PSA_new_scaling_value = 0.0000000353*24*asp_out_fin(1,42)/0.00009;
PSA_install_factor = 2.47; 
PSA_scaling_exp = 0.8; 
PSA_size_ratio = PSA_new_scaling_value/PSA_scaling_value;
PSA_scaled_cost = PSA_price*PSA_size_ratio^PSA_scaling_exp; 
PSA_scaled_cost_proj_yr = PSA_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(15,2)); 
PSA_instal_cost_proj_yr = PSA_scaled_cost_proj_yr*PSA_install_factor;

%Hydrocracker and auxiliaries 
HC_scaling_value = 2200; %bpd
HC_price = 25000000; %2007 dollars 
HC_new_scaling_value = 0.00629*24*asp_out_fin(1,40)/asp_out_fin(1,24);
HC_install_factor = 1.51; 
HC_scaling_exp = 0.75; 
HC_size_ratio = HC_new_scaling_value/HC_scaling_value;
HC_scaled_cost = HC_price*HC_size_ratio^HC_scaling_exp; 
HC_scaled_cost_proj_yr = HC_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(18,2)); 
HC_instal_cost_proj_yr = HC_scaled_cost_proj_yr*HC_install_factor;

UG_INST_EQUIP_CAPEX = HT_instal_cost_proj_yr + HYD_COMP_instal_cost_proj_yr + PSA_instal_cost_proj_yr + HC_instal_cost_proj_yr;

%% Heating Utility Capital Expenses

%Reactor heater
RH_scaling_value = 6032; %sqft
RH_price = 998850; %2012 dollars 
RH_new_scaling_value = HX_A_sqft_E304;
RH_install_factor = 2.2; 
RH_scaling_exp = 0.7; 
RH_size_ratio = RH_new_scaling_value/RH_scaling_value;
RH_scaled_cost = RH_price*RH_size_ratio^RH_scaling_exp; 
RH_scaled_cost_proj_yr = RH_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(23,2)); 
RH_instal_cost_proj_yr = RH_scaled_cost_proj_yr*RH_install_factor;

%Bio-oil heat recovery steam generator
BOHR_scaling_value = 868; %sqft
BOHR_price = 102000; %2012 dollars 
BOHR_new_scaling_value = HX_A_sqft_E305;
BOHR_install_factor = 2.2; 
BOHR_scaling_exp = 0.7; 
BOHR_size_ratio = BOHR_new_scaling_value/BOHR_scaling_value;
BOHR_scaled_cost = BOHR_price*BOHR_size_ratio^BOHR_scaling_exp; 
BOHR_scaled_cost_proj_yr = BOHR_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(23,2)); 
BOHR_instal_cost_proj_yr = BOHR_scaled_cost_proj_yr*BOHR_install_factor;

%Hot oil system package
HOS_scaling_value = 60; %MMBTU/hr
HOS_price = 1200100; %2012 dollars 
HOS_new_scaling_value = asp_out_fin(1,25)*-3.97;
HOS_install_factor = 1.8; 
HOS_scaling_exp = 0.6; 
HOS_size_ratio = HOS_new_scaling_value/HOS_scaling_value;
HOS_scaled_cost = HOS_price*HOS_size_ratio^HOS_scaling_exp; 
HOS_scaled_cost_proj_yr = HOS_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(23,2)); 
HOS_instal_cost_proj_yr = HOS_scaled_cost_proj_yr*HOS_install_factor;

%Hot oil
HO_scaling_value = 1; 
HO_price = 2101710; %2012 dollars 
HO_new_scaling_value = 1;
HO_install_factor = 1; 
HO_scaling_exp = 1; 
HO_size_ratio = HO_new_scaling_value/HO_scaling_value;
HO_scaled_cost = HO_price*HO_size_ratio^HO_scaling_exp; 
HO_scaled_cost_proj_yr = HO_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(23,2)); 
HO_instal_cost_proj_yr = HO_scaled_cost_proj_yr*HO_install_factor;

HEAT_UT_INST_EQUIP_CAPEX = RH_instal_cost_proj_yr + BOHR_instal_cost_proj_yr + HOS_instal_cost_proj_yr + HO_instal_cost_proj_yr;

%% Other Utilities Capital Expenses

%Cooling tower system
CTS_scaling_value = 35631668; %lb/hr
CTS_price = 2000000; %2009 dollars 
CTS_new_scaling_value = 2.205*asp_out_fin(1,47);
CTS_install_factor = 2.95; 
CTS_scaling_exp = 0.6; 
CTS_size_ratio = CTS_new_scaling_value/CTS_scaling_value;
CTS_scaled_cost = CTS_price*CTS_size_ratio^CTS_scaling_exp; 
CTS_scaled_cost_proj_yr = CTS_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(20,2)); 
CTS_instal_cost_proj_yr = CTS_scaled_cost_proj_yr*CTS_install_factor;

%Cooling water pump
CWP_scaling_value = 35631668; %lb/hr
CWP_price = 445700; %2009 dollars 
CWP_new_scaling_value = 2.205*asp_out_fin(1,47);
CWP_install_factor = 2.95; 
CWP_scaling_exp = 0.6; 
CWP_size_ratio = CWP_new_scaling_value/CWP_scaling_value;
CWP_scaled_cost = CWP_price*CWP_size_ratio^CWP_scaling_exp; 
CWP_scaled_cost_proj_yr = CWP_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(20,2)); 
CWP_instal_cost_proj_yr = CWP_scaled_cost_proj_yr*CWP_install_factor;

%Plant air compressor
PAC_scaling_value = 2000; %tonnes per day
PAC_price = 32376; %2002 dollars 
PAC_new_scaling_value = 24/1000*asp_out_fin(1,5);
PAC_install_factor = 2.95; 
PAC_scaling_exp = 0.34; 
PAC_size_ratio = PAC_new_scaling_value/PAC_scaling_value;
PAC_scaled_cost = PAC_price*PAC_size_ratio^PAC_scaling_exp; 
PAC_scaled_cost_proj_yr = PAC_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(13,2)); 
PAC_instal_cost_proj_yr = PAC_scaled_cost_proj_yr*PAC_install_factor;

%Hydraulic dump truck with scale
HDT_scaling_value = 2000; %tonnes per day
HDT_price = 80000; %1998 dollars 
HDT_new_scaling_value = 24/1000*asp_out_fin(1,5);
HDT_install_factor = 2.95; 
HDT_scaling_exp = 0.6; 
HDT_size_ratio = HDT_new_scaling_value/HDT_scaling_value;
HDT_scaled_cost = HDT_price*HDT_size_ratio^HDT_scaling_exp; 
HDT_scaled_cost_proj_yr = HDT_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(9,2)); 
HDT_instal_cost_proj_yr = HDT_scaled_cost_proj_yr*HDT_install_factor;

%Firewater pump
FWP_scaling_value = 2000; %tonnes per day
FWP_price = 184000; %1997 dollars 
FWP_new_scaling_value = 24/1000*asp_out_fin(1,5);
FWP_install_factor = 2.95; 
FWP_scaling_exp = 0.79; 
FWP_size_ratio = FWP_new_scaling_value/FWP_scaling_value;
FWP_scaled_cost = FWP_price*FWP_size_ratio^FWP_scaling_exp; 
FWP_scaled_cost_proj_yr = FWP_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(8,2)); 
FWP_instal_cost_proj_yr = FWP_scaled_cost_proj_yr*FWP_install_factor;

%Instrument air dryer
IAD_scaling_value = 2000; %tonnes per day
IAD_price = 8349; %2002 dollars 
IAD_new_scaling_value = 24/1000*asp_out_fin(1,5);
IAD_install_factor = 2.95; 
IAD_scaling_exp = 0.6; 
IAD_size_ratio = IAD_new_scaling_value/IAD_scaling_value;
IAD_scaled_cost = IAD_price*IAD_size_ratio^IAD_scaling_exp; 
IAD_scaled_cost_proj_yr = IAD_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(13,2)); 
IAD_instal_cost_proj_yr = IAD_scaled_cost_proj_yr*IAD_install_factor;

%Plant air reciever
PAR_scaling_value = 2000; %tonnes per day
PAR_price = 7003; %2002 dollars 
PAR_new_scaling_value = 24/1000*asp_out_fin(1,5);
PAR_install_factor = 2.95; 
PAR_scaling_exp = 0.72; 
PAR_size_ratio = PAR_new_scaling_value/PAR_scaling_value;
PAR_scaled_cost = PAR_price*PAR_size_ratio^PAR_scaling_exp; 
PAR_scaled_cost_proj_yr = PAR_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(13,2)); 
PAR_instal_cost_proj_yr = PAR_scaled_cost_proj_yr*PAR_install_factor;

%Firewater storage tank
FWST_scaling_value = 2000; %tonnes per day
FWST_price = 166100; %1997 dollars 
FWST_new_scaling_value = 24/1000*asp_out_fin(1,5);
FWST_install_factor = 2.95; 
FWST_scaling_exp = 0.51; 
FWST_size_ratio = FWST_new_scaling_value/FWST_scaling_value;
FWST_scaled_cost = FWST_price*FWST_size_ratio^FWST_scaling_exp; 
FWST_scaled_cost_proj_yr = FWST_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(8,2)); 
FWST_instal_cost_proj_yr = FWST_scaled_cost_proj_yr*FWST_install_factor;

%HTL biocrude intermediate storage (3 days)
BIST_scaling_value = 1056846; %gallons
BIST_price = 470000; %2005 dollars 
BIST_new_scaling_value = 3*24*0.264*asp_out_fin(1,39)/asp_out_fin(1,22);
BIST_install_factor = 2.95; 
BIST_scaling_exp = 0.65; 
BIST_size_ratio = BIST_new_scaling_value/BIST_scaling_value;
BIST_scaled_cost = BIST_price*BIST_size_ratio^BIST_scaling_exp; 
BIST_scaled_cost_proj_yr = BIST_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(16,2)); 
BIST_instal_cost_proj_yr = BIST_scaled_cost_proj_yr*BIST_install_factor;

%Naphtha intermediate storage (3 days)
NIST_scaling_value = 558000; %gallons
NIST_price = 320384; %2005 dollars 
NIST_new_scaling_value = 3*24*0.264*asp_out_fin(1,44)/naphtha_dens;
NIST_install_factor = 2.95; 
NIST_scaling_exp = 0.65; 
NIST_size_ratio = NIST_new_scaling_value/NIST_scaling_value;
NIST_scaled_cost = NIST_price*NIST_size_ratio^NIST_scaling_exp; 
NIST_scaled_cost_proj_yr = NIST_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(16,2)); 
NIST_instal_cost_proj_yr = NIST_scaled_cost_proj_yr*NIST_install_factor;

%Biodiesel intermediate storage (3 days)
BDIST_scaling_value = 558000; %gallons
BDIST_price = 320384; %2005 dollars 
BDIST_new_scaling_value = 3*24*0.264*asp_out_fin(1,1)/diesel_dens;
BDIST_install_factor = 2.95; 
BDIST_scaling_exp = 0.65; 
BDIST_size_ratio = BDIST_new_scaling_value/BDIST_scaling_value;
BDIST_scaled_cost = BDIST_price*BDIST_size_ratio^BDIST_scaling_exp; 
BDIST_scaled_cost_proj_yr = BDIST_scaled_cost*(CAPEX_INDEX(CAP_IND_POS,2)/CAPEX_INDEX(16,2)); 
BDIST_instal_cost_proj_yr = BDIST_scaled_cost_proj_yr*BDIST_install_factor;

OTHER_UT_INST_EQUIP_CAPEX = CTS_instal_cost_proj_yr + CWP_instal_cost_proj_yr + PAC_instal_cost_proj_yr + HDT_instal_cost_proj_yr + FWP_instal_cost_proj_yr + IAD_instal_cost_proj_yr + PAR_instal_cost_proj_yr + FWST_instal_cost_proj_yr + BIST_instal_cost_proj_yr + NIST_instal_cost_proj_yr + BDIST_instal_cost_proj_yr; 

%% Total ISBL CAPEX, Direct Costs, and Indirect Costs
Total_ISBL = HTL_INST_EQUIP_CAPEX + UG_INST_EQUIP_CAPEX + HEAT_UT_INST_EQUIP_CAPEX + OTHER_UT_INST_EQUIP_CAPEX;

%Direct costs
warehouse = 0.04*Total_ISBL; 
site_devel = 0.09*Total_ISBL; 
add_piping = 0.045*Total_ISBL; 
HTL_UG_TDC = warehouse + site_devel + add_piping + Total_ISBL; 

%Indirect costs
prorate_exp = 0.10*HTL_UG_TDC; 
field_exp = 0.10*HTL_UG_TDC;
home_off_const = 0.20*HTL_UG_TDC;
contingency = 0.10*HTL_UG_TDC;
other_costs = 0.10*HTL_UG_TDC;

HTL_UG_TINC = prorate_exp + field_exp + home_off_const + contingency + other_costs;

%Fixed capital investment
HTL_UG_FCI = HTL_UG_TDC + HTL_UG_TINC; 
working_cap = 0.05*HTL_UG_FCI; 

%Total capital investment (TCI)
HTL_UP_TCI = HTL_UG_FCI + working_cap; 

%% Operational Expenses

%Variable OPEX - Cost of consumables
H2_OPEX = annual_hydrogen * H_cost * (CHEM_INDEX(CHEM_IND_POS,2)/CHEM_INDEX(33,2));
NG_OPEX = annual_NG * NG_cost * (CHEM_INDEX(CHEM_IND_POS,2)/CHEM_INDEX(32,2));
ELEC_OPEX = annual_elec * Elec_cost * (ELECT_INDEX(ELEC_IND_POS,2)/ELECT_INDEX(6,2));
PW_OPEX = annual_water * Process_water_cost * (CHEM_INDEX(CHEM_IND_POS,2)/CHEM_INDEX(32,2));

Total_var_OPX = H2_OPEX + NG_OPEX + ELEC_OPEX + PW_OPEX; 

%Co-product revenue
Revenue_char = annual_HTL_solids*Char_recovery*Char_value*(CHEM_INDEX(CHEM_IND_POS,2)/CHEM_INDEX(40,2));
Revenue_NH3EQ = annual_NH3EQ*NH3EQ_recovery*NH3EQ_value*(1-nut_rem)*(CHEM_INDEX(CHEM_IND_POS,2)/CHEM_INDEX(40,2)); 
Revenue_DAPEQ = annual_DAPEQ*DAPEQ_recovery*DAPEQ_value*(1-nut_rem)*(CHEM_INDEX(CHEM_IND_POS,2)/CHEM_INDEX(40,2));

%Nutrient Removal Revenue
N_removal_rev = cult_outs(8,1)*0.82*nut_rem*15;
P_removal_rev = cult_outs(9,1)*0.20*nut_rem*100;

Total_co_prod_rev = Revenue_char + Revenue_NH3EQ + Revenue_DAPEQ + N_removal_rev + P_removal_rev; 

%Fixed OPEX - labor, burdens, maintenance, insurance, taxes

%Labor - all salaries quoted in 2009 USD
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

labor_total = (fac_size/5000)*(plant_mngr + plant_engr + maint_sup + maint_tech + lab_mngr + lab_tech + shft_sup + shft_op + yrd_emp + clk_sec);
labor_burden = labor_total * 0.90; 
maintenance = 0.03 * Total_ISBL; 
prop_ins_tax = 0.01 * HTL_UG_FCI;

Total_fixed_OPX = labor_total + labor_burden + maintenance + prop_ins_tax;

%% Output final CAPEX/OPEX/Revenue matrix 
HTL_UG_CAPX_OPX(1,1) = "HTL Equipment"; 
HTL_UG_CAPX_OPX(1,2) = "Installed Cost in Project Year";
HTL_UG_CAPX_OPX(2,1) = "Purge/Reactor Exchanger"; 
HTL_UG_CAPX_OPX(2,2) = PR_instal_cost_proj_yr;
HTL_UG_CAPX_OPX(3,1) = "HTL Reactor";
HTL_UG_CAPX_OPX(3,2) = HTL_instal_cost_proj_yr;
HTL_UG_CAPX_OPX(4,1) = "Reactor Gas KO Drum"; 
HTL_UG_CAPX_OPX(4,2) = KO_instal_cost_proj_yr;
HTL_UG_CAPX_OPX(5,1) = "Solids Filter";
HTL_UG_CAPX_OPX(5,2) = SF_instal_cost_proj_yr;
HTL_UG_CAPX_OPX(6,1) = "Separator"; 
HTL_UG_CAPX_OPX(6,2) = SEP_instal_cost_proj_yr; 
HTL_UG_CAPX_OPX(7,1) = "Purge Water Cooler"; 
HTL_UG_CAPX_OPX(7,2) = PWC_instal_cost_proj_yr;
HTL_UG_CAPX_OPX(8,1) = "HTL Total"; 
HTL_UG_CAPX_OPX(8,2) = HTL_INST_EQUIP_CAPEX;
HTL_UG_CAPX_OPX(9,1) = "Upgrading Equipment"; 
HTL_UG_CAPX_OPX(9,2) = "Installed Cost in Project Year";
HTL_UG_CAPX_OPX(10,1) = "HT Reactor/Vessels/Columns"; 
HTL_UG_CAPX_OPX(10,2) = HT_instal_cost_proj_yr;
HTL_UG_CAPX_OPX(11,1) = "Hydrogen Compressor"; 
HTL_UG_CAPX_OPX(11,2) = HYD_COMP_instal_cost_proj_yr;
HTL_UG_CAPX_OPX(12,1) = "Pressure Swing Absorber";
HTL_UG_CAPX_OPX(12,2) = PSA_instal_cost_proj_yr;
HTL_UG_CAPX_OPX(13,1) = "Hydrocracker and Auxiliaries";
HTL_UG_CAPX_OPX(13,2) = HC_instal_cost_proj_yr;
HTL_UG_CAPX_OPX(14,1) = "Upgrading Total";
HTL_UG_CAPX_OPX(14,2) = UG_INST_EQUIP_CAPEX;
HTL_UG_CAPX_OPX(15,1) = "Heating Utility"; 
HTL_UG_CAPX_OPX(15,2) = "Installed Cost in Project Year";
HTL_UG_CAPX_OPX(16,1) = "Reactor Heater";
HTL_UG_CAPX_OPX(16,2) = RH_instal_cost_proj_yr;
HTL_UG_CAPX_OPX(17,1) = "Bio-Oil Heat Recovery Steam Generator"; 
HTL_UG_CAPX_OPX(17,2) = BOHR_instal_cost_proj_yr;
HTL_UG_CAPX_OPX(18,1) = "Hot Oil System Package";
HTL_UG_CAPX_OPX(18,2) = HOS_instal_cost_proj_yr;
HTL_UG_CAPX_OPX(19,1) = "Hot Oil"; 
HTL_UG_CAPX_OPX(19,2) = HO_instal_cost_proj_yr;
HTL_UG_CAPX_OPX(20,1) = "Heating Utility Total";
HTL_UG_CAPX_OPX(20,2) = HEAT_UT_INST_EQUIP_CAPEX;
HTL_UG_CAPX_OPX(21,1) = "Other Utilities"; 
HTL_UG_CAPX_OPX(21,2) = "Installed Cost in Project Year";
HTL_UG_CAPX_OPX(22,1) = "Cooling Tower System"; 
HTL_UG_CAPX_OPX(22,2) = CTS_instal_cost_proj_yr;
HTL_UG_CAPX_OPX(23,1) = "Cooling Water Pump"; 
HTL_UG_CAPX_OPX(23,2) = CWP_instal_cost_proj_yr;
HTL_UG_CAPX_OPX(24,1) = "Plant Air Compressor";
HTL_UG_CAPX_OPX(24,2) = PAC_instal_cost_proj_yr;
HTL_UG_CAPX_OPX(25,1) = "Hydraulic Dump Truck";
HTL_UG_CAPX_OPX(25,2) = HDT_instal_cost_proj_yr;
HTL_UG_CAPX_OPX(26,1) = "Firewater Pump";
HTL_UG_CAPX_OPX(26,2) = FWP_instal_cost_proj_yr;
HTL_UG_CAPX_OPX(27,1) = "Instrument Air Dryer";
HTL_UG_CAPX_OPX(27,2) = IAD_instal_cost_proj_yr;
HTL_UG_CAPX_OPX(28,1) = "Plant Air Receiver";
HTL_UG_CAPX_OPX(28,2) = PAR_instal_cost_proj_yr;
HTL_UG_CAPX_OPX(29,1) = "Firewater Storage Tank";
HTL_UG_CAPX_OPX(29,2) = FWST_instal_cost_proj_yr;
HTL_UG_CAPX_OPX(30,1) = "HTL Oil Intermediate Storage (3 Days)";
HTL_UG_CAPX_OPX(30,2) = BIST_instal_cost_proj_yr;
HTL_UG_CAPX_OPX(31,1) = "Naphtha Intermediate Storage (3 Days)";
HTL_UG_CAPX_OPX(31,2) = NIST_instal_cost_proj_yr;
HTL_UG_CAPX_OPX(32,1) = "Biodiesel Intermediate Storage (3 Days)";
HTL_UG_CAPX_OPX(32,2) = BDIST_instal_cost_proj_yr;
HTL_UG_CAPX_OPX(33,1) = "Other Utilities Total";
HTL_UG_CAPX_OPX(33,2) = OTHER_UT_INST_EQUIP_CAPEX;
HTL_UG_CAPX_OPX(34,1) = "Total Inside Battery Limits";
HTL_UG_CAPX_OPX(34,2) = Total_ISBL;
HTL_UG_CAPX_OPX(35,1) = "Warehouse";
HTL_UG_CAPX_OPX(35,2) = warehouse;
HTL_UG_CAPX_OPX(36,1) = "Site Development";
HTL_UG_CAPX_OPX(36,2) = site_devel;
HTL_UG_CAPX_OPX(37,1) = "Additional Piping";
HTL_UG_CAPX_OPX(37,2) = add_piping; 
HTL_UG_CAPX_OPX(38,1) = "HTL and UG Total Direct Capital Costs";
HTL_UG_CAPX_OPX(38,2) = HTL_UG_TDC;
HTL_UG_CAPX_OPX(39,1) = "Proratable Expenses";
HTL_UG_CAPX_OPX(39,2) = prorate_exp;
HTL_UG_CAPX_OPX(40,1) = "Field Expenses";
HTL_UG_CAPX_OPX(40,2) = field_exp;
HTL_UG_CAPX_OPX(41,1) = "Home Office and Construction Feed";
HTL_UG_CAPX_OPX(41,2) = home_off_const;
HTL_UG_CAPX_OPX(42,1) = "Project Contingency";
HTL_UG_CAPX_OPX(42,2) = contingency;
HTL_UG_CAPX_OPX(43,1) = "Other Costs";
HTL_UG_CAPX_OPX(43,2) = other_costs;
HTL_UG_CAPX_OPX(44,1) = "Total Indirect Costs";
HTL_UG_CAPX_OPX(44,2) = HTL_UG_TINC;
HTL_UG_CAPX_OPX(45,1) = "Fixed Capital Investment";
HTL_UG_CAPX_OPX(45,2) = HTL_UG_FCI;
HTL_UG_CAPX_OPX(46,1) = "Working Capital";
HTL_UG_CAPX_OPX(46,2) = working_cap;
HTL_UG_CAPX_OPX(47,1) = "Total Capital Investment";
HTL_UG_CAPX_OPX(47,2) = HTL_UP_TCI;
HTL_UG_CAPX_OPX(48,1) = "Variable OPEX";
HTL_UG_CAPX_OPX(48,2) = "Annual OPEX"; 
HTL_UG_CAPX_OPX(49,1) = "H2 OPEX";
HTL_UG_CAPX_OPX(49,2) = H2_OPEX;
HTL_UG_CAPX_OPX(50,1) = "Natural Gas OPEX";
HTL_UG_CAPX_OPX(50,2) = NG_OPEX;
HTL_UG_CAPX_OPX(51,1) = "Grid Electricity";
HTL_UG_CAPX_OPX(51,2) = ELEC_OPEX;
HTL_UG_CAPX_OPX(52,1) = "Process Water";
HTL_UG_CAPX_OPX(52,2) = PW_OPEX;
HTL_UG_CAPX_OPX(53,1) = "Total Variable OPEX";
HTL_UG_CAPX_OPX(53,2) = Total_var_OPX;
HTL_UG_CAPX_OPX(54,1) = "Co-Product Revenue";
HTL_UG_CAPX_OPX(54,2) = "Annual Revenue";
HTL_UG_CAPX_OPX(55,1) = "Char Revenue"; 
HTL_UG_CAPX_OPX(55,2) = Revenue_char;
HTL_UG_CAPX_OPX(56,1) = "NH3 Revenue";
HTL_UG_CAPX_OPX(56,2) = Revenue_NH3EQ;
HTL_UG_CAPX_OPX(57,1) = "DAP Revenue"; 
HTL_UG_CAPX_OPX(57,2) = Revenue_DAPEQ;
HTL_UG_CAPX_OPX(58,1) = "Total Co-Product Revenue";
HTL_UG_CAPX_OPX(58,2) = Total_co_prod_rev;
HTL_UG_CAPX_OPX(59,1) = "Fixed OPEX"; 
HTL_UG_CAPX_OPX(59,2) = "Annual OPEX";
HTL_UG_CAPX_OPX(60,1) = "Total Labor"; 
HTL_UG_CAPX_OPX(60,2) = labor_total;
HTL_UG_CAPX_OPX(61,1) = "Labor Burden";
HTL_UG_CAPX_OPX(61,2) = labor_burden;
HTL_UG_CAPX_OPX(62,1) = "Maintenance";
HTL_UG_CAPX_OPX(62,2) = maintenance;
HTL_UG_CAPX_OPX(63,1) = "Property Insurance and Tax";
HTL_UG_CAPX_OPX(63,2) = prop_ins_tax;
HTL_UG_CAPX_OPX(64,1) = "Total Fixed OPEX"; 
HTL_UG_CAPX_OPX(64,2) = Total_fixed_OPX;

end
