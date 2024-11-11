function [NPV_mfsp, cash_flow_mfsp] = NPV_MFSP(MFSP, GGE_per_yr, land_CAPEX, total_equip_fac_CAPEX, total_working_CAPEX, total_combined_OPEX, IRR, total_coprod_rev)

% Test Inputs
% total_tonnes_AFDW = 188410.14; 
% equip_fac_CAPEX = 457853520;
% land_CAPEX = 24428494; 
% working_CAPEX = 0.05*equip_fac_CAPEX; 
% total_OPEX = 59150012; 
% IRR = 0.10; 
%MFSP = 677.97; 

cost_cap = total_equip_fac_CAPEX; 
cost_op = total_combined_OPEX; 
full_prod = GGE_per_yr; 
cost_rev = full_prod*MFSP; 

% Financial Assumptions
% Time horizon and start-up period
time_horizon = 30; %years
start_up_time = 0.5; 
start_up_prod = 0.5; 

% Loan Assumptions
int_loan = 0.08; % Loan Interest Rate
loan_term = 10; % Term of the loan (years)
loan_amt = 0.60; % Amount of capital financed (%)
equity_amt = 1-loan_amt; % Amount of capital equity (%)   

% Construction Assumptions
con_term = 3; % Construction timeframe (yrs)
con_finish = [0.08,0.60,0.32]; % Percent of construction finished in years -2 to 0 (%)

% Depreciation Assumptions
dep_rate = [14.29, 24.49, 17.49, 12.49, 8.93, 8.92, 8.93, 4.46]/100; % MACRS Depreciation scheme

% Tax Assumptions
tax_rate = 0.35; % Federal + State Tax Rate (%)

% Return Assumptions
irr = IRR; % Internal Rate of Return (%)   

%% Initialize Tables
var_names = ["Year","Capital", "LandCapital", "WorkingCapital", "LoanPayment","InterestPayment","LoanPrinciple","Revenue","OperationCosts","DepreciationRate","CapitalDepreciation","NetRevenue","LossesForward","TaxableIncome","IncomeTax","CashIncome","DiscountFactor","PresentValue","NPVCapitalPlusInterest","NPVTax","NPVRevenue","NPVOperationCosts","NPVCapitalCosts", "NPVLand"];
sz = [con_term+time_horizon,size(var_names,2)]; % Size of the Cash Flow Table
var_types (1:sz(2)) = "double"; % Variable Types in the Cash Flow Table
cash_flow_mfsp = table('Size',sz,'VariableTypes',var_types,'VariableNames',var_names); % Initialize the Cash Flow Table

%% Calculations

cash_flow_mfsp.Year(:) = -con_term+1:1:time_horizon; % Fill in the year column

% Calculate Discount Rates
    for i = 1:size(cash_flow_mfsp,1)
        cash_flow_mfsp.DiscountFactor(i) = 1/(1+irr)^cash_flow_mfsp.Year(i);
    end

% Calculations for Construction Years
for i = 1:con_term
    cash_flow_mfsp.Capital(i) = cost_cap*equity_amt*con_finish(i);
    if isequal(i,1)
        cash_flow_mfsp.LoanPrinciple(i) = cost_cap*con_finish(i)*loan_amt;
        cash_flow_mfsp.LandCapital(i) = land_CAPEX; 
    else
        cash_flow_mfsp.LoanPrinciple(i) = cost_cap*con_finish(i)*loan_amt+cash_flow_mfsp.LoanPrinciple(i-1);
    end
    
    if isequal(i,3)
        cash_flow_mfsp.WorkingCapital(i) = total_working_CAPEX;
    end
    
    cash_flow_mfsp.InterestPayment(i) = cash_flow_mfsp.LoanPrinciple(i)*int_loan;
    cash_flow_mfsp.NPVCapitalPlusInterest(i) = (cash_flow_mfsp.Capital(i) + cash_flow_mfsp.WorkingCapital(i) + cash_flow_mfsp.InterestPayment(i))*cash_flow_mfsp.DiscountFactor(i);
    cash_flow_mfsp.NPVCapitalCosts(i)=cash_flow_mfsp.DiscountFactor(i)*cash_flow_mfsp.Capital(i);
    cash_flow_mfsp.NPVLand(i)=cash_flow_mfsp.DiscountFactor(i)*cash_flow_mfsp.LandCapital(i);
end

% Calculations for Operating Years

loan_pmt = (cost_cap*int_loan*loan_amt)/(1-(1+int_loan)^(-loan_term)); % Calculate Annual Loan Payment
cash_flow_mfsp.LoanPayment(con_term+1:loan_term+con_term) = loan_pmt;

for i = con_term+1:loan_term + con_term
    cash_flow_mfsp.InterestPayment(i) = cash_flow_mfsp.LoanPrinciple(i-1)*int_loan;
    cash_flow_mfsp.LoanPrinciple(i) = cash_flow_mfsp.LoanPrinciple(i-1)-cash_flow_mfsp.LoanPayment(i)+cash_flow_mfsp.InterestPayment(i);
end

for i = con_term+1:time_horizon+con_term
        if isequal(i,4)
            cash_flow_mfsp.Revenue(i) = ((full_prod*MFSP)*(1-start_up_time)+(full_prod*MFSP)*start_up_time*start_up_prod) + (total_coprod_rev*(1-start_up_time) + total_coprod_rev*start_up_time*start_up_prod);
            cash_flow_mfsp.OperationCosts(i) = (cost_op)*(1-start_up_time)+(cost_op)*start_up_time*start_up_prod;
        else
            cash_flow_mfsp.Revenue(i) = cost_rev + total_coprod_rev;
            cash_flow_mfsp.OperationCosts(i) = cost_op;
        end
    
end
    
cash_flow_mfsp.DepreciationRate(con_term+1:con_term+size(dep_rate,2)) = dep_rate;
cash_flow_mfsp.CapitalDepreciation(con_term+1:con_term+size(dep_rate,2)) = cost_cap*dep_rate;

for i = con_term+1:time_horizon+con_term
    cash_flow_mfsp.NetRevenue(i) = cash_flow_mfsp.Revenue(i)-cash_flow_mfsp.LoanPayment(i)-cash_flow_mfsp.OperationCosts(i)-cash_flow_mfsp.CapitalDepreciation(i);
    
    if cash_flow_mfsp.TaxableIncome(i-1) < 0
        cash_flow_mfsp.LossesForward(i) = cash_flow_mfsp.TaxableIncome(i-1);
    end

    cash_flow_mfsp.TaxableIncome(i) = cash_flow_mfsp.NetRevenue(i)+cash_flow_mfsp.LossesForward(i);
    
    if cash_flow_mfsp.TaxableIncome(i)*tax_rate > 0
        cash_flow_mfsp.IncomeTax(i) = cash_flow_mfsp.TaxableIncome(i)*tax_rate;
    end
        
    cash_flow_mfsp.CashIncome(i) = cash_flow_mfsp.Revenue(i)-cash_flow_mfsp.LoanPayment(i)-cash_flow_mfsp.OperationCosts(i)-cash_flow_mfsp.IncomeTax(i);
    cash_flow_mfsp.PresentValue(i) = cash_flow_mfsp.CashIncome(i)*cash_flow_mfsp.DiscountFactor(i);  
    
    cash_flow_mfsp.NPVTax(i)=cash_flow_mfsp.DiscountFactor(i)*cash_flow_mfsp.IncomeTax(i);
    cash_flow_mfsp.NPVRevenue(i)=cash_flow_mfsp.DiscountFactor(i)*cash_flow_mfsp.Revenue(i);
    cash_flow_mfsp.NPVOperationCosts(i)=cash_flow_mfsp.DiscountFactor(i)*cash_flow_mfsp.OperationCosts(i);
    
end
   
NPV_mfsp = sum(cash_flow_mfsp.PresentValue(:)) - (sum(cash_flow_mfsp.NPVCapitalPlusInterest(:)) + sum(cash_flow_mfsp.NPVLand(:)));

end 