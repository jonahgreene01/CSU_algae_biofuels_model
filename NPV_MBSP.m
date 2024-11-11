function [NPV, cash_flow] = NPV_MBSP(MBSP, total_tonnes_AFDW, land_CAPEX, equip_fac_CAPEX, working_CAPEX, total_OPEX, IRR)

% Test Inputs
% total_tonnes_AFDW = 188410.14; 
% equip_fac_CAPEX = 457853520;
% land_CAPEX = 24428494; 
% working_CAPEX = 0.05*equip_fac_CAPEX; 
% total_OPEX = 59150012; 
% IRR = 0.10; 
%MBSP = 677.97; 

cost_cap = equip_fac_CAPEX; 
cost_op = total_OPEX; 
full_prod = total_tonnes_AFDW; 
cost_rev = full_prod*MBSP; 

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
cash_flow = table('Size',sz,'VariableTypes',var_types,'VariableNames',var_names); % Initialize the Cash Flow Table

%% Calculations

cash_flow.Year(:) = -con_term+1:1:time_horizon; % Fill in the year column

% Calculate Discount Rates
    for i = 1:size(cash_flow,1)
        cash_flow.DiscountFactor(i) = 1/(1+irr)^cash_flow.Year(i);
    end

% Calculations for Construction Years
for i = 1:con_term
    cash_flow.Capital(i) = cost_cap*equity_amt*con_finish(i);
    if isequal(i,1)
        cash_flow.LoanPrinciple(i) = cost_cap*con_finish(i)*loan_amt;
        cash_flow.LandCapital(i) = land_CAPEX; 
    else
        cash_flow.LoanPrinciple(i) = cost_cap*con_finish(i)*loan_amt+cash_flow.LoanPrinciple(i-1);
    end
    
    if isequal(i,3)
        cash_flow.WorkingCapital(i) = working_CAPEX;
    end
    
    cash_flow.InterestPayment(i) = cash_flow.LoanPrinciple(i)*int_loan;
    cash_flow.NPVCapitalPlusInterest(i) = (cash_flow.Capital(i) + cash_flow.WorkingCapital(i) + cash_flow.InterestPayment(i))*cash_flow.DiscountFactor(i);
    cash_flow.NPVCapitalCosts(i)=cash_flow.DiscountFactor(i)*cash_flow.Capital(i);
    cash_flow.NPVLand(i)=cash_flow.DiscountFactor(i)*cash_flow.LandCapital(i);
end

% Calculations for Operating Years

loan_pmt = (cost_cap*int_loan*loan_amt)/(1-(1+int_loan)^(-loan_term)); % Calculate Annual Loan Payment
cash_flow.LoanPayment(con_term+1:loan_term+con_term) = loan_pmt;

for i = con_term+1:loan_term + con_term
    cash_flow.InterestPayment(i) = cash_flow.LoanPrinciple(i-1)*int_loan;
    cash_flow.LoanPrinciple(i) = cash_flow.LoanPrinciple(i-1)-cash_flow.LoanPayment(i)+cash_flow.InterestPayment(i);
end

for i = con_term+1:time_horizon+con_term
        if isequal(i,4)
            cash_flow.Revenue(i) = (full_prod*MBSP)*(1-start_up_time)+(full_prod*MBSP)*start_up_time*start_up_prod;
            cash_flow.OperationCosts(i) = (cost_op)*(1-start_up_time)+(cost_op)*start_up_time*start_up_prod;
        else
            cash_flow.Revenue(i) = cost_rev;
            cash_flow.OperationCosts(i) = cost_op;
        end
    
end
    
cash_flow.DepreciationRate(con_term+1:con_term+size(dep_rate,2)) = dep_rate;
cash_flow.CapitalDepreciation(con_term+1:con_term+size(dep_rate,2)) = cost_cap*dep_rate;

for i = con_term+1:time_horizon+con_term
    cash_flow.NetRevenue(i) = cash_flow.Revenue(i)-cash_flow.LoanPayment(i)-cash_flow.OperationCosts(i)-cash_flow.CapitalDepreciation(i);
    
    if cash_flow.TaxableIncome(i-1) < 0
        cash_flow.LossesForward(i) = cash_flow.TaxableIncome(i-1);
    end

    cash_flow.TaxableIncome(i) = cash_flow.NetRevenue(i)+cash_flow.LossesForward(i);
    
    if cash_flow.TaxableIncome(i)*tax_rate > 0
        cash_flow.IncomeTax(i) = cash_flow.TaxableIncome(i)*tax_rate;
    end
        
    cash_flow.CashIncome(i) = cash_flow.Revenue(i)-cash_flow.LoanPayment(i)-cash_flow.OperationCosts(i)-cash_flow.IncomeTax(i);
    cash_flow.PresentValue(i) = cash_flow.CashIncome(i)*cash_flow.DiscountFactor(i);  
    
    cash_flow.NPVTax(i)=cash_flow.DiscountFactor(i)*cash_flow.IncomeTax(i);
    cash_flow.NPVRevenue(i)=cash_flow.DiscountFactor(i)*cash_flow.Revenue(i);
    cash_flow.NPVOperationCosts(i)=cash_flow.DiscountFactor(i)*cash_flow.OperationCosts(i);
    
end
   
NPV = sum(cash_flow.PresentValue(:)) - (sum(cash_flow.NPVCapitalPlusInterest(:))+sum(cash_flow.NPVLand(:)));

end 