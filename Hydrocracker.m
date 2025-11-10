
function [mass_gases,mass_mid_distillates,mass_heavy_cuts,HC_res] = Hydrocracking_FO(FT_C17plus,mm_H)

% Product fractions as wt% of reactant
mass_gases = 0.2215 * FT_C17plus; 
mass_mid_distillates = 0.6102 * FT_C17plus;
mass_heavy_cuts = 0.1684 * FT_C17plus;

% Light ends to be sent to ATR
HC_res = mass_gases;

%H2 requirement
mass_H2_HC = 0.031*FT_C17plus; %kg

%Heat released: 42 KJ/mol of H2 consumed:
Heat_HC = 42 * (mass_H2_HC/(2*mm_H));
end
