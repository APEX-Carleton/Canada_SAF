
function [user_CO2,user_temp,e_prod,mass_b_feed,mass_a_recy,mass_b_recy,mass_c_recy,mass_d_recy,mass_a_out,...
   mass_b_out,mass_c_out,mass_d_out,mass_a_prod,mass_b_prod,mass_c_prod,mass_d_prod,mass_e,...
    mm_a,mm_b,mm_c,mm_d,mm_e,c_prod,d_prod,elec_E,heat_E,syngas_out]... 
    = Co_Electrolysis_FO(user_CO2,user_temp)

% Co-Electrolysis Process: 

% C02 + 2.12H20 -> CO + 2.12H2 + 1.56O2


% a = CO2, b = H2O, c = CO, d = H2, e = O2
% Molar masses
mm_a = 44.0095 * (1/1000); % kg/mol
mm_b = 18.01528 * (1/1000);
mm_c = 28.0101 * (1/1000);
mm_d = 2.01588 * (1/1000);
mm_e = 31.9988 * (1/1000);




% Assigning stoichiometric coefficients 
ratio_b = 2.12; ratio_e = 1.56;



%Defining inlet conditions
a_feed = user_CO2 / mm_a; % in mol
b_feed = ratio_b*a_feed;
feedstock = a_feed + b_feed;
tl_mol_in = feedstock/0.90    
recycled_feedgas = 0.01*tl_mol_in;
feedstock_in = feedstock + recycled_feedgas


% Chemical reaction
feedstock_consumed = 0.9*feedstock_in
feedstock_out = feedstock_in-feedstock_consumed
a_consumed = feedstock_consumed/3.12
b_consumed = feedstock_consumed-a_consumed
c_out = a_consumed
d_out = ratio_b*a_consumed
e_prod = ratio_e * a_consumed

% Outlet conditions
syngas_out = c_out +d_out
recycled_syngas = 0.1*tl_mol_in
c_recy = recycled_syngas/3.12
d_recy = recycled_syngas - c_recy
a_out = feedstock_out/3.12
b_out = feedstock_out - a_out
a_recy = 0.111*a_out
b_recy = 0.111*b_out

% Product Stream components
a_prod = a_out - a_recy
b_prod = b_out - b_recy
c_prod = c_out - c_recy
d_prod = d_out - d_recy

% CHECKS

if  (0.05*recycled_feedgas+recycled_feedgas) >= a_recy + b_recy && (recycled_feedgas-0.02*recycled_feedgas) <= a_recy + b_recy
    disp ('error in feedgas recyle')
end 
    
if recycled_syngas ~= 0.1*tl_mol_in
    disp('error in syngas recycle. Calculation at outlet does not match inlet stream fractions')
end 

if e_prod ~= 0.41*tl_mol_in
    disp ('Error. check Oxygen')
end 

if tl_mol_in ~= feedstock + recycled_feedgas + recycled_syngas
    disp('Error in stream fractions.')
end 

% Calculating component masses:
mass_b_feed =  b_feed * mm_b; 
mass_a_recy = a_recy * mm_a;
mass_b_recy = b_recy * mm_b;
mass_c_recy = c_recy * mm_c;
mass_d_recy = d_recy * mm_d;
mass_a_out = a_out * mm_a;
mass_b_out = b_out * mm_b;
mass_c_out = c_out * mm_c;
mass_d_out = d_out * mm_d;
mass_a_prod = a_prod * mm_a;
mass_b_prod = b_prod * mm_b;
mass_c_prod = c_prod * mm_c;
mass_d_prod = d_prod * mm_d;
mass_e = e_prod * mm_e;




%Electricity demand kJ/mol 
m_E = -0.14 %slope
b_E = 460 %y-int
elec_sp_E = (m_E)*(user_temp) + b_E % kJ/mol
elec_E = elec_sp_E * (syngas_out)/3600 % kWh

%Heat demand kJ/mol 
m_H = 0.14 %slope
b_H = 70 %y-int
heat_sp_E = (m_H)*(user_temp) + b_H % kJ/mol
heat_E = heat_sp_E * (syngas_out)/3600 % kWh
%end