

function [mol_C2to8,mass_c_FT,mass_d_FT,mass_b_FT,FT_C2to8,FT_C9to16,FT_C17plus,...
mm_C,mm_H,mass_c_unreact] = Comidy_FT_FO(mass_c_prod,mass_d_prod,mm_b,mm_c,mm_d)

%Definining molar masses of carbon and hydrogen
mm_C = 12.0107 * (1/1000); %kg/mol
mm_H = 1.00794 * (1/1000); 


 % Reaction based on 1 mol of CO:

 mol_c_prod = mass_c_prod/mm_c; 
 mol_c_react = 0.75*mol_c_prod; % Conversion factor
 mol_c_unreact = 0.25*mol_c_prod;
 mol_C2to8 = 0.089*mol_c_react;
 mol_C9to16 = 0.0304*mol_c_react;
 mol_C17plus = 0.0114*mol_c_react;
        
 mol_d_consumed = 2.15 * mol_c_react;
 
 mol_b_produced = mol_c_react;

%Convert to mass (kg) using average molecular weights
mass_c_FT = mol_c_react * mm_c;
mass_c_unreact = mol_c_unreact * mm_c;
mass_d_FT = mol_d_consumed * mm_d; % kg
mass_b_FT = mol_b_produced * mm_b; % kg
FT_C2to8 = mol_C2to8 * ((4.36*mm_C)+(10.72*mm_H));
FT_C9to16 = mol_C9to16 * ((11.67*mm_C)+(25.34*mm_H));
FT_C17plus = mol_C17plus * ((22.67*mm_C)+(47.34*mm_H));

%Tracking the H2 requirement; 

    if mass_d_prod < mass_d_FT
       H2_count = mass_d_FT - mass_d_prod;
        fprintf('Need %f more kg H2 for full synthesis.\n\n', H2_count);
    end 

      
    %Display the amount of unreacted feed (15% of the reactant)
    fprintf ('Amount of unreacted CO is %f mol\n', mol_c_unreact);

end

