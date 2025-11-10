
% C(n)H(2n+2) + (n/4)O2 + (n/2)H2O <-> (n)CO + ((n*3/2)+1)H2
% C4.36H10.72 + (4.36/4)O2 + (4.36/2)H2O ---> 4.36CO + 7.54H2

function [e_ATR,b_ATR,c_ATR,d_ATR] = Comidy_ATR_FO(mol_C2to8,e_prod)

f = mol_C2to8;
    
% Molar balance
e_ATR = (4.36/4)*f;
b_ATR = (4.36/2)*f; 
c_ATR = 4.36*f; 
d_ATR = 7.54*f;
    

%Checking O2 and H2O requirement:
fprintf('e_ATR: %f\n',e_ATR)
fprintf('e_prod: %f\n',e_prod)
O2_count = e_ATR - e_prod;
fprintf('O2_count: %f',O2_count)
if O2_count > 0 
   fprintf (' O2 required for ATR is insufficient: %f required', O2_count);
else
   fprintf ('O2 for ATR is sufficient: %f leftover',O2_count);
end 

fprintf ('H2O required for ATR is: %f mols',b_ATR);
disp('ATR run complete')
end
