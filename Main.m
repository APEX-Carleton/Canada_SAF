
%Saving outputs to designated file-Prep
date=string(datetime('now','Format',"yyyyMMddHHmmss"));


month_min=85;
startingYear=2019;
refmonth=0;
scenario_Array=["2.50","3.73","2.35","3.10"];
for idx=1
    user_select=scenario_Array(idx);
    
    % Run ATF demand module 
    [monthly_CO2_cons,monthly_ATF_cons,month_prov_ATF,CO2_mass_yearly,ATF_mass_yearly,month_min,month_max_ref,startingYear,user_select] = fuel_demand_FO(user_select,startingYear,month_min);
   
    additional_months=month_max_ref-month_min;
    month_max = month_min + additional_months;
    
    % Sizing the output tables
    header=[2011+1:2011+month_max_ref/12];
    outputcolumnidx=startingYear-2011;
    % Hydrocarbons produced, kg
    Table_sum_useful=monthly_ATF_cons(:,:);
    Table_sum_heavy_sell= monthly_ATF_cons(:,:);
    Table_sum_heavy= monthly_ATF_cons(:,:);
    Table_sum_light_sell= monthly_ATF_cons(:,:);
    % inputs required for the process, kg
    Table_H2_requirement= monthly_ATF_cons(:,:);
    Table_H2O_requirement=monthly_ATF_cons(:,:);
    Table_O2_requirement=monthly_ATF_cons(:,:);
    % Net input required for the process: 
    Table_H2_NetRequirement=monthly_ATF_cons(:,:);
    Table_H2O_NetRequirement=monthly_ATF_cons(:,:);
    Table_O2_NetRequirement=monthly_ATF_cons(:,:);
    % Mass balances
    Table_H2_NetRequired_FT=monthly_ATF_cons(:,:);
    Table_H2_Stock=monthly_ATF_cons(:,:);
    Table_H2O_Produced=monthly_ATF_cons(:,:);
    Table_H2_Produced=monthly_ATF_cons(:,:);
    Table_O2_Produced=monthly_ATF_cons(:,:);
    H2_Stock=0;
    

    %%Sizing unit input tables
    Table_CO2_SOEC_requirement= monthly_ATF_cons(1,:);
    Table_H2O_SOEC_requirement= monthly_ATF_cons(1,:);
    Table_Heat_SOEC_requirement= monthly_ATF_cons(1,:);
    Table_Elec_SOEC_requirement= monthly_ATF_cons(1,:);

    Table_HC_sizingFactor=monthly_ATF_cons(1,:);
    Table_FT_sizingFactor=monthly_ATF_cons(1,:);
    Table_ATR_sizingFactor=monthly_ATF_cons(1,:);
 
    %%

    % Reshape arrays
    local_monthly_CO2_cons = reshape(monthly_CO2_cons,1,numel(monthly_CO2_cons));
    local_ATF_cons = reshape(monthly_ATF_cons,1,numel(monthly_ATF_cons));
    
    % Electrolyzer temperature prompt
    user_temp = input('Enter an electrolyzer temperature between 100 and 1000C: ');
    temp_min = 100; temp_max = 1000; % deg C
    if user_temp < temp_min || user_temp > temp_max
        disp('Sorry, not a valid input! Please restart.');
        return
    end
    check = input('Would you like to keep this temperature constant throughout the specified time period? [Y/N] ', 's');
    if check == 'N' || check == 'n'
        elec_temp = true;
    elseif check == 'Y' || check == 'y'
        elec_temp = false;
    else
        disp('Sorry, not a valid input! Please restart.');
        return
    end
    
    % Array Pre-allocation
    tl_CO2_input_CoElec = zeros(1,additional_months+1); 
    tl_H2O_required_CoElec = zeros(1,additional_months+1); 
    tl_elec_input_CoElec = zeros(1,additional_months+1); 
    tl_heat_input_CoElec = zeros(1,additional_months+1); 
    tl_CO_output_CoElec = zeros(1,additional_months+1); 
    tl_H2_output_CoElec = zeros(1,additional_months+1); 
    tl_O2_output_CoElec = zeros(1,additional_months+1); 
    %First index is 6; the number of complete loops. 
    tl_CO_in_FT = zeros(6,additional_months+1); 
    tl_H2_in_FT = zeros(6,additional_months+1);% ''
    tl_C2to8_FT = zeros(6,additional_months+1); 
    tl_C9to16_FT = zeros(6,additional_months+1);
    tl_C17plus_FT = zeros(6,additional_months+1);
    tl_H2O_out_FT = zeros(6,additional_months+1);
    tl_H2_HC = zeros(6,additional_months+1);
    tl_gases_HC = zeros(6,additional_months+1);
    tl_naphs_HC = zeros(6,additional_months+1);
    tl_hvy_dtls_HC = zeros(6,additional_months+1);
    tl_lightEnds_in_ATR = zeros(6,additional_months+1);
    tl_O2_in_ATR = zeros(6,additional_months+1);
    tl_H2O_in_ATR = zeros(6,additional_months+1);
    tl_CO_out_ATR = zeros(6,additional_months+1);
    tl_H2_out_ATR = zeros(6,additional_months+1);
    tl_unreacted_CO_in_FT = zeros(6,additional_months+1); 
    tl_CO_to_FT = zeros(6,additional_months+1); 
 
    tl_Syngas_out_ATR=zeros(6,additional_months+1);
    sumC17plusForLCOJF = 0;
    sumH2ForLCOJF = 0;
    
    
    
    % Run for designated length of time
    for curr_month = month_min:month_max
        outputrowidx=mod(curr_month,12);
        % Find month and year
        if mod(curr_month,12) == 0
            year = 2011 + (curr_month/12);
            month = "December"; 
            outputrowidx=12;
        else
           year = 2012 + fix(curr_month/12);
           if mod(curr_month,12) == 1
               month = "January";
               if curr_month~=month_min
                    outputcolumnidx=outputcolumnidx+1;
                end
           elseif mod(curr_month,12) == 2
               month = "February";
           elseif mod(curr_month,12) == 3
               month = "March";
           elseif mod(curr_month,12) == 4
               month = "April";
           elseif mod(curr_month,12) == 5
               month = "May";
           elseif mod(curr_month,12) == 6
               month = "June";
           elseif mod(curr_month,12) == 7
               month = "July";
           elseif mod(curr_month,12) == 8
               month = "August";
           elseif mod(curr_month,12) == 9
               month = "September";
           elseif mod(curr_month,12) == 10
               month = "October";
           elseif mod(curr_month,12) == 11
               month = "November";
           end
        end
    
        % Pass input CO2
        user_CO2 = local_monthly_CO2_cons(curr_month);
        
    
        if elec_temp == true
            fprintf('For %s ', month); fprintf('%d\n', year);
            user_temp = input('Enter an electrolyzer temperature between 100 and 1000C: ');
            temp_min = 100; temp_max = 1000; % deg C
            if user_temp < temp_min || user_temp > temp_max
                disp('Sorry, not a valid input! Please restart.');
                return
            end
        end
    
      
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---CO-ELECTROLYSIS---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % Call co-electrolysis function
        [user_CO2,user_temp,e_prod,mass_b_feed,mass_a_recy,mass_b_recy,mass_c_recy,mass_d_recy,mass_a_out,...
            mass_b_out,mass_c_out,mass_d_out,mass_a_prod,mass_b_prod,mass_c_prod,mass_d_prod,mass_e,...
            mm_a,mm_b,mm_c,mm_d,mm_e,c_prod,d_prod,elec_E,heat_E,syngas_out]...
            = Co_Electrolysis(user_CO2,user_temp);
    
    
        % Totals kg
        sum_a_in = user_CO2 + mass_a_recy;
        sum_b_in = mass_b_feed + mass_b_recy;
        sum_mass_in = sum_a_in + sum_b_in + mass_c_recy + mass_d_recy;
        sum_mass_out = mass_a_out + mass_b_out + mass_c_out + mass_d_out + mass_e;
        sum_mass_prod = mass_a_prod + mass_b_prod + mass_c_prod + mass_d_prod + mass_e;
        tl_E = elec_E + heat_E;
    
        % Display co-electrolysis outputs
        fprintf('\nCO-ELECTROLYSIS\n');
        fprintf('Electrolyser temperature: %f deg C\n', user_temp);
        fprintf('Mass of user-entered CO2: %f kg\n', user_CO2);
        fprintf('Mass of recycled CO2: %f kg\n', mass_a_recy);
        fprintf('Total mass of CO2 electrolyser input: %f kg\n', sum_a_in);
        fprintf('Mass of feedstock H2O: %f kg\n', mass_b_feed);
        fprintf('Mass of recycled H2O: %f kg\n', mass_b_recy);
        fprintf('Total mass of H2O electrolyser input: %f kg\n', sum_b_in);
        fprintf('Mass of recycled CO: %f kg\n', mass_c_recy);
        fprintf('Mass of recycled H2: %f kg\n', mass_d_recy);
        fprintf('Total input mass: %f kg\n\n', sum_mass_in);
    
        fprintf('Mass of output CO2: %f kg\n', mass_a_out);
        fprintf('Mass of output H2O: %f kg\n', mass_b_out);
        fprintf('Mass of output CO: %f kg\n', mass_c_out);
        fprintf('Mass of output H2: %f kg\n', mass_d_out);
        fprintf('Mass of output O2: %f kg\n', mass_e);
        fprintf('Total output mass (post-electrolyser): %f kg\n\n', sum_mass_out);
    
        fprintf('Mass of product CO2: %f kg\n', mass_a_prod);
        fprintf('Mass of product H2O: %f kg\n', mass_b_prod);
        fprintf('Mass of product CO: %f kg\n', mass_c_prod);
        fprintf('Mass of product H2: %f kg\n', mass_d_prod);
        fprintf('Mass of output O2: %f kg\n', mass_e);
        fprintf('Total product mass: %f kg\n\n', sum_mass_prod); 
    
        fprintf('Electrical energy required at user-set temp. of %f deg C', user_temp);
        fprintf(' is %f kWh.\n', elec_E)
        fprintf('Heat energy required at user-set temp. of %f deg C', user_temp);
        fprintf(' is %f kWh.\n', heat_E)
        fprintf('Total energy consumed: %f kWh\n\n', tl_E)
    
        % Electrolyzer Mass and Energy Totals (Inlet & Outlet)
        tl_CO2_input_CoElec(curr_month) = sum_a_in; 
        tl_H2O_input_CoElec(curr_month) = sum_b_in; 
        tl_H2O_required_CoElec(curr_month) = mass_b_feed; 
        tl_elec_input_CoElec(curr_month) = elec_E; 
        tl_heat_input_CoElec(curr_month) = heat_E; 
        tl_CO_output_CoElec(curr_month) = mass_c_prod; 
        tl_H2_output_CoElec(curr_month) = mass_d_prod; 
        tl_O2_output_CoElec(curr_month) = mass_e;
    
        % Loop variables
        loop = 1;
        fprintf('LOOP %d: \n', loop);
        sum_light = 0;
        sum_useful = 0;
        sum_heavy = 0;
        sum_hvy_sell = 0;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---FISCHER-TROPSCH SYNTHESIS---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
        % Call FT synthesis function
        [mol_C2to8,mass_c_FT,mass_d_FT,mass_b_FT,FT_C2to8,FT_C9to16,FT_C17plus,...
        mm_C,mm_H,mass_c_unreact] = FischerTropsch(mass_c_prod,mass_d_prod,mm_b,mm_c,mm_d);
        tl_CO_in_FT(1,curr_month) = mass_c_FT; 
        tl_H2_in_FT(1,curr_month) = mass_d_FT; 
        tl_H2O_out_FT(1,curr_month) = mass_b_FT; 
        tl_unreacted_CO_in_FT(1,curr_month) = mass_c_unreact;
        tl_CO_to_FT (1,curr_month) = mass_c_prod;
    
        
    
        % Heat released from synthesis
        Q = ((-2140*(0.85^2) + 15133*0.85 - 49267)/(16 - 2*0.85)) * 0.184; % kJ/kg of product
    
        fprintf(['Total heat released from FT for recapture(- sign indicates ',... 
        'exothermic reaction): %f kJ/kg\n\n'], Q);
    
        loop = 2;
        while loop > 0
            % Initial condition
            if loop == 2
                ini_FT_light_mass = FT_C2to8; 
            end
            if loop ~= 2
                tl_CO_in_FT(loop-1,curr_month) = mass_c_FT; 
                tl_H2_in_FT(loop-1,curr_month) = mass_d_FT; 
                tl_H2O_out_FT(loop-1,curr_month) = mass_b_FT; 
                tl_unreacted_CO_in_FT(loop-1,curr_month) = mass_c_unreact;
                tl_CO_to_FT(loop-1,curr_month) = new_mass_c_react;
            end
    
            % Break loop when light ends produced are less than 1% of light
            % ends produced during first loop
            if FT_C2to8 < 0.01*ini_FT_light_mass
                fprintf('Loop %d not able to be completed.\n', loop-1);
                fprintf('Light cut mass this loop (from FT): %f kg\n', FT_C2to8);
                fprintf('Original light cut mass (from FT): %f kg\n', ini_FT_light_mass);
                fprintf('FT light cuts are less than 1%% of original mass. System ending.\n\n');
           
    
                fprintf('SYSTEM END STATISTICS\n');
                fprintf('Month and year: %s ', month); fprintf('%d\n', year);
                fprintf('Monthly mass of all produced/used light ends: %f kg\n', sum_light);
                fprintf('Monthly mass of all produced useful products: %f kg\n', sum_useful);
                fprintf('Monthly mass of all produced/used heavy ends: %f kg\n', sum_heavy);
                fprintf('Monthly mass of heavy ends to be sold: %f kg\n\n', sum_hvy_sell);
    
                
                FT_H2 = sum(tl_H2_in_FT,2);
                sumH2ForLCOJF = FT_H2(1,1) - sum(tl_H2_output_CoElec,2);

                %% Storing key outputs to tables
    
              
                Table_sum_useful(outputrowidx,outputcolumnidx)=sum_useful;
                Table_sum_heavy(outputrowidx,outputcolumnidx)=sum_heavy;
                Table_sum_heavy_sell(outputrowidx,outputcolumnidx)=sum_hvy_sell;
                Table_sum_light_sell(outputrowidx,outputcolumnidx)=FT_C2to8;
    
                %Process ressources requirements: H2, H2O and O2 
                Table_H2_requirement(outputrowidx,outputcolumnidx)=sum(tl_H2_in_FT(:,curr_month))/0.75+(sum(tl_H2_HC(:,curr_month)));
                Table_H2O_requirement(outputrowidx,outputcolumnidx)=mass_b_feed+sum(tl_H2O_in_ATR(:,curr_month));
                Table_O2_requirement(outputrowidx,outputcolumnidx)=sum(tl_O2_in_ATR(:,curr_month));
                
                %Process net resource consumption: H2O and O2. 
              
                
                Table_H2O_Produced(outputrowidx,outputcolumnidx)=sum(tl_H2O_out_FT(:,curr_month))+mass_b_prod;
                Table_O2_Produced(outputrowidx,outputcolumnidx)=mass_e;
                Table_H2_Produced(outputrowidx,outputcolumnidx)=tl_H2_output_CoElec(curr_month)+sum(tl_H2_out_ATR(:,curr_month));

                Table_H2O_NetRequirement(outputrowidx,outputcolumnidx)=-Table_H2O_requirement(outputrowidx,outputcolumnidx)+ Table_H2O_Produced(outputrowidx,outputcolumnidx);
                Table_O2_NetRequirement(outputrowidx,outputcolumnidx)=mass_e-Table_O2_requirement(outputrowidx,outputcolumnidx);

                % H2 mass balance representing H2 deficit/excess to be bought or sold  
                H2_NetCons=Table_H2_Produced(outputrowidx,outputcolumnidx)-Table_H2_requirement(outputrowidx,outputcolumnidx); 
                H2_Stock=H2_Stock+H2_NetCons;
          
                Table_H2_NetRequirement(outputrowidx,outputcolumnidx)=H2_Stock;
                if H2_Stock<0
                   H2_Stock=0;
                   refmonth=month;
                   % fprintf('Month and year: %s ', month)
                   % fprintf('H2_Stock: %f ', H2_Stock)
                end 
                
               
                H2_Stock=sum(tl_H2_in_FT(:,curr_month))/0.75*0.25+H2_Stock;
                
                
                Table_H2_Stock(outputrowidx,outputcolumnidx)=H2_Stock;

                if outputrowidx==12
                    H2_Stock=0;
                end



                Table_CO2_SOEC_requirement(outputrowidx,outputcolumnidx)=local_monthly_CO2_cons(curr_month);
                Table_H2O_SOEC_requirement(outputrowidx,outputcolumnidx)=tl_H2O_required_CoElec(curr_month); 
                Table_Heat_SOEC_requirement(outputrowidx,outputcolumnidx)=tl_heat_input_CoElec(curr_month);
                Table_Elec_SOEC_requirement(outputrowidx,outputcolumnidx)=tl_elec_input_CoElec(curr_month);

                Table_HC_sizingFactor(outputrowidx,outputcolumnidx)=sum(tl_C17plus_FT(1,curr_month));
                Table_FT_sizingFactor(outputrowidx,outputcolumnidx)= tl_C2to8_FT(1,curr_month)+ tl_C9to16_FT(1,curr_month)+tl_C17plus_FT(1,curr_month); 
                
                Table_ATR_sizingFactor(outputrowidx,outputcolumnidx)=tl_Syngas_out_ATR(1,curr_month);

                
                break
            end
    
        % Display FT outputs
        fprintf('FT\n');
        fprintf('Mass of light products to be sent to ATR (C2-8): %f kg\n', FT_C2to8);
        fprintf('Mass of useful products (C9-16): %f kg\n', FT_C9to16);
        fprintf('Mass of heavy products to be sent to hydrocracker (C17+): %f kg\n\n', FT_C17plus);
    
        % FT Array Totals
        tl_C2to8_FT(loop-1,curr_month) = FT_C2to8; 
        tl_C9to16_FT(loop-1,curr_month) = FT_C9to16; 
        tl_C17plus_FT(loop-1,curr_month) = FT_C17plus; 
        

        % FT Totals
        sum_light = sum_light + FT_C2to8; 
        sum_useful = sum_useful + FT_C9to16; 
        sum_heavy = sum_heavy + FT_C17plus; 
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---HYDROCRAKER & ATR---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % Hydrocracker H2 requirement
        HC_H2 = 0.031 * FT_C17plus;
        tl_H2_HC(loop-1,curr_month) = HC_H2;
    
        
        [mass_gases,mass_mid_distillates,mass_heavy_cuts,HC_res] = Hydrocracker(FT_C17plus,mm_H);
    
        % Display hydrocracker outputs
        fprintf('HYDROCRACKER: Heavy Ends input\n');
        fprintf('Mass of H2 used for hydrocracking: %f kg\n', HC_H2);
        fprintf('Mass of light products to be sent to ATR (C2-8): %f kg\n', mass_gases);
        fprintf('Mass of useful products (C9-16): %f kg\n', mass_mid_distillates);
        fprintf('Mass of heavy products to be sold (C17+): %f kg\n\n', mass_heavy_cuts);
    
        % Hydrocracker Array Totals
        tl_gases_HC(loop-1,curr_month) = mass_gases;
        tl_naphs_HC(loop-1,curr_month) = mass_mid_distillates;
        tl_hvy_dtls_HC(loop-1,curr_month) = mass_heavy_cuts;
    
        %Hydrocracker balance
        sum_light = sum_light + mass_gases; 
        sum_useful = sum_useful + mass_mid_distillates; 
        sum_heavy = sum_heavy + mass_heavy_cuts; 
        sum_hvy_sell = sum_hvy_sell + mass_heavy_cuts;
    
        fprintf('LOOP %d ', loop-1);
        fprintf('cumulative sum of useful (C9-16) products: %f kg\n\n', sum_useful);
    
        ATR_input = mol_C2to8+(mass_gases/(4.36*mm_C)+(10.72*mm_H));   
        
        fprintf('e_prod: %f\n',e_prod)
      
        [e_ATR_1,b_ATR_1,c_ATR_1,d_ATR_1] = ATR(ATR_input,e_prod);
        mass_e_ATR_1 = e_ATR_1*mm_e;
        mass_b_ATR_1 = b_ATR_1*mm_b;
        mass_c_ATR_1 = c_ATR_1*mm_c;
        mass_d_ATR_1 = d_ATR_1*mm_d;
    
   
        % ATR Array Totals
        tl_lightEnds_in_ATR(loop-1,curr_month) = FT_C2to8 + sum(HC_res);
        tl_O2_in_ATR(loop-1,curr_month) = mass_e_ATR_1; 
        tl_H2O_in_ATR(loop-1,curr_month) = mass_b_ATR_1; 
        tl_CO_out_ATR(loop-1,curr_month) = mass_c_ATR_1; 
        tl_H2_out_ATR(loop-1,curr_month) = mass_d_ATR_1; 
    
        
        c_ATR = c_ATR_1;
        mass_c_ATR = c_ATR*mm_c;
        d_ATR = d_ATR_1; 
        mass_d_ATR = d_ATR * mm_d;
        tl_Syngas_out_ATR(loop-1,curr_month) = (mass_c_ATR+mass_d_ATR);
        % Next Loop
        % Display loop round
        fprintf('LOOP %d: \n', loop);
        
        new_mass_c_react = mass_c_ATR + mass_c_unreact; 
  
    [mol_C2to8,mass_c_FT,mass_d_FT,mass_b_FT,FT_C2to8,FT_C9to16,FT_C17plus,...
        mm_C,mm_H,mass_c_unreact] = FischerTropsch(new_mass_c_react,mass_d_ATR,mm_b,mm_c,mm_d);
    
       
    
        loop = loop + 1;
        end
    
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---saving outputs to designated file---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %inserting headers-years number from 1 to n
    Table_sum_useful=[header;Table_sum_useful(1:end,:)];
    Table_sum_heavy=[header;Table_sum_heavy(1:end,:)];
    Table_sum_heavy_sell=[header;Table_sum_heavy_sell(1:end,:)];
    Table_sum_light_sell=[header;Table_sum_light_sell(1:end,:)];

    Table_H2O_requirement=[header;Table_H2O_requirement(1:end,:)];
    Table_H2_requirement=[header;Table_H2_requirement(1:end,:)];
    Table_O2_requirement=[header;Table_O2_requirement(1:end,:)];

    Table_H2O_Produced=[header;Table_H2O_Produced(1:end,:)];
    Table_O2_Produced=[header;Table_O2_Produced(1:end,:)];
    Table_H2_Produced=[header;Table_H2_Produced(1:end,:)];


    Table_H2O_NetRequirement=[header;Table_H2O_NetRequirement(1:end,:)];
    Table_H2_NetRequirement=[header;Table_H2_NetRequirement(1:end,:)];
    Table_O2_NetRequirement=[header;Table_O2_NetRequirement(1:end,:)];
    Table_H2_Stock=[header;Table_H2_Stock(1:end,:)];

    Table_CO2_SOEC_requirement=[header;Table_CO2_SOEC_requirement(1:end,:)];
    Table_Heat_SOEC_requirement=[header;Table_Heat_SOEC_requirement(1:end,:)];
    Table_Elec_SOEC_requirement=[header;Table_Elec_SOEC_requirement(1:end,:)];

    Table_HC_sizingFactor=[header;Table_HC_sizingFactor(1:end,:)];
    Table_FT_sizingFactor=[header;Table_FT_sizingFactor(1:end,:)];
    Table_ATR_sizingFactor=[header;Table_ATR_sizingFactor(1:end,:)];

    
    folder = fileparts(which(mfilename)); 
    addpath(genpath(folder));

    subfolder=fullfile(pwd,'Output/FedFuelProduction/');
    ProductionResultsFile = strcat('FedResults',num2str(startingYear,'%d'),'_',user_select,'xGrowth','_',date,'.xlsx');

    % % Fuel product fraction results
    writematrix(Table_sum_useful,ProductionResultsFile,"Sheet","Monthly_SumUseful_Kg");
    writematrix(Table_sum_heavy,ProductionResultsFile,"Sheet","Monthly_SumHeavy_Kg");
    writematrix(Table_sum_heavy_sell,ProductionResultsFile,"Sheet","Monthly_SumHeavySell_Kg");
    writematrix(Table_sum_light_sell,ProductionResultsFile,"Sheet","Monthly_Sumlight_Kg");

    % % Fuel production component requirements
    writematrix(Table_H2O_requirement,ProductionResultsFile,"Sheet","Monthly_H2ORequirement_Kg");
    writematrix(Table_H2_requirement,ProductionResultsFile,"Sheet","Monthly_H2Requirement_Kg");
    writematrix(Table_O2_requirement,ProductionResultsFile,"Sheet","Monthly_O2Requirement_Kg");


    %%Produced components
    writematrix(Table_H2O_Produced,ProductionResultsFile,"Sheet","Monthly_H2OProd_FT_Kg");
    writematrix(Table_O2_Produced,ProductionResultsFile,"Sheet","Monthly_O2Prod_FT_Kg");
    writematrix(Table_H2_Produced,ProductionResultsFile,"Sheet","Monthly_H2Prod_FT_Kg");
  
    % % Fuel production required net input, Negative value = deficit
   
    writematrix(Table_H2O_NetRequirement,ProductionResultsFile,"Sheet","Monthly_NetH2ORequirement_Kg");
    writematrix(Table_H2_NetRequirement,ProductionResultsFile,"Sheet","Monthly_NetH2Requirement_Kg");
    writematrix(Table_O2_NetRequirement,ProductionResultsFile,"Sheet","Monthly_NetO2Requirement_Kg");

    writematrix(Table_CO2_SOEC_requirement,ProductionResultsFile,"Sheet","CO2_input_SOEC_Kg");
    writematrix(Table_Heat_SOEC_requirement,ProductionResultsFile,"Sheet","Heat_input_SOEC_kWh");
    writematrix(Table_Elec_SOEC_requirement,ProductionResultsFile,"Sheet","Elec_input_SOEC_kWh");

    % %Unit sizing factors
    writematrix(Table_HC_sizingFactor,ProductionResultsFile,"Sheet","HC_sizingFactor_kg");
    writematrix(Table_FT_sizingFactor,ProductionResultsFile,"Sheet","FTS_sizingFactor_kg");
    writematrix(Table_ATR_sizingFactor,ProductionResultsFile,"Sheet","ATR_sizingFactor_kg");



    movefile(ProductionResultsFile,fullfile(subfolder,ProductionResultsFile));
 
end

