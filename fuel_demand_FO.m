

function [monthly_CO2_cons,monthly_ATF_cons,month_prov_ATF,CO2_mass_yearly,ATF_mass_yearly,month_min,month_max,startingYear,user_select] = fuel_demand_FO(user_select,startingYear,month_min)

% Retrieve input data on provincial fuel allocation and monthly consumption
addpath('Inputs');
file_airports = 'ProvAirportMvts_FO.txt';
fileID_airports = fopen(file_airports,'r');
prov_list = fscanf(fileID_airports,'%f');
file_months = 'MonthlyATFConsumption_FO.txt';
fileID_months = fopen(file_months, 'r');
monthly_list = fscanf(fileID_months,'%f');

% Constants
dens_ATF = 807.5; % kg/m^3
month_max = 468; %


%SCENARIOS [2.5x,3.73x,2.35x,3.1x]
slope_ATF_arr = [3019 5494.14; 2716.88 4226.26]; 
y_int_arr = [519226.399 311294.339; 544583.97 417796.13];

slope_prodVsCO2 = 0.24;
num_provinces = 6;

if user_select == "2.50"
    slope_ATF_cons = slope_ATF_arr(1,1);
    y_int = y_int_arr(1,1);
elseif user_select == "3.73"
    slope_ATF_cons = slope_ATF_arr(1,2);
    y_int = y_int_arr(1,2);
elseif user_select == "2.35"
    slope_ATF_cons = slope_ATF_arr(2,1);
    y_int = y_int_arr(2,1);
elseif user_select == "3.10"
   slope_ATF_cons = slope_ATF_arr(2,2);
    y_int = y_int_arr(2,2);
else
    disp('Sorry, not a valid input! Please restart.');
    return
end
% Pre-allocate size for arrays
trendline_CO2_mass_res = zeros(1,(month_max-month_min));
trendline_ATF_mass_res = zeros(1,(month_max-month_min));
month_prov_ATF = zeros(6,(month_max-month_min+1));

if mod(month_max,12) == 0
    monthly_ATF_cons = zeros(12,month_max/12);
    monthly_CO2_cons = zeros(12,month_max/12);
    CO2_mass_yearly = zeros(1,month_max/12);
    ATF_mass_yearly = zeros(1,month_max/12);
    
else
    monthly_ATF_cons = zeros(12,(month_max + (12-mod(month_max,12)))/12);
    monthly_CO2_cons = zeros(12,(month_max + (12-mod(month_max,12)))/12);
    CO2_mass_yearly = zeros(1,(month_max + (12-mod(month_max,12)))/12);
    ATF_mass_yearly = zeros(1,(month_max + (12-mod(month_max,12)))/12);
    
end

count = 1;
year = 0;
for month = month_min:month_max
    
    % Volume and mass of monthly ATF production based on trendline
    trendline_month_vol_ATF = (slope_ATF_cons * month) + y_int; 
    trendline_month_mass_ATF = trendline_month_vol_ATF * dens_ATF;
    
    % Calculate corresponding CO2 requirement 
    trendline_month_mass_CO2 = trendline_month_mass_ATF / slope_prodVsCO2;
    
    % Store trendline-based results in array
    trendline_CO2_mass_res(month) = trendline_month_mass_CO2;
    trendline_ATF_mass_res(month) = trendline_month_mass_ATF;
    

    if mod(month,12) == 0
       year = year + 1;
% Calculate mass of CO2 produced every December
       ATF_mass_yearly(year)= trendline_ATF_mass_res(month)/monthly_list(12);

       % Calculate total mass of CO2 produced annually        
       CO2_mass_yearly(year) = ATF_mass_yearly(year)/slope_prodVsCO2; 
    end

    
    % Monthly ATF allocation
    if month >= 72 && mod(month,12) == 0
        for temp = (month-11):month
            if mod(temp,12) == 0
                monthly_ATF_cons(temp) = monthly_list(12) * ATF_mass_yearly(year)
                monthly_CO2_cons(temp) = monthly_list(12) * CO2_mass_yearly(year);%check that this is right
            else
                monthly_ATF_cons(temp) = monthly_list(mod(temp,12)) * ATF_mass_yearly(year)
                monthly_CO2_cons(temp) = monthly_list(mod(temp,12)) * CO2_mass_yearly(year);
            end
            
            % Provincial ATF allocation (1-6: BC,ON,QC,AB,MB,NB)
            for prov_ID = (num_provinces*(count-1)+1):num_provinces*count
                
                if count >= 2 && mod(prov_ID,num_provinces) ~= 0
                    month_prov_ATF(prov_ID) = prov_list(mod(prov_ID,6)) * monthly_ATF_cons(temp);
                elseif count >= 2 && mod(prov_ID,num_provinces) == 0
                    month_prov_ATF(prov_ID) = prov_list(6) * monthly_ATF_cons(temp);
                else 
                    month_prov_ATF(prov_ID) = prov_list(prov_ID) * monthly_ATF_cons(temp);
                end
            end
        count = count + 1;
        end
    end
end

% Close textfiles
fclose(fileID_airports);
fclose(fileID_months);
