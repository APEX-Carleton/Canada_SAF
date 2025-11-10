
%% Main inputs: 

% General input and techno-economic parameters and arrays
provinces=["BC" "ON" "QC" "AB" "MB" "NB"];
numProv = 6;
H2color=["Blue","Green"];
Co2ProductionMethod  ="DAC";
CombustionSource="Natural Gas";
ftr_CAPEXdistr=[3.28e9 4.56e9 ;4.17e9 5.79e9;3.16e9 4.39e9;3.73e9 5.18e9];
scenario_Array=["2.50","3.73","2.35","3.10"];
carbonCost=170;% $cad/ton of co2
startingYear=2025;
discountRate = [0.03 0.07]; 
capacityFactor = [0.7 0.8];
hoursPerYear = 8760*capacityFactor; % capacity hr/year
hoursPerMonth = 24*30*capacityFactor; % capacity hr/month
densATF = 0.00305672; % tonne/gal
kmolSyngas=0.087; % converting the ATR syngas from kg to kmol (July01-Fawzi); previousy 0.10985857, why?;
n = 5000;% Number of samples

% Currency and cost indexes
USDtoCad=1.3;
CEPCI_2023=800.8;
CEPCI_2006=499.6;
CEPCI_2009=521.9;
CEPCI_2019=607.5;
    

% Capital and O&M parameters (all dollar values in CAD)
electrolyzerCap = [500 1000 2500]*USDtoCad; % $/kW
electrolyzerOM = [14 45 158]; % $/kWyr
HC_RefCost=1.70e9*USDtoCad;%$
HC_RefCap=50000;%bbl/day
FT_RefCost=[15e6 70e6]*USDtoCad;%$
otherOM = [5 8 18]; % $/kW-yr
averageDens=887;%km/m3, average crude density
ATR_RefCost=16.70e6*USDtoCad;%$
ATR_RefCap=31000;%kmol/hr

% Provincial electricity costs and emissions (BC, ON, QC, AB, MB, NB)
electricityCost = [0.0659 0.0994 0.0526 0.2374 0.0503 0.084]; % $/kWh
electricityEmis = [0.015 0.03 0.0017 0.540 0.002 0.3]; % kgCO2/kWh

% Provincial DAC cost and emissions intensity: BC, ON, QC, AB, MB, NB
DACcost=[175 190;270 320;210 260;230 325;230 325;330 670];
DACElecIntensity=1.93*1000;%kwh/tCO2;

% Hydrogen: costs and emissions by source
blueH2cost = [1.85 1.34];% $/kgH2 
greenH2cost = [3.10 5.01];% $/kgH2 
blueH2emis = [2.3 4.1];% kgCO2/kgH2 
greenH2emis = [0 0.6];%kgCO2/kgH2 

% Heat: costs and emissions for natural gas
combustHeatCost = 36/1000*USDtoCad; % USD per Mwh converted to USD$/kWh
combustHeatEmis = (54431*35.3)/(1000000*37.5); % kgCO2/MJ heat, based off 54431 kgCO2/scf natural gas.0.0512 kgCO2/MJ heat


% Timestamp for input/output file tracking
date=string(datetime('now','Format',"yyyyMMddHHmmss"));
% Timestamp marking when the fuel production results file was generated in main.m.
% This timestamp is appended to the end of the output file name.
% Update this value to match the timestamp of the results files for the scenarios to be analyzed.
timestamp="20240702160026";

% Output folders
outputLcojfFolder =fullfile(pwd,'Output/LCOJF/');
outputCDFFolder=fullfile(pwd,'Output/CDF/');

% Input folder
fedFuelProdFolder =fullfile(pwd,'Output/FedFuelProduction/');


% Add working folder paths
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));


%% Loop through each fuel demand growth scenario
for rootidx=1:4
    GrowthRate=scenario_Array(rootidx);
    ProductionFilename=strcat('FedResults2019_',GrowthRate,'xGrowth','_',timestamp,'.xlsx');
    FileInSubfolder=fullfile(fedFuelProdFolder ,ProductionFilename);

    % Temporarily move the scenario input file to the current working directory for processing
    movefile(FileInSubfolder,fullfile(pwd,ProductionFilename));
    
    % Read process input data required for LCOJF calculations
    tl_CO2_input_CoElec=readmatrix(ProductionFilename,'Sheet',"CO2_input_SOEC_Kg"); % CO2 input for solid oxide electrolyzer cell (SOEC) process [kg]
    sum_useful=readmatrix(ProductionFilename,'Sheet',"Monthly_SumUseful_Kg");% ATF produced
    Table_H2_NetRequirement=readmatrix(ProductionFilename,'Sheet',"Monthly_NetH2Requirement_Kg");% Net hydrogen requirement each month [kg] (negative values indicate external H2 purchase)
    tl_elec_input_CoElec=readmatrix(ProductionFilename,'Sheet',"Elec_input_SOEC_kWh");% Electricity input for SOEC process [kWh]
    tl_heat_input_CoElec = 0;% Assume all heat supply for SOEC is satisfied by Joule heating


    % Unit sizing for core plant sections [kg]
    HC_sizing=readmatrix(ProductionFilename,'Sheet',"HC_sizingFactor_kg");  % Hydrocracker
    FT_sizing=readmatrix(ProductionFilename,'Sheet',"FTS_sizingFactor_kg"); % Fischer-Tropsch Synthesis
    ATR_sizing=readmatrix(ProductionFilename,'Sheet',"ATR_sizingFactor_kg"); % Autothermal Reformer

    % Move input file back to its output folder after processing
    movefile(ProductionFilename,fullfile(fedFuelProdFolder ,ProductionFilename));
    
    
    % Prepare output file names for LCOJF results and distributions
    OFile_LCOJF_MainResults=strcat('LCOJF',num2str(startingYear,'%d'),'_',GrowthRate,'xGrowth','_',date,'.xlsx');
    OFile_LCOJF_CDFPerProvince=strcat('DistribuctionsMC',num2str(startingYear,'%d'),'_',GrowthRate,'xGrowth','_',date,'.xlsx');
    
    % Calculate annual sums of resources (CO2, ATF, heat, electricity, and hydrogen) for the selected starting scenario year
    inputcolumnidx=startingYear-2011;
    yearCO2 = sum(tl_CO2_input_CoElec(2:size(tl_CO2_input_CoElec,1),inputcolumnidx),1); % kg
    yearATF=sum(sum_useful(2:size(sum_useful,1),inputcolumnidx),1);%kg
    yearHeat = 0; %see tl_heat_input_CoElec comment above
    tonneUserCO2 = yearCO2/1000; % convert CO2 from kg into tonnes 
    tonneATF = yearATF/1000;
    MJheat = yearHeat * 3.6; % convert kWh to MJ
    yearElectricity = sum(tl_elec_input_CoElec(2:size(tl_elec_input_CoElec,1),inputcolumnidx),1); % kWh
    yearH2=sum(Table_H2_NetRequirement(2:size(Table_H2_NetRequirement,1),inputcolumnidx),1); % kg

    
    %Select reference Fischer-Tropsch plant size for CAPEX estimation
    FT_RefCap=ftr_CAPEXdistr(rootidx,:);%bbl/day
    
    
   
    
    % Compute the plant lifetime from the starting year
    if startingYear~=2050
        lifetime = 2050-startingYear+1; % years
    else 
        lifetime=20;
    end
    
    % Other Variables
    elecYear = startingYear-2010;
    
    
    %Pre-allocation of arrays

    f_H2Cost = cell(2,1); f_H2Emis = cell(2,1);
   
    % Scenario Distributions
    % Each letter in the sequence below represents a variable input for the Monte Carlo simulation:
    % d: discountRate, i: electrolyzerCap, j: BOPcap, k: DACcost (BC), l: DACcost (ON), m: DACcost (QC),
    % q: DACcost (AB), o: DACcost (MB), p: DACcost (NB), u: hoursPerMonth, v: FT_RefCapex,
    % x: blueH2cost, y: blueH2emis, z: greenH2cost, a: greenH2emis
 
    %d,i,k,l,m,q,o,p,u,v,x,y,z,a
    distsScn = {
        {'Uniform',discountRate(1),discountRate(2)};%d     
        {'Triangular',electrolyzerCap(1),electrolyzerCap(2),electrolyzerCap(3)};%i
        {'Uniform',DACcost(1,1),DACcost(1,2)};%Prov1: BC %k
        {'Uniform',DACcost(2,1),DACcost(2,2)};%Prov2: ON %l
        {'Uniform',DACcost(3,1),DACcost(3,2)};%Prov3: QC %m
        {'Uniform',DACcost(4,1),DACcost(4,2)};%Prov1: AB %q
        {'Uniform',DACcost(5,1),DACcost(5,2)};%Prov1: MC %o
        {'Uniform',DACcost(6,1),DACcost(6,2)};%Prov6: NB %p
        {'Uniform',hoursPerMonth(1),hoursPerMonth(2)}; %u
        {'Uniform',FT_RefCap(1),FT_RefCap(2)}; %v


        {'Uniform',blueH2cost(1),blueH2cost(2)}; %X
        {'Uniform',blueH2emis(1),blueH2emis(2)}; %Y
        {'Uniform',greenH2cost(1),greenH2cost(2)}; %z
        {'Uniform',greenH2emis(1),greenH2emis(2)}; %a
        };
  
    disp('Running');
    
%% Capital Cost Estimation at Year 0
    % Sized to meet the load at the end of the stack lifetime:
    capitalidx=inputcolumnidx+lifetime;
    if capitalidx>=size(tl_CO2_input_CoElec,2)
        capitalidx=size(tl_CO2_input_CoElec,2);
    end
    % Extract maximum monthly load values for sizing calculations
    monthElec=max(tl_elec_input_CoElec(:,capitalidx));
    HCsizingFactor=(max(HC_sizing(:,capitalidx))/averageDens)*6.3/30;%Input must be in kg, divide by 6.2898107704 to convert in barrel, average over 30 days to get barrel per day. Assuming that the heavy end has a density of 887 kg/m3
    FTSsizingFactor=(max(FT_sizing(:,capitalidx))/averageDens)*6.3/30;%Input must be in kg, divide by 6.2898107704 to convert in barrel, average over 30 days to get barrel per day. Assuming that the heavy end has a density of 887 kg/m3
    ATRsizingFactor=max(ATR_sizing(:,capitalidx))*kmolSyngas/(30*24);%must be in kg of syngas produced by the ATR. Convert to kmol/hr 

    % Row indices for writing results in final output sheet
    r_s2 = 4; r_s3 = 4;

    %% Monte Carlo: Electrolyzer Capital Cost
    f_elecCap = @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)((monthElec)./u).*i;
    [~,stdev_elecCap,CI_elecCap,MCfuncVals_elecCap,~,largestBin_lb,largestBin_ub] = simulateMC(f_elecCap,distsScn,n,'hist');clear simulateMC;
    ModeElecCap = (largestBin_lb+largestBin_ub)/2; 
    MinElecCap = min(MCfuncVals_elecCap,[],'all'); 
    MaxElecCap = max(MCfuncVals_elecCap,[],'all');
    writematrix(MinElecCap,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('B',int2str(r_s2)));
    writematrix(ModeElecCap,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('C',int2str(r_s2)));
    writematrix(MaxElecCap,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('D',int2str(r_s2)));
    
    %% Monte Carlo: Hydrocracker Capital Cost
    f_HCcap = @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)(CEPCI_2023/CEPCI_2009).*(HC_RefCost.*((HCsizingFactor/HC_RefCap).^(0.6))); 
    [~,stdev_HCcap,CI_HCcap,MCfuncVals_HCcap,~,largestBin_lb,largestBin_ub] = simulateMC(f_HCcap,distsScn,n,'hist');clear simulateMC;
    ModeHCcap = (largestBin_lb+largestBin_ub)/2; 
    MinHCCap = min(MCfuncVals_HCcap,[],'all'); 
    MaxHCCap = max(MCfuncVals_HCcap,[],'all');
    writematrix(MinHCCap,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('BK',int2str(r_s2)));
    writematrix(ModeHCcap,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('BL',int2str(r_s2)));
    writematrix(MaxHCCap,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('BM',int2str(r_s2)));


    %% Monte Carlo: Fischer-Tropsch Reactor Capital Cost
    f_FTRcap = @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)(v);
    [~,stdev_FTRcap,CI_FTRcap,MCfuncVals_FTRcap,~,largestBin_lb,largestBin_ub] = simulateMC(f_FTRcap,distsScn,n,'hist');clear simulateMC;
    ModeFTRcap = (largestBin_lb+largestBin_ub)/2; 
    MinFTRcap = min(MCfuncVals_FTRcap,[],'all'); 
    MaxFTRcap = max(MCfuncVals_FTRcap,[],'all');
    writematrix(MinFTRcap,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('BN',int2str(r_s2)));
    writematrix(ModeFTRcap,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('BO',int2str(r_s2)));
    writematrix(MaxFTRcap,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('BP',int2str(r_s2)));

    %% Monte Carlo: ATR Capital Cost
    f_ATRcap = @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)(CEPCI_2023/CEPCI_2019).*(ATR_RefCost*((ATRsizingFactor/ATR_RefCap).^(0.6))); 
    [~,stdev_ATRcap,CI_ATRcap,MCfuncVals_ATRcap,~,largestBin_lb,largestBin_ub] = simulateMC(f_ATRcap,distsScn,n,'hist');clear simulateMC;
    ModeATRcap = (largestBin_lb+largestBin_ub)/2; 
    MinATRCap = min(MCfuncVals_ATRcap,[],'all'); 
    MaxATRCap = max(MCfuncVals_ATRcap,[],'all');
    writematrix(MinATRCap,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('BQ',int2str(r_s2)));
    writematrix(ModeATRcap,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('BR',int2str(r_s2)));
    writematrix(MaxATRCap,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('BS',int2str(r_s2)));
    
    %% TOTAL CAPITAL INCLUDING INDIRECT CAPITAL COSTS
    f_totalCapex= @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)(1.4*(4.7*(f_elecCap(d,i,k,l,m,q,o,p,u,v,x,y,z,a) +f_HCcap(d,i,k,l,m,q,o,p,u,v,x,y,z,a)+f_FTRcap(d,i,k,l,m,q,o,p,u,v,x,y,z,a)+f_ATRcap(d,i,k,l,m,q,o,p,u,v,x,y,z,a)))); 
    [~,stdev_totalcap,CI_totalcap,MCfuncVals_totalcap,~,largestBin_lb,largestBin_ub] = simulateMC(f_totalCapex,distsScn,n,'hist');clear simulateMC;
    Modetotalcap = (largestBin_lb+largestBin_ub)/2; 
    MintotalRCap = min(MCfuncVals_totalcap,[],'all'); 
    MaxtotalCap = max(MCfuncVals_totalcap,[],'all');
    writematrix(MintotalRCap,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('BT',int2str(r_s2)));
    writematrix(Modetotalcap,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('BU',int2str(r_s2)));
    writematrix(MaxtotalCap,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('BV',int2str(r_s2)));

    disp("Done with capital costs");

%% Variable Costs Calculation

    % Loop through years from Year 0 onwards. 
    
    inputcolumnidx=startingYear-2011;

    % End-of-life (EOL) values: used to normalize O&M cost allocation for units operating below peak capacity
    EOL_yearElectricity=sum(tl_elec_input_CoElec(2:size(tl_elec_input_CoElec,1),capitalidx),1);
    EOL_HCsizingFactor=sum(HC_sizing(2:size(HC_sizing,1),capitalidx),1); %must be in kg
    EOL_FTRsizingFactor=sum(FT_sizing(2:size(FT_sizing,1),capitalidx),1); %must be in kg
    EOL_ATRsizingFactor=sum(ATR_sizing(2:size(ATR_sizing,1),capitalidx),1); %must be in kg of syngas produced by the ATR
    
    tl_ATF = @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)0; % Placeholder for annual ATF production

    for year = 1:lifetime
        % Annual resource requirements (per year, per scenario)
        tonneUserCO2=sum(tl_CO2_input_CoElec(2:size(tl_CO2_input_CoElec,1),inputcolumnidx),1)/1000;
        yearElectricity= sum(tl_elec_input_CoElec(2:size(tl_elec_input_CoElec,1),inputcolumnidx),1);
        %yearHeat = sum(tl_heat_input_CoElec(2:size(tl_heat_input_CoElec,1),inputcolumnidx),1)/0.8; % kWh. Assume 80% efficiency of the boiler
        %yearHeat = 0;
        %MJheat = yearHeat * 3.6; % convert kWh to MJ
        yearHCProd=sum(HC_sizing(2:size(HC_sizing,1),inputcolumnidx),1)/0.8;
        yearFTRProd=sum(FT_sizing(2:size(FT_sizing,1),inputcolumnidx),1)/0.8;
        yearATRProd=sum(ATR_sizing(2:size(ATR_sizing,1),inputcolumnidx),1)/0.8;
        H2purchased=sum(Table_H2_NetRequirement(2:size(Table_H2_NetRequirement,1),inputcolumnidx),1);

        % If excess H2, set purchased amount to 0 (no external buy)
        if H2purchased<0
            H2purchased=-H2purchased;
        else
            H2purchased=0;
        end

        %% Monte Carlo: O&M Costs for Electrolyzer (scaled by normalized electricity use)
        f_elecOM = @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)(f_elecCap(d,i,k,l,m,q,o,p,u,v,x,y,z,a).*(yearElectricity/EOL_yearElectricity).*0.05)./((1+d).^(year-1));   
        [~,stdev_elecOM,CI_elecOM,MCfuncVals_elecOM,~,largestBin_lb,largestBin_ub] = simulateMC(f_elecOM,distsScn,n,'hist');clear simulateMC;
        ModeElecOM = (largestBin_lb+largestBin_ub)/2; 
        MinElecOM = min(MCfuncVals_elecOM,[],'all'); 
        MaxElecOM = max(MCfuncVals_elecOM,[],'all');
        writematrix(MinElecOM,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('H',int2str(r_s2)));
        writematrix(ModeElecOM,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('I',int2str(r_s2)));
        writematrix(MaxElecOM,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('J',int2str(r_s2)));
        
        %% Monte Carlo: O&M Costs for FTR, HC, ATR (scaled by EOL production values)
        f_otherOM = @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)(f_HCcap(d,i,k,l,m,q,o,p,u,v,x,y,z,a).*(yearHCProd/EOL_HCsizingFactor).*0.03...
            +f_FTRcap(d,i,k,l,m,q,o,p,u,v,x,y,z,a).*(yearFTRProd/EOL_FTRsizingFactor).*0.03+f_ATRcap(d,i,k,l,m,q,o,p,u,v,x,y,z,a).*(yearATRProd/EOL_ATRsizingFactor).*0.02)./((1+d).^(year-1));
        [~,stdev_otherOM,CI_otherOM,MCfuncVals_otherOM,~,largestBin_lb,largestBin_ub] = simulateMC(f_otherOM,distsScn,n,'hist');clear simulateMC;
        ModeOtherOM = (largestBin_lb+largestBin_ub)/2; 
        MinOtherOM = min(MCfuncVals_otherOM,[],'all'); 
        MaxOtherOM = max(MCfuncVals_otherOM,[],'all');
        writematrix(MinOtherOM,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('K',int2str(r_s2)));
        writematrix(ModeOtherOM,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('L',int2str(r_s2)));
        writematrix(MaxOtherOM,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('M',int2str(r_s2)));
        
        %% Monte Carlo: Hydrogen Costs and Emissions (Blue and Green)
        % 1=Blue H2
        f_H2Cost{1} = @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)(H2purchased.*x)./((1+d).^(year-1));
        [~,stdev_BlueH2Cost,CI_BlueH2Cost,MCfuncVals_BlueH2Cost,~,largestBin_lb,largestBin_ub] = simulateMC(f_H2Cost{1},distsScn,n,'hist');clear simulateMC;
        ModeBlueH2Cost = (largestBin_lb+largestBin_ub)/2; 
        MinBlueH2Cost = min(MCfuncVals_BlueH2Cost,[],'all'); 
        MaxBlueH2Cost = max(MCfuncVals_BlueH2Cost,[],'all');
        writematrix(MinBlueH2Cost,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('AX',int2str(r_s2)));
        writematrix(ModeBlueH2Cost,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('AY',int2str(r_s2)));
        writematrix(MaxBlueH2Cost,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('AZ',int2str(r_s2)));

        f_H2Emis{1} = @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)(H2purchased.*y);
        [~,stdev_BlueH2Emis,CI_BlueH2Emis,MCfuncVals_BlueH2Emis,~,largestBin_lb,largestBin_ub] = simulateMC(f_H2Emis{1},distsScn,n,'hist');clear simulateMC;
        ModeBlueH2Emis = (largestBin_lb+largestBin_ub)/2; MinBlueH2Emis = min(MCfuncVals_BlueH2Emis,[],'all'); MaxBlueH2Emis = max(MCfuncVals_BlueH2Emis,[],'all');
        writematrix(MinBlueH2Emis,OFile_LCOJF_MainResults,'Sheet','Emissions','Range',strcat('X',int2str(r_s2)));
        writematrix(ModeBlueH2Emis,OFile_LCOJF_MainResults,'Sheet','Emissions','Range',strcat('Y',int2str(r_s2)));
        writematrix(MaxBlueH2Emis,OFile_LCOJF_MainResults,'Sheet','Emissions','Range',strcat('Z',int2str(r_s2)));

        % 2=Green H2
        f_H2Cost{2} = @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)(H2purchased.*z)./((1+d).^(year-1));
        [~,stdev_GreenH2Cost,CI_GreenH2Cost,MCfuncVals_GreenH2Cost,~,largestBin_lb,largestBin_ub] = simulateMC(f_H2Cost{2},distsScn,n,'hist');clear simulateMC;
        ModeGreenH2Cost = (largestBin_lb+largestBin_ub)/2; 
        MinGreenH2Cost = min(MCfuncVals_GreenH2Cost,[],'all'); 
        MaxGreenH2Cost = max(MCfuncVals_GreenH2Cost,[],'all');
        writematrix(MinGreenH2Cost,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('BA',int2str(r_s2)));
        writematrix(ModeGreenH2Cost,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('BB',int2str(r_s2)));
        writematrix(MaxGreenH2Cost,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('BC',int2str(r_s2)));

        f_H2Emis{2} = @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)(H2purchased.*a);
        [testH2,stdev_GreenH2Emis,CI_GreenH2Emis,MCfuncVals_GreenH2Emis,MCH2_Samples,largestBin_lb,largestBin_ub] = simulateMC(f_H2Emis{2},distsScn,n,'hist');clear simulateMC;
        ModeGreenH2Emis = (largestBin_lb+largestBin_ub)/2; MinGreenH2Emis = min(MCfuncVals_GreenH2Emis,[],'all'); MaxGreenH2Emis = max(MCfuncVals_GreenH2Emis,[],'all');
        writematrix(MinGreenH2Emis,OFile_LCOJF_MainResults,'Sheet','Emissions','Range',strcat('AA',int2str(r_s2)));
        writematrix(ModeGreenH2Emis,OFile_LCOJF_MainResults,'Sheet','Emissions','Range',strcat('AB',int2str(r_s2)));
        writematrix(MaxGreenH2Emis,OFile_LCOJF_MainResults,'Sheet','Emissions','Range',strcat('AC',int2str(r_s2)));


        %% Heat Costs and Emissions
        f_combustCost = @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)(combustHeatCost*yearHeat)./((1+d).^(year-1));
        [funcVal_combustHeatCost,stdev_combustHeatCost,CI_combustHeatCost,MCfuncVals_combustHeatCost] = simulateMC(f_combustCost,distsScn,n,'median');clear simulateMC;  
        MinCombustCost = min(MCfuncVals_combustHeatCost,[],'all'); 
        MaxCombustCost = max(MCfuncVals_combustHeatCost,[],'all');
        writematrix(MinCombustCost,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('BD',int2str(r_s2)));
        writematrix(funcVal_combustHeatCost,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('BE',int2str(r_s2)));
        writematrix(MaxCombustCost,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('BF',int2str(r_s2)));
    
        f_combustEmis = @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)(combustHeatEmis*MJheat);
        [funcVal_combustHeatEmis,~,~,~] = simulateMC(f_combustEmis,distsScn,n,'median');clear simulateMC;  
        writematrix(funcVal_combustHeatEmis,OFile_LCOJF_MainResults,'Sheet','Emissions','Range',strcat('P',int2str(r_s3)));
        writecell({year},OFile_LCOJF_MainResults,'Sheet','Emissions','Range',strcat('A',int2str(r_s3)));

        

        %% Electricity and DAC CO2 Costs and Emissions by Province
        emissionsCosts_col_idx=30;
        for elecIdx = 1:numProv

            % Monte Carlo: Electricity cost per province for the year's total electricity use
            f_ElecCost{elecIdx} = @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)(electricityCost(elecIdx)*(yearElectricity))./((1+d).^(year-1));
            [funcVal_elecCost,stdev_elecCost,CI_elecCost,MCfuncVals_elecCost] = simulateMC(f_ElecCost{elecIdx},distsScn,n,'median');clear simulateMC;
            MinElecCost = min(MCfuncVals_elecCost,[],'all'); 
            MaxElecCost = max(MCfuncVals_elecCost,[],'all');

            costCols = {'N','Q','T','W','Z','AC'};  % Minimum
            meanCols = {'O','R','U','X','AA','AD'}; % Median
            maxCols  = {'P','S','V','Y','AB','AE'}; % Maximum
            writematrix(MinElecCost, OFile_LCOJF_MainResults, 'Sheet', 'Costs', 'Range', strcat(costCols{elecIdx}, int2str(r_s2)));
            writematrix(funcVal_elecCost, OFile_LCOJF_MainResults, 'Sheet', 'Costs', 'Range', strcat(meanCols{elecIdx}, int2str(r_s2)));
            writematrix(MaxElecCost, OFile_LCOJF_MainResults, 'Sheet', 'Costs', 'Range', strcat(maxCols{elecIdx}, int2str(r_s2)));

            % Monte Carlo: Electricity emissions per province for the year's total electricity use
            f_ElecEmis{elecIdx} = @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)(electricityEmis(elecIdx)*(yearElectricity));
            [funcVal_elecEmis,stdev_elecEmis,CI_elecEmis,~] = simulateMC(f_ElecEmis{elecIdx},distsScn,n,'median');clear simulateMC;
            emisCols = {'B','C','D','E','F','G'};
            dacVars = {'k','l','m','q','o','p'};
            writematrix(funcVal_elecEmis, OFile_LCOJF_MainResults, 'Sheet', 'Emissions', 'Range', strcat(emisCols{elecIdx}, int2str(r_s3)));
            f_CO2cost{elecIdx} = @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)(dacVars{elecIdx}.*tonneUserCO2)./((1+d).^(year-1));
            [funcVal_CCScost,stdev_CCScost,CI_CCScost,MCfuncVals_CCScost] = simulateMC(f_CO2cost{elecIdx},distsScn,n,'median');clear simulateMC;
            MinCCSCost = min(MCfuncVals_CCScost,[],'all'); 
            MaxCCSCost = max(MCfuncVals_CCScost,[],'all');
            ccsCostCols  = {'AF','AI','AL','AO','AR','AU'};
            ccsMeanCols  = {'AG','AJ','AM','AP','AS','AV'};
            ccsMaxCols   = {'AH','AK','AN','AQ','AT','AW'};
            
            writematrix(MinCCSCost, OFile_LCOJF_MainResults, 'Sheet', 'Costs', 'Range', strcat(ccsCostCols{elecIdx}, int2str(r_s2)));
            writematrix(funcVal_CCScost, OFile_LCOJF_MainResults, 'Sheet', 'Costs', 'Range', strcat(ccsMeanCols{elecIdx}, int2str(r_s2)));
            writematrix(MaxCCSCost, OFile_LCOJF_MainResults, 'Sheet', 'Costs', 'Range', strcat(ccsMaxCols{elecIdx}, int2str(r_s2)));
            
            % Monte Carlo: DAC CO2 emissions per province (from DAC electricity consumption)
            dacEmisCols = {'H','I','J','K','L','M'};
            f_DACemis{elecIdx} = @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)(electricityEmis(elecIdx)*(tonneUserCO2*DACElecIntensity));
            [funcVal_CCSemis,~,~,~] = simulateMC(f_DACemis{1},distsScn,n,'median');clear simulateMC;
            writematrix(funcVal_CCSemis,OFile_LCOJF_MainResults,'Sheet','Emissions','Range',strcat(dacEmisCols{elecIdx},int2str(r_s3)));

            % Monte Carlo: DAC electricity demand per province
            dacElecCols = {'R','S','T','U','V','W'};
            f_DACElec{elecIdx} = @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)tonneUserCO2*DACElecIntensity;
            [funcVal_DACElec,~,~,~] = simulateMC(f_DACElec{elecIdx},distsScn,n,'median');clear simulateMC;
            writematrix(funcVal_DACElec,OFile_LCOJF_MainResults,'Sheet','Emissions','Range',strcat(dacElecCols{elecIdx},int2str(r_s3)));

            % !!Reporting only the cost of green hydrogen emissions
            for H2Idx = [1 2]
                f_CO2emisCost{elecIdx,H2Idx} = @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)(((f_ElecEmis{elecIdx}(d,i,k,l,m,q,o,p,u,v,x,y,z,a)+f_combustEmis(d,i,k,l,m,q,o,p,u,v,x,y,z,a))./1000).*carbonCost)./((1+d).^(year-1));
                        [funcVal_CO2emisCost,~,CI_CO2emisCost,MCfuncVals_CO2emisCost,~,largestBin_lb,largestBin_ub] = simulateMC(f_CO2emisCost{elecIdx,H2Idx},distsScn,n,'hist');clear simulateMC;
                
                ModeCO2emisCost = (largestBin_lb+largestBin_ub)/2; 
                MinCO2emisCost = min(MCfuncVals_CO2emisCost,[],'all'); 
                MaxCO2emisCost = max(MCfuncVals_CO2emisCost,[],'all');
                % writematrix(MinCO2emisCost,OFile_LCOJF_MainResults,'Sheet','Emissions','Range',strcat('AD',int2str(r_s2)));
                % writematrix(ModeCO2emisCost,OFile_LCOJF_MainResults,'Sheet','Emissions','Range',strcat('AE',int2str(r_s2)));
                % writematrix(MaxCO2emisCost,OFile_LCOJF_MainResults,'Sheet','Emissions','Range',strcat('AF',int2str(r_s2)));
                writecell({ "provinces(elecIdx)" "H2color(H2Idx)"},OFile_LCOJF_MainResults,'Sheet','Emissions','Range',strcat(num2xlcol(emissionsCosts_col_idx),int2str(2)));
                writecell({ "Min" "Peak (Mode)" "Max"},OFile_LCOJF_MainResults,'Sheet','Emissions','Range',strcat(num2xlcol(emissionsCosts_col_idx),int2str(3)))
                writematrix(MinCO2emisCost,OFile_LCOJF_MainResults,'Sheet','Emissions','Range',strcat(num2xlcol(emissionsCosts_col_idx),int2str(r_s2)));
                writematrix(ModeCO2emisCost,OFile_LCOJF_MainResults,'Sheet','Emissions','Range',strcat(num2xlcol(emissionsCosts_col_idx+1),int2str(r_s2)));
                writematrix(MaxCO2emisCost,OFile_LCOJF_MainResults,'Sheet','Emissions','Range',strcat(num2xlcol(emissionsCosts_col_idx+2),int2str(r_s2)));
                emissionsCosts_col_idx=emissionsCosts_col_idx+3;

                
            end
      
        end    
        

   
        %% Target ATF production increases while still being discounted to reflect a lesser value
        tonneATF=sum(sum_useful(2:size(sum_useful,1),inputcolumnidx),1)/1000;
        f_ATF = @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)(tonneATF)./((1+d).^(year-1)); % tonnes
        [funcVal_ATF,stdev_ATF,CI_ATF,MCfuncVals_ATF] = simulateMC(f_ATF,distsScn,n,'median');clear simulateMC;
        MinATF = min(MCfuncVals_ATF,[],'all'); MaxATF = max(MCfuncVals_ATF,[],'all');
        writematrix(MinATF,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('BH',int2str(r_s2)));
        writematrix(funcVal_ATF,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('BI',int2str(r_s2)));
        writematrix(MaxATF,OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('BJ',int2str(r_s2)));
        writecell({year},OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('A',int2str(r_s2)));
        % Lifetime ATF totals
        tl_ATF = @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)f_ATF(d,i,k,l,m,q,o,p,u,v,x,y,z,a) + tl_ATF(d,i,k,l,m,q,o,p,u,v,x,y,z,a);
        
        
        
        r_s2 = r_s2 + 1; r_s3 = r_s3 + 1;
    
    
        % Total costs per scenario each year
        for elecIdx = [1 2 3 4 5 6]
            % for CO2Idx = [1 2] 
                for H2Idx = [1 2]              
                   
                    f_province{year,elecIdx,H2Idx} = @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)f_elecOM(d,i,k,l,m,q,o,p,u,v,x,y,z,a) + f_otherOM(d,i,k,l,m,q,o,p,u,v,x,y,z,a)...
                        + f_ElecCost{elecIdx}(d,i,k,l,m,q,o,p,u,v,x,y,z,a) + f_CO2cost{elecIdx}(d,i,k,l,m,q,o,p,u,v,x,y,z,a)...
                        + f_combustCost(d,i,k,l,m,q,o,p,u,v,x,y,z,a)+f_CO2emisCost{elecIdx,H2Idx}(d,i,k,l,m,q,o,p,u,v,x,y,z,a)...
                    +f_H2Cost{H2Idx}(d,i,k,l,m,q,o,p,u,v,x,y,z,a);
                    if year == 1
                       f_province{year,elecIdx,H2Idx} = @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)f_province{year,elecIdx,H2Idx}(d,i,k,l,m,q,o,p,u,v,x,y,z,a) + f_totalCapex(d,i,k,l,m,q,o,p,u,v,x,y,z,a);
                    end
                    [~,stdev_prov1235,CI_prov1235,MCfuncVals_prov1235,~,largestBin_lb,largestBin_ub] = simulateMC(f_province{year,elecIdx,H2Idx},distsScn,n,'hist');clear simulateMC;
                    ModeScn1235(year,elecIdx,H2Idx) = (largestBin_lb+largestBin_ub)/2;
                    MinScn1235(year,elecIdx,H2Idx) = min(MCfuncVals_prov1235,[],'all');
                    MaxScn1235(year,elecIdx) = max(MCfuncVals_prov1235,[],'all');                    
                end
        end
       
        elecYear = elecYear + 1;
        
    
        if inputcolumnidx>=size(tl_CO2_input_CoElec,2)
            inputcolumnidx=size(tl_CO2_input_CoElec,2);
        else
            inputcolumnidx=inputcolumnidx+1;
        end
    end
    [~,stdev_tlATF,CI_tlATF,MCfuncVals_tlATF,~]=simulateMC(tl_ATF,distsScn,n,'median');clear simulateMC;
    disp("Done with OPEX") 
   
    
    %% Lifetime Total Costs and LCOJF Calculation by Province & Scenario
    r_s1 = 4; % Starting row for results
    for elecIdx = [1 2 3 4 5 6]
          for H2Idx = [1 2]
                tl_1235costs = @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)0; % Initialize a cumulative cost function for the scenario
                for year = 1:lifetime % Sum annual costs across scenario years
                    tl_1235costs = @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)f_province{year,elecIdx,H2Idx}(d,i,k,l,m,q,o,p,u,v,x,y,z,a) + tl_1235costs(d,i,k,l,m,q,o,p,u,v,x,y,z,a);
                  
                end               
                [~,stdev_tl1235costs,CI_tl1235costs,MCfuncVals_tl1235costs,~,largestBin_lb,largestBin_ub] = simulateMC(tl_1235costs,distsScn,n,'hist');clear simulateMC;
          
                % Monte Carlo: Levelized Cost of Jet Fuel ($/GAL) for this province and H2 case
                LCOJF1235dPg = @(d,i,k,l,m,q,o,p,u,v,x,y,z,a)densATF*((tl_1235costs(d,i,k,l,m,q,o,p,u,v,x,y,z,a)./tl_ATF(d,i,k,l,m,q,o,p,u,v,x,y,z,a)));
                [~,stdev_LCOJF1235,CI_LCOJF1235,MCfuncVals_LCOJF1235,MCsamples_LCOJF,largestBin_lb,largestBin_ub] = simulateMC(LCOJF1235dPg,distsScn,n,'hist');clear simulateMC;

                ModeLCOJF(elecIdx,H2Idx) = (largestBin_lb+largestBin_ub)/2;
                MinLCOJF(elecIdx,H2Idx) = min(MCfuncVals_LCOJF1235,[],'all');
                MaxLCOJF(elecIdx,H2Idx) = max(MCfuncVals_LCOJF1235,[],'all');

               % Output results to the LCOJF sheet
                Header={provinces(elecIdx),Co2ProductionMethod  ,H2color(H2Idx),CombustionSource};
                writecell(Header,OFile_LCOJF_MainResults,'Sheet','LCOJF','Range',strcat('A',int2str(r_s1)));
                writematrix(stdev_LCOJF1235,OFile_LCOJF_MainResults,'Sheet','LCOJF','Range',strcat('E',int2str(r_s1)));
                writematrix(CI_LCOJF1235(1),OFile_LCOJF_MainResults,'Sheet','LCOJF','Range',strcat('F',int2str(r_s1)));
                writematrix(CI_LCOJF1235(2),OFile_LCOJF_MainResults,'Sheet','LCOJF','Range',strcat('G',int2str(r_s1)));               
                writematrix(MinLCOJF(elecIdx,H2Idx),OFile_LCOJF_MainResults,'Sheet','LCOJF','Range',strcat('H',int2str(r_s1)));
                writematrix(ModeLCOJF(elecIdx,H2Idx),OFile_LCOJF_MainResults,'Sheet','LCOJF','Range',strcat('I',int2str(r_s1)));
                writematrix(MaxLCOJF(elecIdx,H2Idx),OFile_LCOJF_MainResults,'Sheet','LCOJF','Range',strcat('J',int2str(r_s1)));

                % Output CDF MC results for each province to separate sheet for distribution analysis
                writecell({"LCOJF ($/GAL)"},OFile_LCOJF_CDFPerProvince,'Sheet',provinces(elecIdx),'Range',strcat('A',int2str(1)));
                writematrix(MCfuncVals_LCOJF1235,OFile_LCOJF_CDFPerProvince,'Sheet',provinces(elecIdx),'Range',strcat('A',int2str(2)));
                r_s1 = r_s1 + 1;
               
          end
      
    end 

    
 

   %% Add headings to all generated sheets in the LCOJF output files 
    % Add headings to the "LCOJF" sheet
    writecell({"Lifetime LCOJF Scenarios (including capital and O&M costs)" "" "" "" "Std Dev" "Confidence Interval" "" "Triangular" ""},OFile_LCOJF_MainResults,'Sheet','LCOJF','Range',strcat('A',int2str(2)));
    writecell({"Province" "CO2 Production" "H2 Production" "Heat Production" "" "Lowerbound" "Upperbound" "Min" "Peak (mode of frequency)" "Max"},OFile_LCOJF_MainResults,'Sheet','LCOJF','Range',strcat('A',int2str(3)));
    
    
    % Add  headings to the "Costs" sheet
    writecell({"Main Capital"},OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('B',int2str(1)));
    writecell({"Electricity Cost"},OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('N',int2str(1)));
    writecell({"CO2 Cost DAC"},OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('AF',int2str(1)));
    writecell({"Hydrogen Cost" "" "" "" "" "" "Heat Production Cost" "" "" "" "ATF(tonnes)" "" "" "Other Capex Components" "" "" "" "" "" " " "" "" "Total Capex"},OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('AX',int2str(1)));

    writecell({"Electrolyzer Capital" "" "" "BOP Capital" "" "" "Electrolyzer O&M" "" "" "Other O&M"},OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('B',int2str(2)));
    writecell({"BC" "" "" "ON" "" "" "QC" "" "" "AB" "" "" "MB" "" "" "NB"},OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('N',int2str(2)));
    writecell({"BC" "" "" "ON" "" "" "QC" "" "" "AB" "" "" "MB" "" "" "NB"},OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('AF',int2str(2)));
    writecell({"Blue" "" "" "Green" "" "" "Combustion" "" "" "" "ATF(tonnes)" "" "" "Hydrocracker Capex" "" "" "FTS Capex" "" "" "ATR Capex" "" "" "Total Capex" "" ""},OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('AX',int2str(2)));

    writecell({"Year" "Min" "Peak (Mode)" "Max" "" "" "" "Min" "Peak (Mode)" "Max" "Min" "Peak (Mode)" "Max"},OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('A',int2str(3)));
    writecell({"Min" "Peak (Mode)" "Max" "Min" "Peak (Mode)" "Max" "Min" "Peak (Mode)" "Max" "Min" "Peak (Mode)" "Max" "Min" "Peak (Mode)" "Max" "Min" "Peak (Mode)" "Max"},OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('N',int2str(3)));
    writecell({"Min" "Peak (Mode)" "Max" "Min" "Peak (Mode)" "Max" "Min" "Peak (Mode)" "Max" "Min" "Peak (Mode)" "Max" "Min" "Peak (Mode)" "Max" "Min" "Peak (Mode)" "Max"},OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('AF',int2str(3)));
    writecell({"Min" "Peak (Mode)" "Max" "Min" "Peak (Mode)" "Max" "Min" "Peak (Mode)" "Max" "" "Min" "Peak (Mode)" "Max" "Min" "Peak (Mode)" "Max" "Min" "Peak (Mode)" "Max" "Min" "Peak (Mode)" "Max" "Min" "Peak (Mode)" "Max"},OFile_LCOJF_MainResults,'Sheet','Costs','Range',strcat('AX',int2str(3)));

    % Add headings to the "Emissions" sheet
    writecell({"" "ElectricityEmissions" "" "" ""  "" "" "DAC Emissions" "" "" ""  "" "" "" "" "CombustionEmissions" "" "DACElectricity" "" "" ""  "" "" "BlueH2Emissions" "" "" "GreenH2Emissions" "" ""},OFile_LCOJF_MainResults,'Sheet','Emissions','Range',strcat('A',int2str(2)));
    writecell({"EmissionsCosts" "" ""},OFile_LCOJF_MainResults,'Sheet','Emissions','Range',strcat('AD',int2str(1)));
    writecell({"Year" "BC" "ON" "QC" "AB"  "MB" "NB" "BC" "ON" "QC" "AB"  "MB" "NB" "" "" "Combustion Emissions" "" "BC" "ON" "QC" "AB"  "MB" "NB"  "Min" "Peak (Mode)" "Max" "Min" "Peak (Mode)" "Max" "Min" "Peak (Mode)" "Max"},OFile_LCOJF_MainResults,'Sheet','Emissions','Range',strcat('A',int2str(3)));
    
    % Move generated files to their respective output subfolders
    movefile(OFile_LCOJF_CDFPerProvince,fullfile(outputCDFFolder,OFile_LCOJF_CDFPerProvince))
    movefile(OFile_LCOJF_MainResults,fullfile(outputLcojfFolder ,OFile_LCOJF_MainResults))


    disp("Done adding headings")
end 

function xlcol_addr=num2xlcol(col_num)
% col_num - positive integer greater than zero
    n=1;
    while col_num>26*(26^n-1)/25
        n=n+1;
    end
    base_26=zeros(1,n);
    tmp_var=-1+col_num-26*(26^(n-1)-1)/25;
    for k=1:n
        divisor=26^(n-k);
        remainder=mod(tmp_var,divisor);
        base_26(k)=65+(tmp_var-remainder)/divisor;
        tmp_var=remainder;
    end
    xlcol_addr=char(base_26); % Character vector of xlcol address
end