% Supporting file for combining per LTE data with environment, generating
% the data-series to be used in the main statistical module

% Functions include:
% Main function: PerLTEcombine()
% Daughters:
% - NormMax(EnvIn, YearChar, Starts)
% - NormMaxTot(EnvIn, YearChar)
% - NormMaxFuture(EnvIn, RefYield, MinRefArray)
% - OpenArray(File, Type)
% - FirstTraitTypeMention(Array)
% - RainCutOffs()


function [CombinedYields,PerYearFunctions,FutureEnvironments] = PerLTEcombine
warning off
Current = pwd;

OutPutFile = 'CombinedYields';
load('InputData.mat','LTEstoRun')
LTEList = (LTEstoRun(:,1)); 
clear LTEstoRun
%
%Calculate yield arrays
disp('Running per LTE arrays')
[~,LTEFullSelected] = GenerateLTEMatFiles;
cd(Current)
Types = [{'MaxYieldTreat'},'AverageYield','MinimumYieldTreat'];

for Var = 1:3
    clc
    disp('Running per Type combination')
    disp(Types(Var))
    ParameterNameList =  [{'DailyTempRange'},'MinTemperatureYear','MaxTemperatureYear','RadiationMax',...
    'CO2ManuaLoa','AnnualPET','AnnualRainfall','CWDMAXperCallenderYear','PopPressure',...
    'MeanOzone','AnnualTempRange','PrecipitationSeasonality','ETOSeasonality','Isothermality','DailyMeanTemperatureYear',...
    'CWDMAXperCallenderYearGlobal','GrowthSeasonLength','DailyMeanTempChange']; 
    % Combine in one file
    Type = Types(Var);

    for LTEname = 1: length(LTEList)
        %disp(List(i))
        File = (char(LTEList(LTEname)));
        cd('MatlabLTEFiles')
        [Array,LTETesterSelect] = OpenArray(File,Type);
        cd(Current)
        if LTEname == 1
            Combined = Array;
            PerLTEYearFunction = dataset({File},'varnames', {'LTE'}); %#ok<*DTSET>    
        else
            Combined = [Combined;Array]; %#ok<*AGROW>
        end
        PerLTEYearFunction.LTE(LTEname,1) = {File};
        ToTest = FirstTraitTypeMention(Array);
        Tester = find(strcmp(LTETesterSelect.Treatment,ToTest.Treatment)==1);
        PerLTEYearFunction.TreatType(LTEname,1) =ToTest.Treatment;
        PerLTEYearFunction.Crop(LTEname,1) =ToTest.CType;
        PerLTEYearFunction.YearConstant(LTEname,1) = LTETesterSelect.YearConstant(Tester);
        PerLTEYearFunction.YearCoef(LTEname,1) = LTETesterSelect.YearCoef(Tester);
        PerLTEYearFunction.Year_FValue(LTEname,1) = LTETesterSelect.Year_FValue(Tester);
        PerLTEYearFunction.Year_PValue(LTEname,1) = LTETesterSelect.Year_PValue(Tester);
        PerLTEYearFunction.Years(LTEname,1) = LTETesterSelect.Years(Tester);
        PerLTEYearFunction.StartYear(LTEname,1) = LTETesterSelect.StartYear(Tester);
        PerLTEYearFunction.EndYear(LTEname,1) = LTETesterSelect.EndYear(Tester);
        clear Array LTETesterSelect
    end
    disp('Adding Environmental info')
    load('Environment.mat');
    load('InputData.mat','GeoLocationsStart')
    LTETesterList = PerLTEYearFunction.LTE;
    AnnualStart = find(strcmp(AnnualPET.Properties.VariableNames,'Y_1970')==1);
    Annual = AnnualPET.Properties.VariableNames(AnnualStart:end);
    B = cellstr(string(Combined.Year));
    for i = 1:size(Combined,1)
        ObservationsList(i,1) = {[char(Combined.FullLTE(i)),'_',char(B(i))]};
    end
    BarFig = waitbar(0, 'Starting');
    for intry = 1:size(Combined,1)
        Year =  Combined.Year(intry);
        YearChar = {['Y_', mat2str(Year)]}; 
        LTETest = Combined.LTE(intry);
        LTETestFull = Combined.FullLTE(intry);
        if intry/25 == round(intry/25)
            waitbar(intry/size(Combined,1), BarFig, sprintf('Adding Environments for %d rows: %d %%', size(Combined,1),floor(intry/size(Combined,1)*100)));
        end
        %pause(0.1);
            % clc
            % Mesesage = ['Adding Environments for row: ',mat2str(intry), ' of total ', mat2str(size(Combined,1)), ' with LTE: ',char(LTETest) ];
            % disp(Mesesage)
        LTENum = find(strcmp(PerLTEYearFunction.LTE,LTETestFull));
        Starts(1) = PerLTEYearFunction.StartYear(LTENum,1);
        Starts(2) = PerLTEYearFunction.EndYear(LTENum,1);
        if Starts(1) < 1970
            Starts(1) = 1970;
        end
        if Starts(2) > 2022
            Starts(2) = 2022;
        end
        % Make Loop based on List
        Combined.Longitude_X(intry,1) = GeoLocationsStart.Longitude_X(LTETest,1);
        Combined.Lattitude_Y(intry,1) = GeoLocationsStart.Lattitude_Y(LTETest,1);
        Combined.Elevation(intry,1) = GeoLocationsStart.Elevation(LTETest,1);
        Combined.Continent(intry,1) = GeoLocationsStart.Continent(LTETest,1);
        Combined.LatRegion(intry,1) = GeoLocationsStart.LatRegion(LTETest,1);
        Combined.Cluster(intry,1) = GeoLocationsStart.Cluster(LTETest,1);
        Combined.CWDChange(intry,1) = GeoLocationsStart.CWDChangeCurrent(LTETest,1);
        Combined.RainFallChange(intry,1) = GeoLocationsStart.RainFallChange(LTETest,1);
        Combined.ID_LTE(intry,1) = GeoLocationsStart.SeqNum(LTETest,1);
        %LTENum = Combined.ID_LTE(intry);

        % Normalisation per LTE
        for x = 1: length(ParameterNameList)
            Combined.(genvarname(char(ParameterNameList(x))))(intry,1) = table2array(eval([char(ParameterNameList(x)),'(LTETest,YearChar)']));
            [Combined.(genvarname(char({[char(ParameterNameList(x)),'Norm']})))(intry,1),NormRef.RefYieldNorm(intry,x),MinValueNorm(intry,x),NormRefLTE(LTENum,x),MinRefLTE(LTENum,x)] =...
                                eval(['NormMax(',char(ParameterNameList(x)),'(LTETest,:),YearChar,Starts);']); %#ok<*ASGLU,*NASGU>
            NormalisedParameterListToDo(x) = {[char(ParameterNameList(x)),'Norm']};
            [Combined.(genvarname(char({[char(ParameterNameList(x)),'NormFullTime']})))(intry,1),NormRef.RefYieldNormTot(intry,x),MinValueNormTot(intry,x)] =...
                                eval(['NormMaxTot(',char(ParameterNameList(x)),'(LTETest,:),YearChar);']); 
            NormalisedParameterListToDoFullTime(x) = {[char(ParameterNameList(x)),'NormFullTime']};
            
            AnnualStartParameter = find(strcmp(eval([char(ParameterNameList(x)),'.Properties.VariableNames']),'Y_1970')==1); 
            ValuesToAddRange = (eval(['NormMaxFuture(',char(ParameterNameList(x)),'(LTETest,AnnualStartParameter:end),NormRef.RefYieldNorm(intry,x),MinValueNorm(intry,x));']));
            ValuesToFullTime =  (eval(['NormMaxFuture(',char(ParameterNameList(x)),'(LTETest,AnnualStartParameter:end),NormRef.RefYieldNormTot(intry,x),MinValueNormTot(intry,x));']));
            PresentEnvironment.Range.(genvarname(char(ParameterNameList(x))))(ObservationsList(intry),:) = array2table(ValuesToAddRange,'VariableNames',Annual);
            PresentEnvironment.FullTime.(genvarname(char(ParameterNameList(x))))(ObservationsList(intry),:) = array2table(ValuesToFullTime,'VariableNames',Annual);
        end
        Combined.RadiationNormPLS(intry,1) = GeoLocationsStart.RadiationNormPLS(LTETest,1);
        NormalisedParameterListToDo(x+1) = {'RadiationNormPLS'};

    end
    delete(BarFig)
    clear Annual B ObservationsList    
    % Future Environment
    Annual = Future_AnnualPET.Properties.VariableNames;
    B = cellstr(string(Combined.Year));
    for i = 1:size(Combined,1)
        ObservationsList(i,1) = {[char(Combined.FullLTE(i)),'_',char(B(i))]};
    end
    BarFig = waitbar(0, 'Starting');
    for intry = 1:size(Combined,1)
        if intry/25 == round(intry/25)
            waitbar(intry/size(Combined,1), BarFig, sprintf('Adding Future Environments for %d rows: %d %%', size(Combined,1),floor(intry/size(Combined,1)*100)));
        end
        %     clc
        %     Mesesage = ['Adding Future Environments for row: ',mat2str(intry), ' of total ', mat2str(size(Combined,1)), ' with LTE: ',char(LTETest) ];
        %     disp(Mesesage)
        % end
        LTETest = Combined.LTE(intry);
        for x = 1: length(ParameterNameList)
            FutureParaName = ['Future_',char(ParameterNameList(x))];
            ValuesToAddRange = (eval(['NormMaxFuture(',char(FutureParaName),'(LTETest,:),NormRef.RefYieldNorm(intry,x),MinValueNorm(intry,x));']));
            ValuesToFullTime =  (eval(['NormMaxFuture(',char(FutureParaName),'(LTETest,:),NormRef.RefYieldNormTot(intry,x),MinValueNormTot(intry,x));'])); 
            FutureEnvironment.Range.(genvarname(FutureParaName))(ObservationsList(intry),:) = array2table(ValuesToAddRange,'VariableNames',Annual);
            FutureEnvironment.FullTime.(genvarname(FutureParaName))(ObservationsList(intry),:) = array2table(ValuesToFullTime,'VariableNames',Annual);
        end
    end
    delete(BarFig)
    NormRef.ReferenceMaxEnvperLTE(LTETesterList,:) = array2table(NormRefLTE,'VariableNames',ParameterNameList);
    NormRef.ReferenceMinEnvperLTE = array2table(MinRefLTE,'VariableNames',ParameterNameList);  
    CombinedYields.(genvarname(char(Types(Var)))) = Combined;
    FutureEnvironments.(genvarname(char(Types(Var)))) = FutureEnvironment;
    PresentEnvironments.(genvarname(char(Types(Var)))) = PresentEnvironment;
    PerYearFunctions.(genvarname(char(Types(Var)))) = PerLTEYearFunction;
    ObservationsNameList.(genvarname(char(Types(Var)))) = ObservationsList;
    NormalisationRefs.(genvarname(char(Types(Var)))) = NormRef;
    RainCWDCutOffs = RainCutOffs;
    save(OutPutFile,'CombinedYields','LTEFullSelected','FutureEnvironments','PresentEnvironments','PerYearFunctions','ObservationsNameList', ...
        'NormalisationRefs','RainCWDCutOffs','NormalisedParameterListToDo','NormalisedParameterListToDoFullTime')  
     clearvars -except OutPutFile LTEList Var Types Current  LTEFullSelected CombinedYields FutureEnvironments PresentEnvironments PerYearFunctions   ObservationsNameList NormalisationRefs
end
end
%%
% Normalisation per LTE!!!
function [OutValue,RefEnv,MinRefArray,RefEnvCopy,MinRefArrayCopy] = NormMax(EnvIn,YearChar,Starts)
Start = find(strcmpi(EnvIn.Properties.VariableNames, ['Y_', mat2str(Starts(1))]));
End = find(strcmpi(EnvIn.Properties.VariableNames, ['Y_', mat2str(Starts(2))]));

% Remove spurious negative values
RefArray = table2array(EnvIn);
RefArray = RefArray(1,Start:End);
MinTest = RefArray((RefArray<0));
if isempty(MinTest) ~= 1
    MinRefArray = min(RefArray);
    RefArrayNorm = RefArray - MinRefArray;
else
     MinRefArray = 0;
     RefArrayNorm = RefArray;
end
% MinRefArray = min(RefArray);
% RefArrayNorm = RefArray - MinRefArray;
% Normalised for the range (min - max)!!, not between 0 and max. So Low values
% are unlikely (e.g. Co2 ranges from 0.76 to 1)
RefEnv = max(RefArrayNorm); 
OutValue =((table2array(EnvIn(1,YearChar)))- MinRefArray)./RefEnv;
if OutValue > 1
  OutValue = 1;
end
RefEnvCopy = RefEnv;
MinRefArrayCopy = MinRefArray;
end

%%
function [OutValue,RefYield,MinRefArray] = NormMaxTot(EnvIn,YearChar)
Start = find(strcmpi(EnvIn.Properties.VariableNames, 'Y_1970'));
End = find(strcmpi(EnvIn.Properties.VariableNames, 'Y_2022'));

% Remove spurious negative values
RefArray = table2array(EnvIn);
RefArray = RefArray(1,Start:End);
MinTest = RefArray((RefArray<0));
if isempty(MinTest) ~= 1
    MinRefArray = min(RefArray);
    RefArrayNorm = RefArray - MinRefArray;
else
     MinRefArray = 0;
     RefArrayNorm = RefArray;
end
% Normalise
RefYield = max(RefArrayNorm); 
OutValue =((table2array(EnvIn(1,YearChar)))- MinRefArray)./RefYield;
if OutValue > 1 %#ok<BDSCI>
    OutValue = 1;
end
end
%% 
function [OutArray] = NormMaxFuture(EnvIn,RefYield,MinRefArray)
RefArray = table2array(EnvIn);
OutArray = (RefArray - MinRefArray)./RefYield;
OutArray(OutArray<0) = 0;
end
%%
function [Array,LTETesterSelect] = OpenArray(File,Type) %#ok<STOUT>
Filename = [File,'.mat'];
load(Filename)
Array = eval(char(Type));
end
%%
function ToTest = FirstTraitTypeMention(Array)
count = 1;
ToTest = [];
while isempty(ToTest)==1
    if strcmp(Array.TreatType(count,1),'No values')~=1
        ToTest.Treatment = Array.TreatType(count,1);
        ToTest.CType = Array.Ctype(count,1);
    end
    count = count + 1;
end
end
%%
% Note since comparing Min and Max the same cut-off is used for both
% Calcuated as 4 way average for their Tstat, taken the highest mean
% 'AverageYield' is seperate
function RainCWDCuttOffs = RainCutOffs
RainCWDCuttOffs.MaxYieldTreat.RainDecrease = 0;
RainCWDCuttOffs.MaxYieldTreat.RainIncrease = 0;
RainCWDCuttOffs.MaxYieldTreat.C3RainDecrease = 0;
RainCWDCuttOffs.MaxYieldTreat.C3RainIncrease = 0;
RainCWDCuttOffs.MaxYieldTreat.C4RainDecrease = 0;
RainCWDCuttOffs.MaxYieldTreat.C4RainIncrease = 0;
RainCWDCuttOffs.MaxYieldTreat.CWDDecrease = 0;
RainCWDCuttOffs.MaxYieldTreat.CWDIncrease = 0;
RainCWDCuttOffs.MaxYieldTreat.C3CWDDecrease =0;
RainCWDCuttOffs.MaxYieldTreat.C3CWDIncrease = 0;
RainCWDCuttOffs.MaxYieldTreat.C4CWDDecrease = 0;
RainCWDCuttOffs.MaxYieldTreat.C4CWDIncrease = 0;

% RainCWDCuttOffs.MaxYieldTreat.RainDecrease = 53.48;
% RainCWDCuttOffs.MaxYieldTreat.RainIncrease = 53.48;
% RainCWDCuttOffs.MaxYieldTreat.C3RainDecrease = 107;
% RainCWDCuttOffs.MaxYieldTreat.C3RainIncrease = 107;
% RainCWDCuttOffs.MaxYieldTreat.C4RainDecrease = 74.83;
% RainCWDCuttOffs.MaxYieldTreat.C4RainIncrease = 74.83;
% RainCWDCuttOffs.MaxYieldTreat.CWDDecrease = -0.05;
% RainCWDCuttOffs.MaxYieldTreat.CWDIncrease = -0.05;
% RainCWDCuttOffs.MaxYieldTreat.C3CWDDecrease = 0.04;
% RainCWDCuttOffs.MaxYieldTreat.C3CWDIncrease = 0.04;
% RainCWDCuttOffs.MaxYieldTreat.C4CWDDecrease = -0.06;
% RainCWDCuttOffs.MaxYieldTreat.C4CWDIncrease = -0.06;

% Note Minimum cut-off is maximum cutt-oof
RainCWDCuttOffs.MinimumYieldTreat = RainCWDCuttOffs.MaxYieldTreat;
RainCWDCuttOffs.AverageYield  = RainCWDCuttOffs.MaxYieldTreat;

% RainCWDCuttOffs.AverageYield.RainDecrease = 107;
% RainCWDCuttOffs.AverageYield.RainIncrease = 107;
% RainCWDCuttOffs.AverageYield.C3RainDecrease = 107;
% RainCWDCuttOffs.AverageYield.C3RainIncrease = 107;
% RainCWDCuttOffs.AverageYield.C4RainDecrease = 61.3;
% RainCWDCuttOffs.AverageYield.C4RainIncrease = 61.3;
% RainCWDCuttOffs.AverageYield.CWDDecrease = -0.06;
% RainCWDCuttOffs.AverageYield.CWDIncrease = -0.06;
% RainCWDCuttOffs.AverageYield.C3CWDDecrease = 0.05;
% RainCWDCuttOffs.AverageYield.C3CWDIncrease = 0.05;
% RainCWDCuttOffs.AverageYield.C4CWDDecrease = -0.03;
% RainCWDCuttOffs.AverageYield.C4CWDIncrease = -0.03;
end