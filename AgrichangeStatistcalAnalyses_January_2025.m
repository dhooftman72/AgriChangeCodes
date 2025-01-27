% Main statistical running file, calcuating and outputting all NLM and PLS
% results. Loads the data generated via the 'PerLTEcombine' code

% Functions iclude:
% Main function: AgrichangeStatistcalAnalyses_January_2025
% Daughters:
% - Mainprocedure(TypesDoing, Types, DataFile)
% - PLSandNLMRegres(X_array, Y_values, PerLTEWeights, NrParameters, AxNum, ~)
% - nlmRegression(XVarIn, YVarIN, Weights, Name, NLMOut, Env, Loop, FirstTime)
% - WhatToOutputPLS(NLM, PLS, Parameters, FullAxName, AxNum, Selected, NameListPLS)
% - NormalisationArrayFunc(Array, InNameList, OutNameList, Observations)
% - FutureNormalisationArrayFun(EnvIn, MinRange, OutNameList, InNameList)
% - ForwardRegression(XVar, InVar, EnvironmentNLM, NameList, Weights, AxNum)
% - PLSFinalRun(PLS, EnvironmentalArray, YieldforPLS, NumbersFixed, NameList, Weights, AxNum, Extra)
% - LinearHeteroskedacityCorrection(XVarIn, YVarIn, Predictions)
% - MeanPerYear(EnvIn, ~, Name, InVar, YearsIn)
% - MeanPerYearYield(Yield, InVar, YearsIn, Name)
% - MeanPerYearFuture(EnvIn, Name, InVar, LTEs)
% - YieldDevelopment(PLS, PerYear, FutureYears, AxNum, ~, ActualYield)
% - YieldDevelopmentCombine(PLSIn, ~, PerYear, NLMNeed, ActualYield)
% - CorrelationMatrixGeneration(VarIn, NameList)
% - RainCWDCutOffFunc(OnetoFind, InPutArray, NumberofObservations, TypesYield, YieldToRun, RainCWDCutOffs, Input)
% - HoghbergFunc(VarIn)
% - ExtractNLMResults(VarIn)
% - ExtractPLSResults(StructIn)
% - Hoghberg(pvals, q, method, report)


function AgrichangeStatistcalAnalyses_January_2025
warning off
%[~,~,~] = PerLTEcombine;

for TheTypeofYieldToRun = [1:3] %#ok<*NBRAK2>
    TypesYield = [{'MaxYieldTreat'},'AverageYield','MinimumYieldTreat'];
    load('CombinedYields.mat','CombinedYields','FutureEnvironments','PerYearFunctions','NormalisedParameterListToDo','RainCWDCutOffs') %#ok<*LOAD>
    % CHANGE if needed 
    RunningPars.TypesDoing = [1:53]; % MonoLimited: [1:19,21:23,25:31,33:34,36:49]; Limited 1:53

    RunningPars.Types = [{'AllData'};{'C3'};{'C4'}; 'Temperate';'Tropical';...%5
        'CWDDecrease';'CWDIncrease';'C3Temperate';'C3Tropical';'C3CWDDecrease';'C3CWDIncrease';... %11
         'C4Temperate';'C4Tropical';'C4CWDDecrease';'C4CWDIncrease';... %15
        'EuropeGlobal';'SubSaharanAfrica';'USMexico'; 'OtherParts';... %19
        'Europe';'Africa';'NorthAmericas';'SouthAmericas';... %23
        'Maize';'Wheat';'Oats';'Barley';'Rice';... %28
        'MaizeTropical'; 'MaizeTemperate';'WheatTropical'; 'WheatTemperate';... %32
        'MaizeSubSaharanAfrica';'MaizeEuropeGlobal';'MaizeUSMexico';'MaizeOtherParts';... %36
        'WheatEuropeGlobal';'WheatUSMexico';'WheatOtherParts';... %39
        'C3SubSaharanAfrica'; 'C3EuropeGlobal';'C3USMexico';'C3OtherParts';... %43
        'C4SubSaharanAfrica'; 'C4EuropeGlobal';'C4USMexico';'C4OtherParts';... %47
        'RainDecrease';'RainIncrease';'C3RainDecrease';'C3RainIncrease';... %51
        'C4RainDecrease';'C4RainIncrease']; %53
    % Types added per LTE
    Types2 = unique(CombinedYields.(genvarname(char(TypesYield(TheTypeofYieldToRun)))).LTE); %54:163
    RunningPars.Types = [RunningPars.Types;Types2];
    clear Types2

    % Name used in outputs; Needs aligment with InputList(NormalisedParameteListToDo)
    NameList = [{'DailyTempRange'},'MinTemp','MaxTemp','Radiation','CO2Norm',... % Remove 4 from Global PLS
        'PET','Rainfall','CWD','PopPressure','Ozone','AnnualTempRange',...   % Remove 9 from all PLS; remove 8 and 10 from Global PLS
        'PrecipitationSeasonality','ETOSeasonality','Isothermality','DailyTemperature',...
        'CWDGlobalPLS','GrowthSeasonLength','RadiationGlobalPLS','DailyMeanTempChange']; % remove 17 & 19 from all PLS; 

    % Clean out PopPressure (9), Ozone (10), Temp change (19), GrowthSeasonLength (17)
    % and potentially others for PLS
    % Note I changed everything to the global parameters, so CWD (8) and
    % Radiation (4) are the ones for global in all cases 
    %RunningPars.ToRemove = [4,8,9,10,17,19];
    
    % Test without ranges and seasonality (except isotheminality 
    % - 1,11,12,13
    %RunningPars.ToRemove = [1,4,8,9,10:13,17,19]; % No ranges test

     % Test without PET, since it is covered in full CWD - 6 & 13
    RunningPars.ToRemove = [4,6,9,10,16,13,17,19]; % NO PET test 2401_v2
   
    RunningPars.closeparalel = 0;
    RunningPars.NLM_report = 1;
    RunningPars.FirstTime = 1;
    % Close the parallel pool if present
    poolobj = gcp('nocreate');
    delete(poolobj);
    clear poolobj 
    delete('StartingData.mat')

    OutputFileName = ['Outputs_',char(TypesYield(TheTypeofYieldToRun))];
    OutputFileNameInBetween = [OutputFileName,'Store.mat'];
    clc
    disp(TypesYield(TheTypeofYieldToRun))

    WhichYield = 'MedianYieldVariety' ;

    InPutArray = CombinedYields.(genvarname(char(TypesYield(TheTypeofYieldToRun))));
    FutureEnvironmentsToUse = FutureEnvironments.(genvarname(char(TypesYield(TheTypeofYieldToRun)))).Range;
    PerLTETimesSeries = PerYearFunctions.(genvarname(char(TypesYield(TheTypeofYieldToRun))));
    Tester = unique(InPutArray.FullLTE);
    for i = 1:length(Tester)
        Values = find(strcmp(InPutArray.FullLTE,Tester(i))==1);
        WhichOne = find(strcmp(PerLTETimesSeries.LTE,Tester(i))==1);
        InPutArray.YearsRunning(Values,1) = PerLTETimesSeries.Years(WhichOne);
    end
    Fields = fieldnames(FutureEnvironmentsToUse);
    clear Values Tester WhichOne i CombinedYields  FutureEnvironments PerYearFunctions 

    % Remove No yield data values (environmental data are prechecked for
    % data containment and corrected)
    IsnanList = find(isnan(InPutArray.(genvarname(WhichYield)))==1); %#ok<*COMPNOP>
    CropSpecies = cellstr(InPutArray.Crop);
    IsnanList = unique([IsnanList; (find(strcmp(CropSpecies,'Soybean')))]);
    IsnanList = unique([IsnanList; (find(strcmp(CropSpecies,'Casava')))]);
    IsnanList = unique([IsnanList; (find(strcmp(CropSpecies,'SugarBeet')))]);
    IsnanList = unique([IsnanList; (find(strcmp(CropSpecies,'Rapeseed')))]);
    IsnanList = unique([IsnanList; (find(InPutArray.(genvarname(WhichYield))<=0))]);
    InPutArray(IsnanList,:) = [];
    for slen = 1: length(Fields)
        FutureEnvironmentsToUse.(genvarname(char(Fields(slen))))(IsnanList,:) = [];
    end
    Years = InPutArray.Year - (min(InPutArray.Year)-1);
    Observations = FutureEnvironmentsToUse.(genvarname(char(Fields(1)))).Properties.RowNames;
    clear IsnanLis* slen CropSpecies Fields

    NrNamesList = length(NameList);
    for i = 1:NrNamesList
        Test = char(NormalisedParameterListToDo(i));
        Test = Test(1:(end-4));
        NameListFuture(i) = {['Future_',Test]};
    end
    clear i Test

    disp('Renormalising Present arrays')
    [~,EnvironmentalArray,NormalisationArrayGlobal,~] =...
        NormalisationArrayFunc(InPutArray,NormalisedParameterListToDo(1:end),NameList(1:end),Observations);  %EnvironmentalOrgArray -2
    %EnvironmentalArray.(genvarname(char(NameList(NrNamesList-1)))) = InPutArray.ElevationNorm;
    EnvironmentalArray.(genvarname(char(NameList(NrNamesList)))) = InPutArray.RadiationNormPLS;

    %EnvironmentalArray = EnvironmentalOrgArray;
    % Since this is already normalised in PetLTECombine with max = 1,
    % this function only stretches the minimum to 0: so range is 0-1
    % Note that the names are changed from incl. 'Norm' in
    % NormalisedParameterListToDo to Without following Namelist

    disp('Renormalising and smoothening Future arrays')
    FutureEnvironmentalArray = FutureNormalisationArrayFun(FutureEnvironmentsToUse,NormalisationArrayGlobal.Min,NameListFuture,NameList(1:end-1)); 
    Annual = FutureEnvironmentalArray.(genvarname(char(NameListFuture(1)))).Properties.VariableNames(1:end);
    % RepeatArray = repmat(InPutArray.ElevationNorm,1,length(Annual));
    % FutureEnvironmentalArray.(genvarname(char(NameListFuture(NrNamesList-1)))) = array2table(RepeatArray,'VariableNames',Annual,'RowNames',Observations);
    RepeatArray = repmat(InPutArray.RadiationNormPLS,1,length(Annual));
    FutureEnvironmentalArray.(genvarname(char(NameListFuture(NrNamesList)))) = array2table(RepeatArray,'VariableNames',Annual,'RowNames',Observations);
    % Future is normalised with same minumum as present data, but keeps its
    % own maximum (as it may exceed the range of present). It is also
    % smoothend with a 5-year average to generate better trends
    clear NormalisedParameterListToDo FutureEnvironmentsToUse Annual RepeatArray

    Yield = InPutArray.(genvarname(WhichYield));
    YieldNormMedian =  InPutArray.MedianYieldNorm;
    NumberofObservations = length(Yield);
    if NumberofObservations ~= length(Observations)
        disp(' ')
        disp(' ')
        disp('FATAL ERROR')
        disp('Yield and Environmental Arrays are not aligned')
        disp('Check lengths in PerLTECombine outputs and Correct first')
        save('AllParametersatTimeofFailure')
        return
    end
    CropType = cellstr(InPutArray.Ctype); %#ok<*NASGU>
    CropSpecies = cellstr(InPutArray.Crop);
    Continent = cellstr(InPutArray.Continent);
    Cluster =  cellstr(InPutArray.Cluster);
    Test = find(strcmp(Cluster,'Europe') ==1);
    Cluster(Test) = {'EuropeGlobal'};
    clear Test

    LTEsIndividual = cellstr(InPutArray.LTE);

    % Reset Area names where needed
    List = unique([(find(strcmp(Continent,'WestAfrica'))); (find(strcmp(Continent,'EastAfrica')))]);
    Continent(List) = {'Africa'};
    clear List
    LatRegion = cellstr(InPutArray.LatRegion);
    List = find(strcmp(LatRegion,'(sub)Tropical'));
    LatRegion(List) = {'Tropical'};
    clear List
    List = unique([(find(strcmp(Cluster,'Asia'))); (find(strcmp(Cluster,'Australia')));(find(strcmp(Cluster,'Brazil'))) ]);
    Cluster(List) = {'OtherParts'};
    clear List
    List = find(strcmp(Cluster,'US&Mexico'));
    Cluster(List) = {'USMexico'};
    clear List
    

    Weights = InPutArray.YearsRunning;
    YearsOrg = InPutArray.Year;
    LTEs = InPutArray.FullLTE;    
    load('CombinedYields.mat','NormalisationRefs')
    NormalisationRefsPerLTE.Max = NormalisationRefs.(genvarname(char(TypesYield(TheTypeofYieldToRun)))).ReferenceMaxEnvperLTE;
    NormalisationRefsPerLTE.Min = NormalisationRefs.(genvarname(char(TypesYield(TheTypeofYieldToRun)))).ReferenceMinEnvperLTE;

    save('StartingData')
    [NLM, PLS,HeteroCorrected, PerYear,FutureYears,CorrelationMatrix,FullExtractedStats,ExtractedPLS,MedianYield] = Mainprocedure(RunningPars.TypesDoing,RunningPars.Types,'StartingData.mat');
    save(OutputFileName, 'NLM', 'PLS','HeteroCorrected', 'PerLTETimesSeries','InPutArray', 'PerYear','FutureYears','CorrelationMatrix',...
        'NormalisationArrayGlobal','NormalisationRefsPerLTE','FullExtractedStats','ExtractedPLS','MedianYield');
    delete(OutputFileNameInBetween)
    clearvars -except TheTypeofYieldToRun
end % for crop yield typesTheTypeofYieldToRun
end % End of main function
%%
function [NLM, PLS,HeteroCorrected, PerYear,FutureYears,CorrelationMatrix,FullExtractedStats,ExtractedPLS,MedianYield] = Mainprocedure(TypesDoing,Types,DataFile)
load(DataFile)
for Env = 1:1:(NrNamesList) % +2 if also years are to be included
    load('StartingData.mat','Yield')
    if Env == (size(EnvironmentalArray,2) + 1)  % So years is done after all the other ones.
        EnvironVar = Years; % Note Year is not Normalised
        Name = {'Years'};
    elseif Env == (size(EnvironmentalArray,2) + 2)
         EnvironVar = Years; % Note Year is not Normalised
         Yield = YieldNormMedian;
         Name = {'YearsNotVarityCorrected'};   
    else
        EnvironVar =table2array(EnvironmentalArray(:,Env));   
        FutEnviron =(FutureEnvironmentalArray.(genvarname(char(NameListFuture(Env)))));
        Name = NameList(Env);
    end
    if RunningPars.NLM_report == 1
        disp(['Doing NLM variable: ',char(Name)])
    end
    for Type = TypesDoing  
        OnetoFind = Types(Type);
        TypeString = char(OnetoFind);
        % Add the increase/decrease cut-offs here as function
        if RunningPars.CutOffDetection ~= 1
            if (Type >=6 && Type <= 7) ||(Type >=10 && Type <= 11)||(Type >=14 && Type <= 15)||(Type >=48 && Type <= 53)
                ChangeType = RainCWDCutOffFunc(OnetoFind,InPutArray,NumberofObservations,TypesYield,TheTypeofYieldToRun,RainCWDCutOffs,'Run it with preset cut-offs'); %#ok<*USENS>
            end
        end
        StringCutoff = 2;
        if isempty(find(strcmp(TypeString(1:2),'Ma') ==1)) ~= 1 || isempty(find(strcmp(TypeString(1:2),'Wh') ==1)) ~= 1
            StringCutoff = 5;
        end
        Typea = TypeString(1:StringCutoff);
        Typeb = TypeString((StringCutoff+1):end);
        Listb = [];

        if Type == 1
            Lista = 1:NumberofObservations;
        elseif Type <= 3
            Lista = find(strcmp(CropType,OnetoFind)==1);
        elseif Type <= 5
            Lista = find(strcmp(LatRegion,OnetoFind)==1);
        elseif Type <= 7
            Lista = find(strcmp(ChangeType,OnetoFind)==1); % CWD
        elseif (Type >= 8 && Type <= 9) || (Type >= 12 && Type <= 13)
            Lista = find(strcmp(CropType,{Typea})==1);
            Listb = find(strcmp(LatRegion,{Typeb})==1);
        elseif (Type >= 10 && Type <= 11) || (Type >= 14 && Type <= 15)    
            Lista = find(strcmp(CropType,{Typea})==1);
            Listb = find(strcmp(ChangeType,{Typeb})==1); % CWD
       elseif Type <= 19
            Lista = find(strcmp(Cluster,OnetoFind)==1);
        elseif Type <= 23
            Lista = find(strcmp(Continent,OnetoFind)==1);
        elseif Type <=28
            Lista = find(strcmp(CropSpecies,OnetoFind)==1);
        elseif Type <=32
            Lista = find(strcmp(CropSpecies,{Typea})==1);
            Listb = find(strcmp(LatRegion,{Typeb})==1);
        elseif Type <=39
            Lista = find(strcmp(CropSpecies,{Typea})==1);
            Listb = find(strcmp(Cluster,{Typeb})==1);
        elseif Type <= 47
            Lista = find(strcmp(CropType,{Typea})==1);
            Listb = find(strcmp(Cluster,{Typeb})==1);
        elseif Type <=49
            Lista = find(strcmp(ChangeType,OnetoFind)==1); % Rain
        elseif Type <= 53
            Lista = find(strcmp(CropType,{Typea})==1);
            Listb = find(strcmp(ChangeType,{Typeb})==1); % Rain
        elseif Type > 53
            Lista = find(strcmp(LTEsIndividual,OnetoFind)==1);
        end
        if isempty(Listb)==1
            Listb = Lista;
        end
        List = intersect(Lista,Listb);
        EnvironType = EnvironVar(List);
        FutEnvironType = FutEnviron(List,:);
        LTEsType = LTEs(List);
        YieldType = Yield(List);
        WeightType = Weights(List);
        YearsType = YearsOrg(List);
        clear List Typea Typeb Lista Listb
        
        % Those need to be used later, so store under Type name
        YieldperType.(genvarname(char(OnetoFind))).(genvarname(char(Name))) = YieldType;
        WeightperType.(genvarname(char(OnetoFind))) = WeightType;
        YearsperType.(genvarname(char(OnetoFind))) = YearsType;
        clear YieldType WeightType YearsType
        %PerYear.(genvarname(char(OnetoFind))).MedianYield = MeanPerYearYield(YieldperType.(genvarname(char(OnetoFind))),YearsperType.(genvarname(char(OnetoFind))));
        if Env == 1
            % Initiate datasets for NLM, PerYear, HeteroCorrected, FutureYears. Note no data is written
            % here
            NLM.UnCorrected.(genvarname(char(OnetoFind))) = dataset({'dummy'}, 'varnames',{'Parameters'});
            NLM.HeteroCorrected.(genvarname(char(OnetoFind))) = NLM.UnCorrected.(genvarname(char(OnetoFind)));
            NLM.HeteroCorrectedRunning.(genvarname(char(OnetoFind))) = NLM.UnCorrected.(genvarname(char(OnetoFind)));
            HeteroCorrected.(genvarname(char(OnetoFind))).Yield = dataset(1, 'varnames',Name);
            HeteroCorrected.(genvarname(char(OnetoFind))).Environment =  HeteroCorrected.(genvarname(char(OnetoFind))).Yield;
            HeteroCorrected.(genvarname(char(OnetoFind))).RunningYield =  HeteroCorrected.(genvarname(char(OnetoFind))).Yield;
            HeteroCorrected.(genvarname(char(OnetoFind))).RunningEnvironment =  HeteroCorrected.(genvarname(char(OnetoFind))).Yield;
            PerYear.(genvarname(char(OnetoFind))).UnCorrectedEnvironment.Median = dataset(([1970:2022])','varnames','Year'); 
            PerYear.(genvarname(char(OnetoFind))).UnCorrectedEnvironment.STD = PerYear.(genvarname(char(OnetoFind))).UnCorrectedEnvironment.Median;
            PerYear.(genvarname(char(OnetoFind))).UnCorrectedEnvironment.NPerYear =  PerYear.(genvarname(char(OnetoFind))).UnCorrectedEnvironment.Median;
            PerYear.(genvarname(char(OnetoFind))).HeteroCorrectedEnvironment = PerYear.(genvarname(char(OnetoFind))).UnCorrectedEnvironment;
            MedianYield.(genvarname(char(OnetoFind))) = PerYear.(genvarname(char(OnetoFind))).UnCorrectedEnvironment;
            FutureYears.(genvarname(char(OnetoFind))).Environment.Median = dataset(([2023:2099])','varnames','Year'); 
            FutureYears.(genvarname(char(OnetoFind))).Environment.STD =  FutureYears.(genvarname(char(OnetoFind))).Environment.Median;
            NLM.YearCorrelationsRegress =  dataset({'dummy'}, 'varnames',{'Parameters'});
        end
        % Write data to PerYear, FutureYears, HeterCorrected, NLM
        MedianYield.(genvarname(char(OnetoFind))) = MeanPerYearYield(YieldperType.(genvarname(char(OnetoFind))).(genvarname(char(Name))),MedianYield.(genvarname(char(OnetoFind))),...
        YearsperType.(genvarname(char(OnetoFind))),Name);
        PerYear.(genvarname(char(OnetoFind))).UnCorrectedEnvironment =...
            MeanPerYear(EnvironType,YieldperType.(genvarname(char(OnetoFind))).(genvarname(char(Name))),Name,PerYear.(genvarname(char(OnetoFind))).UnCorrectedEnvironment,YearsperType.(genvarname(char(OnetoFind))));
        FutureYears.(genvarname(char(OnetoFind))).Environment = MeanPerYearFuture(FutEnvironType,Name,FutureYears.(genvarname(char(OnetoFind))).Environment,LTEsType);

        [~, PredictedValues, NLM.UnCorrected.(genvarname(char(OnetoFind)))] =...
            nlmRegression(EnvironType,YieldperType.(genvarname(char(OnetoFind))).(genvarname(char(Name))),WeightperType.(genvarname(char(OnetoFind))),Name, NLM.UnCorrected.(genvarname(char(OnetoFind))),Env,'No Loop',RunningPars.FirstTime);
        RunningPars.FirstTime = 2;
        [EnvironVarCorType, YieldCorType,EnvironVarCorRunningType,YieldCorRunningType] = LinearHeteroskedacityCorrection(EnvironType,YieldperType.(genvarname(char(OnetoFind))).(genvarname(char(Name))),PredictedValues);
        % Save the values corrected for later depictions
        HeteroCorrected.(genvarname(char(OnetoFind))).Yield.(genvarname(char(Name))) = YieldCorType;
        HeteroCorrected.(genvarname(char(OnetoFind))).Environment.(genvarname(char(Name))) = EnvironVarCorType;
        PerYear.(genvarname(char(OnetoFind))).HeteroCorrectedEnvironment =...
            MeanPerYear(EnvironVarCorType,YieldperType.(genvarname(char(OnetoFind))).(genvarname(char(Name))),Name,PerYear.(genvarname(char(OnetoFind))).HeteroCorrectedEnvironment,YearsperType.(genvarname(char(OnetoFind))));
        if Env <= NrNamesList
            HeteroCorrected.(genvarname(char(OnetoFind))).EnvironmentForPLS(:,Env) = EnvironVarCorType;
        end

        % Base on 10-point running averages
        HeteroCorrected.(genvarname(char(OnetoFind))).RunningYield.(genvarname(char(Name))) = YieldCorRunningType;
        HeteroCorrected.(genvarname(char(OnetoFind))).RunningEnvironment.(genvarname(char(Name))) = EnvironVarCorRunningType;
        RunningWeights = ones(length(YieldCorRunningType),1);
        if Env <= NrNamesList
            HeteroCorrected.(genvarname(char(OnetoFind))).RunningEnvironmentForPLS(:,Env) = EnvironVarCorRunningType;
        end
        [~, ~,NLM.HeteroCorrected.(genvarname(char(OnetoFind)))] =...
            nlmRegression(EnvironVarCorType, YieldCorType,WeightperType.(genvarname(char(OnetoFind))),Name,NLM.HeteroCorrected.(genvarname(char(OnetoFind))),Env,RunningPars.FirstTime);
        [~, ~, NLM.HeteroCorrectedRunning.(genvarname(char(OnetoFind)))] =...
            nlmRegression(EnvironVarCorRunningType,YieldCorRunningType,RunningWeights,Name, NLM.HeteroCorrectedRunning.(genvarname(char(OnetoFind))),Env,RunningPars.FirstTime);

        clear EnvironType LTEsType
    end
    clear list
    % % Regression between the environmental parameters and time!
    [~, ~,NLM.YearCorrelationsRegress] = nlmRegression((Years-min(Years)),EnvironVar,Weights,Name,NLM.YearCorrelationsRegress,Env,RunningPars.FirstTime);
end
% Hogberg corrections 
for Type = TypesDoing
    OnetoFind = Types(Type);
    NLM.HeteroCorrected.(genvarname(char(OnetoFind))).Hochberg = HoghbergFunc(NLM.HeteroCorrected.(genvarname(char(OnetoFind))).PValue); %Per type
    NLM.HeteroCorrected.(genvarname(char(OnetoFind))).HochbergRestricted([2,3,4,5,7,8,9,10,14,15,19],1) = HoghbergFunc(NLM.HeteroCorrected.(genvarname(char(OnetoFind))).PValue([2,3,4,5,7,8,9,10,14,15,19])); %Per type
    NLM.HeteroCorrected.(genvarname(char(OnetoFind))).HochbergRestricted([1,6,11,12,13,16,17,18,20,21],1) = -1;
end
NLM.YearCorrelationsRegress.Hochberg = HoghbergFunc(NLM.YearCorrelationsRegress.PValue); %between the environmental parameters and time!
[FullExtractedStats] = ExtractNLMResults(NLM);

%%
% Close the parallel pool if open
if RunningPars.closeparalel == 1
    disp(' ')
    poolobj = gcp('nocreate');
    delete(poolobj);
end
if RunningPars.NLM_report == 1
    disp(' ')
    disp('Doing CorrelationMatrices')
end
% Correlation Matrices
NameListCorr = NameList;
NameListCorr(end+1) = {'Time'};
for Type = TypesDoing
    OnetoFind = Types(Type);
    CorrEnter = HeteroCorrected.(genvarname(char(OnetoFind))).EnvironmentForPLS;
    CorrEnter(:,end+1) = YearsperType.(genvarname(char(OnetoFind)))-min(YearsperType.(genvarname(char(OnetoFind))));
    CorrelationMatrix.(genvarname(char(OnetoFind))) = CorrelationMatrixGeneration(CorrEnter,NameListCorr);
    clear CorrEnter
end
if RunningPars.NLM_report == 1
    disp(' ')
    disp('Starting PLS calculations')
end

for Type = TypesDoing
    OnetoFind = Types(Type);
    NLMTypesPLS = NLM.HeteroCorrected.(genvarname(char(OnetoFind)));
    EnvironmentforPLS = HeteroCorrected.(genvarname(char(OnetoFind))).EnvironmentForPLS;
    NameListPLS = NameList;
    if Type > 53 % It is just one LTE location so no Radiation differences (18) per site among years
        RunningPars.ToRemove = sort([RunningPars.ToRemove,18]); 
    end
    NLMTypesPLS(RunningPars.ToRemove,:) = [];
    EnvironmentforPLS(:,RunningPars.ToRemove) = [];
    NameListPLS(RunningPars.ToRemove) = [];

     if RunningPars.NLM_report == 1
        disp(' ')
        disp(['Doing Type: ', char(OnetoFind)])
     end
    YieldNeed = YieldperType.(genvarname(char(OnetoFind))).(genvarname(char(NameList(1)))); % Variety corrected Yield
    EnvironNeed =EnvironmentforPLS;
    WeightNeed = WeightperType.(genvarname(char(OnetoFind)));
    NLMNeed = NLMTypesPLS;
    EqualWeightsNeed = ones(length(YieldNeed),1);
    for AxNum = 1:2
        % Full regression PLS
        [NLMOuts,PLSOuts,~] = PLSandNLMRegres(EnvironNeed,YieldNeed,EqualWeightsNeed,3,1,0); 
        [TmpPLS.Overview,TmpPLS.PerEnv.FullAxis1,EstimatedValues,~] = WhatToOutputPLS(NLMOuts,PLSOuts,NameListPLS,{'Full_Axis'},AxNum ,'Full',NameListPLS);
        [~,TmpPLS.PerEnv.FullAxis2,~,~] = WhatToOutputPLS(NLMOuts,PLSOuts,NameListPLS,{'Full_Axis'},AxNum ,'Full',NameListPLS);
        TmpPLS.Values.XFull = EstimatedValues.X;
        TmpPLS.Values.YFull = EstimatedValues.Y;
        % Forward regression PLS
        [NumbersFixed, ~,~] = ForwardRegression(EnvironNeed,YieldNeed,NLMNeed,NameListPLS,WeightNeed,AxNum);
        TmpPLS = PLSFinalRun(TmpPLS,EnvironNeed,YieldNeed,NumbersFixed,NameListPLS,WeightNeed,AxNum,'No Extra');
        TmpPLS = YieldDevelopment(TmpPLS, PerYear.(genvarname(char(OnetoFind))).HeteroCorrectedEnvironment,FutureYears.(genvarname(char(OnetoFind))).Environment,AxNum,'No Extra',MedianYield.(genvarname(char(OnetoFind))).Median(:,2));
        clear NumbersFixed NLMOuts PLSOuts EstimatedValues
    end
    clear *Need
    TmpPLS = YieldDevelopmentCombine(TmpPLS,'No Extra',PerYear.(genvarname(char(OnetoFind))).HeteroCorrectedEnvironment,NLMTypesPLS,MedianYield.(genvarname(char(OnetoFind))).Median(:,2));
    PLS.(genvarname(char(OnetoFind))) = TmpPLS;
    clear TmpPLS
    save(OutputFileNameInBetween, 'NLM', 'PLS','HeteroCorrected', 'PerYear','FutureYears','PerLTETimesSeries','CorrelationMatrix','NormalisationArrayGlobal','NormalisationRefsPerLTE','InPutArray','FullExtractedStats','MedianYield'); 
end
ExtractedPLS = ExtractPLSResults(PLS);
end
%%


%%
function [VarOut,PLS,VIPNum] = PLSandNLMRegres(X_array,Y_values,PerLTEWeights,NrParameters,AxNum,~)
nrAxis = NrParameters;
if nrAxis > 3
    nrAxis = 3;
end
if nrAxis >= 2 % not only is this two axes, it is rotated based on maximising explaining Y-variation
    [XL,yl,PLS.Axis,PLS.YAxis,PLS.BETA,PLS.PCTVAR,~,stats] = plsregress(X_array,Y_values,nrAxis);
    PLS.Weights = stats.W;
    absAxisWeights = abs(PLS.Weights(:,AxNum));
    SumAxisweights = sum(absAxisWeights);
    PLS.AxisWeightsNormalised = absAxisWeights./SumAxisweights;
else
    PLS.Axis = X_array;
end

[~, PredictedValues,~] = nlmRegression(PLS.Axis(:,1),Y_values,PerLTEWeights,{'Dummy'},{'Dummy'},1,2);
[EnvironVarCorType, YieldCorType,~,~] = LinearHeteroskedacityCorrection(PLS.Axis(:,1),Y_values,PredictedValues);
[VarOut, ~,~] = nlmRegression(EnvironVarCorType, YieldCorType,PerLTEWeights,{'Dummy'},{'Dummy'},1,2);

% Calculate Variable Importance in Projection for PLS Regression. NOT USED
W0 = stats.W ./ sqrt(sum(stats.W.^2,1));
p = size(XL,1);
sumSq = sum(PLS.Axis.^2,1).*sum(yl.^2,1);
vipScore = sqrt(p* sum(sumSq.*(W0.^2),2) ./ sum(sumSq,2));
VIPNum = find(vipScore >= 1);
end

%%
function [VarOut, PredictedValues,NLMOut] = nlmRegression(XVarIn,YVarIN,Weights,Name,NLMOut,Env,FirstTime)
% This function is not only used as factor regression but also as
% heteroskedacity correlation function
modelFun = @(b,x) b(1)+ b(2).*x;
start = [1; 0];

VarOut = fitnlm(XVarIn,YVarIN,modelFun,start,'Weight',Weights);
Constant = table2array(VarOut.Coefficients(1,1));
Rico = table2array(VarOut.Coefficients(2,1));
PredictedValues = Constant + Rico.*XVarIn;
if isa(NLMOut,'dataset')== 1
    NLMOut.Parameters(Env,1) = Name;
    NLMOut.PValue(Env,1) = table2array(VarOut.Coefficients(2,4));
    NLMOut.Rsquared(Env,1) = VarOut.Rsquared.Ordinary(1,1);
    NLMOut.Coefficient(Env,1) = Rico;
    NLMOut.Constant(Env,1) = Constant;
    NLMOut.CoefficientSE(Env,1) = table2array(VarOut.Coefficients(2,2));
    NLMOut.ConstantSE(Env,1) = table2array(VarOut.Coefficients(1,2));
    NLMOut.Tstat(Env,1) = table2array(VarOut.Coefficients(2,3));
    NLMOut.LogLikelihood(Env,1) = VarOut.LogLikelihood(1,1);
    NLMOut.AIC(Env,1) = (-2*VarOut.LogLikelihood(1,1)) - (2);
    NLMOut.N(Env,1) = VarOut.NumObservations(1,1);
end
end



%%
function [OutOverall,OutParameter,EstimatedValues,OutParameterFull] = WhatToOutputPLS(NLM,PLS,Parameters,FullAxName,AxNum,Selected,NameListPLS)
EstimatedValues.X = PLS.Axis(:,1);
EstimatedValues.Y = PLS.YAxis(:,1);
LengthArray = length(EstimatedValues.X);
if strcmp(Selected,'Selected')
    if LengthArray > 10
        ToSorted(:,1) = 1:LengthArray;
        ToSorted(:,2) = EstimatedValues.X;
        ToSorted(:,3) = EstimatedValues.Y;
        Sorted = sortrows(ToSorted,2);
        Count = 1;
        for i = 10:1:LengthArray
            Running(Count,1) = Sorted(i,1);
            Running(Count,2) = nanmean(Sorted(i,2)); %#ok<*AGROW,*SAGROW> % X value
            Running(Count,3) = nanmean(Sorted(i-9:i,3)); %Running mean new Values
            Count = Count + 1;
        end
        ToRunning = sortrows(Running,1);
        EstimatedValues.XRunning = ToRunning (:,2);
        EstimatedValues.YRunning = ToRunning (:,3);
        clear ToSorted Sorted Running ToRunning

        EstimatedValues.XHetero = PLS.XVarHetero;
        EstimatedValues.YHetero = PLS.YVarHetero;
        ToSorted(:,1) = 1:LengthArray;
        ToSorted(:,2) = EstimatedValues.XHetero;
        ToSorted(:,3) = EstimatedValues.YHetero;
        Sorted = sortrows(ToSorted,2);
        Count = 1;
        for i = 10:1:LengthArray
            Running(Count,1) = Sorted(i,1);
            Running(Count,2) = nanmean(Sorted(i,2)); %#ok<*AGROW,*SAGROW> % X value
            Running(Count,3) = nanmean(Sorted(i-9:i,3)); %Running mean new Values
            Count = Count + 1;
        end
        ToRunning = sortrows(Running,1);
        EstimatedValues.XHeteroRunning = ToRunning (:,2);
        EstimatedValues.YHeteroRunning = ToRunning (:,3);
    else
        EstimatedValues.XHeteroRunning = EstimatedValues.X;
        EstimatedValues.YHeteroRunning = EstimatedValues.Y;
        EstimatedValues.XHetero = EstimatedValues.X;
        EstimatedValues.YHetero = EstimatedValues.Y;
        EstimatedValues.XRunning = EstimatedValues.X;
        EstimatedValues.YRunning = EstimatedValues.Y;
    end
end
OutOverall = dataset(FullAxName, 'varnames',{'Parameters'});
OutOverall.Tstat(1,1) = table2array(NLM.Coefficients(2,3));
OutOverall.PValue(1,1) = table2array(NLM.Coefficients(2,4));
OutOverall.ExplainedAxis1X(1,1) = PLS.PCTVAR(1,1);
OutOverall.ExplainedAxis1Y(1,1) = PLS.PCTVAR(2,1);
OutOverall.ExplainedAxis2X(1,1) = PLS.PCTVAR(1,2);
OutOverall.ExplainedAxis2Y(1,1) =  PLS.PCTVAR(2,2);
OutOverall.Rsquared(1,1) = NLM.Rsquared.Ordinary(1,1);
OutOverall.LogLikelihood(1,1) = NLM.LogLikelihood(1,1);
OutOverall.AIC(1,1) = (-2*NLM.LogLikelihood(1,1)) - (2 *(length(Parameters)));
OutOverall.Coefficient(1,1) = table2array(NLM.Coefficients(2,1));
OutOverall.Constant(1,1) = table2array(NLM.Coefficients(1,1));
OutOverall.CoefficientSE(1,1) = table2array(NLM.Coefficients(2,2));
OutOverall.ConstantSE(1,1) = table2array(NLM.Coefficients(1,2));
OutOverall.N(1,1) = LengthArray;

OutParameter = dataset((Parameters(1:length(Parameters)))', 'varnames',{'Parameters'});
OutParameter.Weightfactors = PLS.Weights(:,AxNum);
absWeights = abs(OutParameter.Weightfactors);
SumWeights = sum(absWeights);
OutParameter.RelativeWeights = absWeights./SumWeights;
OutParameter.Betas = PLS.BETA(2:end,1); 
% Note there is only one Beta list. That means that Betas for ax 1 and 2
% are identical, unless stochasticty sends them differently among
% calculations. This because there is a small margin of error built into
% the foward regression loop.
OutParameter.Betas(length(OutParameter.Betas)+1) =PLS.BETA(1); % Constant
OutParameter.Parameters(length(OutParameter.Betas)) = {'Constant'};
OutParameter.Weightfactors(length(OutParameter.Betas)) = OutParameter.Betas(length(OutParameter.Betas));

% Include zeros for non selected parameters
OutParameterFull = dataset({'dummty'},'Varnames','Parameters');
for x = 1:length(NameListPLS)
    OutParameterFull.Parameters(x,1) = NameListPLS(x);
    Test = find(strcmp(OutParameter.Parameters,NameListPLS(x)) ==1);
    if isempty(Test) ~= 1
        OutParameterFull.Weightfactors(x,1) = OutParameter.Weightfactors(Test);
        OutParameterFull.RelativeWeights(x,1)= OutParameter.RelativeWeights(Test);
        OutParameterFull.Betas(x,1)= OutParameter.Betas(Test);
    else
        OutParameterFull.Weightfactors(x,1) = 0;
        OutParameterFull.RelativeWeights(x,1)= 0;
        OutParameterFull.Betas(x,1)= 0;
    end
end
OutParameterFull(x+1,:) = OutParameter(length(OutParameter.Parameters),:);

end
%%
function  [VarsOut,EnvironmentalArrayOut,NormalisationArrayGlobal,NonChangeArray] = NormalisationArrayFunc(Array,InNameList,OutNameList,Observations)
% Assuming Array = dataset
% This stretches all between 0 and 1. It does not change the PER LTE
% relative values
NormalisationArrayGlobal.Max = dataset(1,'varnames',OutNameList(1));
NormalisationArrayGlobal.Min = dataset(1,'varnames',OutNameList(1));
NormalisationArrayGlobal.MaxUpdated = dataset(1,'varnames',OutNameList(1));
EnvironmentalArray = [];
NonChangeArray = [];
if isa(Array,'dataset')~= 1
    disp(' ')
    display('Wrong Input data set type in Normalisation Function; correct this ERROR')
    disp(' ')
    save('AllParametersatTimeofFailure')
    ccc
else
    for x = 1:length(OutNameList)
        NormalisationArrayGlobal.Max.(genvarname(char(OutNameList(x)))) = max(Array.(genvarname(char(InNameList(x)))));
        NormalisationArrayGlobal.Min.(genvarname(char(OutNameList(x)))) = min(Array.(genvarname(char(InNameList(x)))));
        Tmp.(genvarname(char(OutNameList(x)))) =  (Array.(genvarname(char(InNameList(x)))) - NormalisationArrayGlobal.Min.(genvarname(char(OutNameList(x)))));
        NormalisationArrayGlobal.MaxUpdated.(genvarname(char(OutNameList(x)))) = max(Tmp.(genvarname(char(OutNameList(x)))));
        VarsOut.(genvarname(char(OutNameList(x)))) = Tmp.(genvarname(char(OutNameList(x))))./max(Tmp.(genvarname(char(OutNameList(x)))));
        NonChangeVarsOut.(genvarname(char(OutNameList(x)))) = Array.(genvarname(char(InNameList(x))));
        EnvironmentalArray = [EnvironmentalArray,VarsOut.(genvarname(char(OutNameList(x))))];
        NonChangeArray =  [NonChangeArray,NonChangeVarsOut.(genvarname(char(OutNameList(x))))];
        %NormalisationArray.MaxCheck.(genvarname(char(OutNameList(x)))) = max(VarsOut.(genvarname(char(OutNameList(x)))));
        %NormalisationArray.MinCheck.(genvarname(char(OutNameList(x)))) = min(VarsOut.(genvarname(char(OutNameList(x)))));
    end
end
EnvironmentalArrayOut(Observations,:) = array2table(EnvironmentalArray,'VariableNames',OutNameList);
end
%%
function EnvironmentOut = FutureNormalisationArrayFun(EnvIn,MinRange,OutNameList,InNameList)
Observations = EnvIn.(genvarname(char(OutNameList(1)))).Properties.RowNames;
Annual = EnvIn.(genvarname(char(OutNameList(1)))).Properties.VariableNames(5:end);
for x = 1:length(InNameList)
    %disp(['Doing Parameter: ',char(ParameterNameList(x))])
    Environment = table2array(EnvIn.(genvarname(char(OutNameList(x)))));
    for row = 1:size(Environment,1)
        MaxRow = max(Environment(row,:),[],2);
        Environtmp = Environment(row,:) - MinRange.(genvarname(char(InNameList(x))));
        Environtmp(Environtmp<0) = 0;
        MaxRowedited = max(Environtmp,[],2);
        Environ = (Environtmp./MaxRowedited).*MaxRow;
        Environ(Environ<0) = 0;
        % Smoothen with 5-year average
        Count = 0;
        for col = 5:length(Environ)
            Count = Count + 1;
            EnvironSmooth(Count) = nanmean(Environ(col-4:col));
        end
        OutArray(row,:) = EnvironSmooth;
        clear Enviro MaxRow MaxRowedited EnvironSmooth
    end
    EnvironmentOut.(genvarname(char(OutNameList(x))))(Observations,:) = array2table(OutArray,'VariableNames',Annual);
end
end
%%
function [NumbersFixedStore, IterationListStore,PLSOutStore] =  ForwardRegression(XVar, InVar, EnvironmentNLM,NameList,Weights,AxNum) 
count = 0;
LookAtStat = 'Tstat';

for testX = 1:size(XVar,2)
    Vtest = XVar(:,testX);
    if all(Vtest == Vtest(1))
        count = count + 1;
        Equals(count) = testX;
    end
end
if count ~= 0
    NameList(Equals) = [];
    XVar(:,Equals) = [];
    EnvironmentNLM(Equals,:) =[];
end
LengthTodo = length(NameList);
SingleFactor = EnvironmentNLM;
YearTest1 = find(strcmp(SingleFactor.Parameters,'Years')==1);
YearTest2 = find(strcmp(SingleFactor.Parameters,'YearsNotVarityCorrected')==1);
YearTest = [YearTest1,YearTest2];
SingleFactor(YearTest,:) = [];
ShortendEnvArray = XVar;

for StartPoint = 1:LengthTodo
    StopLike = 0;
    loop = 1;
    FactorNumbersToRun = 1:LengthTodo;
    Same = 0;
    clear IterationList
    IterationList.Variables(1,1:LengthTodo) = {'empty'};

    MaxLike = SingleFactor.(genvarname(LookAtStat))(StartPoint);
    NumbersFixed = StartPoint;
    MaxLikeStore = MaxLike;
    NrParameterStore = 1;
    VariableIteration = NameList(NumbersFixed);
    IterationList.Variables(loop,1:length(NumbersFixed)) = VariableIteration;
    IterationList.Tstat(loop) = (SingleFactor.Tstat(NumbersFixed));
    IterationList.Rsquared(loop) = SingleFactor.Rsquared(NumbersFixed);
    IterationList.PValue(loop) =SingleFactor.PValue(NumbersFixed);
    IterationList.AIC(loop) =SingleFactor.AIC(NumbersFixed);

    while StopLike == 0
        clear TmpShort
        loop = loop + 1;
        NrParameters = length(NumbersFixed)+1;
        for Env = 1:LengthTodo
            Alreadythere = find(NumbersFixed == Env); %#ok<*EFIND>
            if isempty(Alreadythere)==1
                List = [NumbersFixed, FactorNumbersToRun(Env)];
                Tmp_Array = ShortendEnvArray(:,List);
                clear NLMOuts
                [NLMOuts,PLSOut,~] = PLSandNLMRegres(Tmp_Array,InVar,Weights,NrParameters,AxNum,0);
                TmpShort.Rsquared(Env,1) = NLMOuts.Rsquared.Ordinary(1,1);
                TmpShort.LogLikelihood(Env,1) = NLMOuts.LogLikelihood(1,1);
                TmpShort.Tstat(Env,1) = abs(table2array(NLMOuts.Coefficients(2,3)));
                TmpShort.PValue(Env,1) = table2array(NLMOuts.Coefficients(2,4));
                TmpShort.AIC(Env,1) = (-2*NLMOuts.LogLikelihood(1,1)) - (2 *NrParameters);
            else
                TmpShort.Rsquared(Env,1) = NaN;
                TmpShort.LogLikelihood(Env,1) = NaN;
                TmpShort.Tstat(Env,1) = NaN;
                TmpShort.PValue(Env,1) = NaN;
                TmpShort.AIC(Env,1) = NaN;
            end
        end
        MaxLIkeliHood = max(TmpShort.(genvarname(LookAtStat)));
        MaxLike = find(TmpShort.(genvarname(LookAtStat))== MaxLIkeliHood);
        if length(MaxLike)>1
            Pick = randi(length(MaxLike));
            NrMax = MaxLike(Pick);
            MaxLike = NrMax;
            clear NrMax Pick
        end 
        if isempty(MaxLike) ==1
            StopLike = 1;
        elseif (TmpShort.(genvarname(LookAtStat))(MaxLike) >= (MaxLikeStore-0.01)) || (loop <=5)
            NumbersFixed = [NumbersFixed,MaxLike]; %#ok<*AGROW>
            VariableIteration = NameList(NumbersFixed);
            IterationList.Variables(loop,1:length(NumbersFixed)) = VariableIteration;
            IterationList.Tstat(loop) = TmpShort.Tstat(MaxLike);
            IterationList.LogLikelihood(loop) = TmpShort.LogLikelihood(MaxLike);
            IterationList.Rsquared(loop) = TmpShort.Rsquared(MaxLike);
            IterationList.PValue(loop) = TmpShort.PValue(MaxLike);
            IterationList.AIC(loop) = TmpShort.AIC(MaxLike);
            MaxLikeStore =MaxLIkeliHood; % only here to allow message
            % Remove one loop until it becomes the same
            if NrParameters > NrParameterStore
                Same = 0;
                StoreRemoveVariable = inf;
                Iterations = 0;
            end         
            if  NrParameters >= 3 && Same == 0 && Iterations < 10
                Tmp_Array = ShortendEnvArray(:,NumbersFixed);
                [~,PLS,~] = PLSandNLMRegres(Tmp_Array,InVar,Weights,NrParameters,AxNum,0);
                minAxisWeight = find(PLS.AxisWeightsNormalised(:,1) ==min(PLS.AxisWeightsNormalised(:,1)));
                RemoveVariable = NumbersFixed(minAxisWeight); %#ok<*FNDSB>
                if length(RemoveVariable)>1
                    Pick = randi(length(RemoveVariable));
                    Nrrem = RemoveVariable(Pick);
                    RemoveVariable = Nrrem;
                    clear Nrrem Pick
                end
                NrParameterStore = NrParameters;
                if RemoveVariable == StoreRemoveVariable
                    Same = 1;
                else
                    Same = 0;
                    StoreRemoveVariable = RemoveVariable;
                    NrRemove = find(NumbersFixed == RemoveVariable);
                    NumbersFixed(NrRemove) = [];
                    Iterations = Iterations + 1;
                    NrParameters = NrParameters -1; 
                end
            end
        elseif loop > 50
            StopLike = 1;
        else
            StopLike = 1;
        end
    end % while
    % Save the best one
    NextOne = 0;
    if StartPoint == 1
        NumbersFixedStore = NumbersFixed;
        IterationListStore = IterationList;
        PLSOutStore = PLSOut;
        TMaxStore = MaxLikeStore;
        StartPointStore = 1;
    else
        % Test whether there is a better startingPoint.
        % in case starting points have equal Tmax, take the one with the
        % highets initial value (rounded to avoid equals not being equal
        % problems)
        if round(MaxLikeStore,6) == round(TMaxStore,6)
            WhichOneArray = [abs(SingleFactor.(genvarname(LookAtStat))(StartPointStore)),abs(SingleFactor.(genvarname(LookAtStat))(StartPoint))];
            WhichOneMaxSingle = find(WhichOneArray== (max(WhichOneArray)));
            if  WhichOneMaxSingle == 2 % so the new one was better
                NextOne = 1;
            end
        elseif MaxLikeStore > TMaxStore
            NextOne = 1;
        end
    end
    if NextOne == 1
        NumbersFixedStore = NumbersFixed;
        IterationListStore = IterationList;
        PLSOutStore = PLSOut;
        TMaxStore = MaxLikeStore;
        StartPointStore = StartPoint;
    else
    end
end % StartPoint loop
end % end Forward regression function

%%
function PLS = PLSFinalRun(PLS,EnvironmentalArray,YieldforPLS,NumbersFixed,NameList,Weights,AxNum,Extra)
% Final Run
Final_Array = EnvironmentalArray(:,NumbersFixed);
NrParameters = length(NumbersFixed);
VariablesIncluded = NameList(NumbersFixed);
[NLMOutsFinal,PLSOut,~] = PLSandNLMRegres(Final_Array,YieldforPLS, Weights, NrParameters,AxNum,1);
% for Values
XVarIn = PLSOut.Axis(:,AxNum);
[~, PredictedValues,~] = nlmRegression(XVarIn,YieldforPLS,Weights,{'Dummy'},{'Dummy'},1,2);
[PLSOut.XVarHetero, PLSOut.YVarHetero,~,~] = LinearHeteroskedacityCorrection(XVarIn,YieldforPLS,PredictedValues);
%[NLMOutsFinal,~] = PLSandNLMRegres(XVarHetero,YVarHetero,Weights, 1,AxNum);

if strcmp(Extra,'Extra') ==1
    AxName = {['FixedAxis',num2str(AxNum)]};
    [NLM_Clean,PLS.PerEnv.(genvarname(char(AxName))),EstimatedValues,~] = WhatToOutputPLS(NLMOutsFinal,PLSOut,VariablesIncluded,AxName,AxNum,'Selected',NameList);
else
    AxName = {['CleanedAxis',num2str(AxNum)]};
    AxNameFull = {['CleanedAxisFull',num2str(AxNum)]};
    [NLM_Clean,PLS.PerEnv.(genvarname(char(AxName))),EstimatedValues,PLS.PerEnv.(genvarname(char(AxNameFull)))] = WhatToOutputPLS(NLMOutsFinal,PLSOut,VariablesIncluded,AxName,AxNum,'Selected',NameList);
end

%[OutOverall,OutParameter,EstimatedValues,OutParameterFull] = WhatToOutputPLS(NLM,PLS,Parameters,FullAxName,AxNum,Selected,NameListPLS)
if AxNum == 1 && strcmp(Extra,'Extra') ~=1
    PLS.Values.XSelected = EstimatedValues.X;
    PLS.Values.YSelected = EstimatedValues.Y;
    PLS.Values.XSelectedRunning = EstimatedValues.XRunning;
    PLS.Values.YSelectedRunning = EstimatedValues.YRunning;
    PLS.Values.XSelectedHetero = EstimatedValues.XHetero;
    PLS.Values.YSelectedHetero = EstimatedValues.YHetero;
    PLS.Values.XSelectedHeteroRunning = EstimatedValues.XHeteroRunning;
    PLS.Values.YSelectedHeteroRunning = EstimatedValues.YHeteroRunning;
end

if strcmp(Extra,'Extra') ==1
    PLS.Overview(3,:) = NLM_Clean(1,:);
    PLS.Overview.Parameters(3,1) = {'FixedAxis'};
else
    PLS.Overview(2,:) = NLM_Clean(1,:);
    PLS.Overview.Parameters(2,1) = {'CleanedAxis'};
end

%LEAVEONE OUT PLS
LeaveOutRegressHistoric =  dataset({'dummy'}, 'varnames',{'ParameterRemoved'});
LeaveOutRegressHistoric.ParameterRemoved(1,1) = {'CleanedAxisFull'};
LeaveOutRegressHistoric.PValue(1,1) =  PLS.Overview.PValue(2);
LeaveOutRegressHistoric.Rsquared(1,1) =  PLS.Overview.Rsquared(2);
LeaveOutRegressHistoric.Coefficient(1,1) = PLS.Overview.Coefficient(2);
LeaveOutRegressHistoric.Constant(1,1) = PLS.Overview.Constant(2);
LeaveOutRegressHistoric.Tstat(1,1) = PLS.Overview.Tstat(2);
LeaveOutRegressHistoric.LogLikelihood(1,1) = PLS.Overview.LogLikelihood(2);
LeaveOutRegressHistoric.AIC(1,1) = PLS.Overview.AIC(2);

for ParaNum = 1:length(VariablesIncluded)
    ArrayLeaveOut = Final_Array;
    ArrayLeaveOut(:,ParaNum) = [];
    VariableExcluded = VariablesIncluded(ParaNum);
    NrParametersNew = NrParameters-1;
    [VarLeave,~,~] = PLSandNLMRegres(ArrayLeaveOut,YieldforPLS, Weights, NrParametersNew,AxNum,1);
    LeaveOutRegressHistoric.ParameterRemoved((ParaNum+1),1) = {VariableExcluded};
    LeaveOutRegressHistoric.PValue((ParaNum+1),1) = table2array(VarLeave.Coefficients(2,4));
    LeaveOutRegressHistoric.Rsquared((ParaNum+1),1) = VarLeave.Rsquared.Ordinary(1,1);
    LeaveOutRegressHistoric.Coefficient((ParaNum+1),1) = table2array(VarLeave.Coefficients(2,1));
    LeaveOutRegressHistoric.Constant((ParaNum+1),1) = table2array(VarLeave.Coefficients(1,1));
    LeaveOutRegressHistoric.Tstat((ParaNum+1),1) = table2array(VarLeave.Coefficients(2,3));
    LeaveOutRegressHistoric.LogLikelihood((ParaNum+1),1) = VarLeave.LogLikelihood(1,1);
    LeaveOutRegressHistoric.AIC((ParaNum+1),1) = (-2*VarLeave.LogLikelihood(1,1)) - (2);
end
PLS.(genvarname(['LeaveParametersOutHistoric_',mat2str(AxNum)])) = LeaveOutRegressHistoric;
end
%%
function [XVarOut, YVarOut,XVarOutRunning,YVarOutRunning] = LinearHeteroskedacityCorrection(XVarIn,YVarIn,Predictions)
LengthArray = length(YVarIn);
if LengthArray>=10
    Deviance = YVarIn-Predictions;
    AbsDeviance = abs(Deviance);
    CorrectionFact = 0;
    if min(XVarIn)<0
        CorrectionFact = min(XVarIn);
    end
    XVarCor = XVarIn-CorrectionFact;

    Combined(:,1) = 1:LengthArray;
    Combined(:,2) = XVarCor;
    Combined(:,3) = Predictions;
    Combined(:,4) = Deviance;
    Combined(:,5) = AbsDeviance;
    Sorted = sortrows(Combined,2);
    Count = 1;
    for i = 10:1:LengthArray
        MaxDev(Count,1) = Sorted(i,3); %#ok<*SAGROW> % Predicted Value
        MaxDev(Count,2) = nanmax(Sorted(i-9:i,5)); %#ok<NANMAX> %Absoulte Deviance
        Count = Count + 1;
    end
    HeteromodelFun = @(b,x) b(1)+ b(2).*x;
    HeteroStart = [1; 0];
    HeteroRel = fitnlm(MaxDev(:,1),MaxDev(:,2),HeteromodelFun,HeteroStart);
    HeteroRico = abs(table2array(HeteroRel.Coefficients(2,1))); % between max deviance and predicted
    if HeteroRico < 0
        HeteroRico = -HeteroRico;
    end
    CorOnDeviance = HeteroRico.* (Sorted(:,3)-(Sorted(1,3))); % On Predcited Value

    for i = 1:LengthArray % Deviance Change
        if Sorted(i,4) <0
            Sorted(i,6) = Sorted(i,4) + CorOnDeviance(i);
        else
            Sorted(i,6) = Sorted(i,4) - CorOnDeviance(i);
        end
    end
    % Heterosedacity corrected value
    Sorted(:,7) = Sorted(:,3) + Sorted(:,6); % new Value
    Count = 1;
    for i = 10:1:LengthArray
        XVarOutRunning(Count,1) = nanmean(Sorted(i,2)); %#ok<*AGROW,*SAGROW> % X value
        YVarOutRunning(Count,1) = nanmean(Sorted(i-9:i,7)); %Running mean new Values
        Count = Count + 1;
    end
    ReSorted = sortrows(Sorted,1);
    XVarOut = ReSorted (:,2);
    XVarOut = XVarOut + CorrectionFact;
    YVarOut = ReSorted (:,7);
else
    XVarOut = XVarIn;
    YVarOut = YVarIn;
    XVarOutRunning = XVarIn;
    YVarOutRunning = YVarIn;
end
end

%%
function InVar = MeanPerYear(EnvIn,~,Name,InVar,YearsIn)
YearStart = 1970;
YearEnd = 2022;
for yearrun =  YearStart:YearEnd
    List = find(YearsIn == yearrun);
    if isempty(List)~=1
        EnvPerYear = nanmean(EnvIn(List));
        STDPerYear = nanstd(EnvIn(List));
        NPerYear = length(EnvIn(List));
    else
        EnvPerYear = NaN;
        STDPerYear = NaN;
        YieldPerYear = NaN;
        NPerYear = NaN;
    end
    WhichYear = find(InVar.Median.Year == yearrun);
    InVar.Median.(genvarname(char(Name)))(WhichYear,1) = EnvPerYear;
    InVar.STD.(genvarname(char(Name)))(WhichYear,1) = STDPerYear;
    InVar.NPerYear.(genvarname(char(Name)))(WhichYear,1) = NPerYear;
    clear EnvPerYear STDPerYear NPerYear
end
end
%%
function InVar = MeanPerYearYield(Yield,InVar,YearsIn,Name)
YearStart = 1970;
YearEnd = 2022;
for yearrun =  YearStart:YearEnd
    List = find(YearsIn == yearrun);
    if isempty(List)~=1
        PerYear = nanmean(Yield(List));
        StdPerYear =  nanstd(Yield(List));
        NPerYear = length(Yield(List));
    else
        PerYear = NaN;
        StdPerYear =  NaN;
        NPerYear = NaN;
    end
    WhichYear = find(InVar.Median.Year == yearrun);
    InVar.Median.(genvarname(char(Name)))(WhichYear,1) = PerYear;
    InVar.STD.(genvarname(char(Name)))(WhichYear,1) = StdPerYear;
    InVar.NPerYear.(genvarname(char(Name)))(WhichYear,1) = NPerYear;
    clear PerYear StdPerYear NPerYear
end
end
%%
function InVar = MeanPerYearFuture(EnvIn,Name,InVar,LTEs)
% FutureYears.(genvarname(char(OnetoFind))).Environment = MeanPerYearFuture(FutEnviron,Name,FutureYears.(genvarname(char(OnetoFind))).Environment,LTEsType);
UniList = unique(LTEs);
% Make one value per LTEs (avoiding unweighted averaging)
Start = find(strcmp(EnvIn.Properties.VariableNames,'Y_2023')==1);
for i = 1:length(UniList)
    Tester = find(strcmp(LTEs,UniList(i))==1);
    EnviromentSingle(i,:) = table2array(EnvIn(Tester(1),Start:end));
end
if size(EnviromentSingle,1) == 1
    EnviromentSingle(2,:) = NaN;
end
EnvPerYear = nanmean(EnviromentSingle,1);
STDPerYear = nanstd(EnviromentSingle,1);
InVar.Median.(genvarname(char(Name)))(:,1) = EnvPerYear';
InVar.STD.(genvarname(char(Name)))(:,1) = STDPerYear';
end

%%
function PLS = YieldDevelopment(PLS,PerYear,FutureYears,AxNum,~,ActualYield)
AxName = ['CleanedAxis',num2str(AxNum)];

Yield = double(ActualYield(:,1));

Transfer.Betas = PLS.PerEnv.(genvarname(AxName)).Betas;
Transfer.RelativeWeights =PLS.PerEnv.(genvarname(AxName)).RelativeWeights;
Transfer.ParametersNames = PLS.PerEnv.(genvarname(AxName)).Parameters;
ConstantNum = find(strcmp(char(Transfer.ParametersNames),"Constant"));
Transfer.Constant = double(PLS.PerEnv.(genvarname(AxName)).Betas(ConstantNum));

ConstantNum = find(strcmp(char(Transfer.ParametersNames),"Constant"));
Years = double(PerYear.Median.Year);
RangeYears = (max(Years) -  min(Years))+1;
Transfer.NumberParameters = ConstantNum-1;

YieldDevelopments = zeros(RangeYears,3);
Zvalue = 1; % %1.645; % Corresponding to 10% and 90% Confidence Interval
NMax = max(PerYear.NPerYear.(genvarname(char(Transfer.ParametersNames(1)))));
% Historic predictions
for yearstodo = 1:RangeYears
    YieldDevelopments(yearstodo,1) = Years(yearstodo);
    Value = Transfer.Constant;
    ValueSTD = 0;
    for Paranum = 1:Transfer.NumberParameters
        YearNum = find(Years == Years(yearstodo));
        STDFact = 1./(sqrt((PerYear.NPerYear.(genvarname(char(Transfer.ParametersNames(Paranum))))(YearNum))./NMax));
        Value = Value + double(Transfer.Betas(Paranum)).* PerYear.Median.(genvarname(char(Transfer.ParametersNames(Paranum))))(YearNum);
        ValueSTD = ValueSTD + (((PerYear.STD.(genvarname(char(Transfer.ParametersNames(Paranum))))(YearNum))).*STDFact).*double(Transfer.RelativeWeights(Paranum));
    end
    YieldDevelopments(yearstodo,2) = Value;
    if ValueSTD == 0 && yearstodo ~=1
        ValueSTD = nanmean(YieldDevelopments(1:yearstodo,3));
    end
    YieldDevelopments(yearstodo,3) = ValueSTD;
end

%Scale factors for the differences to 1, so above will be scaled to above;
% below to below
Transfer.MaxFactor = ((prctile(Yield,95))-1)./((prctile(YieldDevelopments(:,2),95))-1);
Transfer.MinFactor = ((prctile(Yield,5))-1)./((prctile(YieldDevelopments(:,2),5))-1);

FutureYearsArray = double(FutureYears.Median.Year);
FutureRangeYears = (max(FutureYearsArray) -  min(FutureYearsArray))+1;
Transfer.StartFutureYears = yearstodo + 1;
Transfer.EndFutureYears = yearstodo + FutureRangeYears;
Count = 0;

%Correction = nanmean(YieldDevelopments(:,2)) - Transfer.Constant;
Correction =  nanmean(Yield) - Transfer.Constant;
% Build in step for actual values

CorNum = 4;
IsnanListLastFive = find(isnan(Yield((end-CorNum):end)));
if isempty(IsnanListLastFive) == 1
    Tester = Yield((end-CorNum):end);
else
    NonNan = find(isfinite(Yield));
    Tester = Yield(NonNan(end-CorNum:end));
end
FutureConstant = nanmedian(Tester)-Correction;

% Future Predictions
for yearstodoFut = Transfer.StartFutureYears:Transfer.EndFutureYears
    Count = Count + 1;
    YieldDevelopments(yearstodoFut,1) = FutureYearsArray(Count);
    Value = FutureConstant;
    ValueSTD = 0;
    for Paranum = 1:Transfer.NumberParameters
        YearNum = find(FutureYearsArray == FutureYearsArray(Count));
        Value =    Value + double(Transfer.Betas(Paranum)).* FutureYears.Median.(genvarname(char(Transfer.ParametersNames(Paranum))))(YearNum);
        ValueSTD = ValueSTD + (FutureYears.STD.(genvarname(char(Transfer.ParametersNames(Paranum))))(YearNum)).*double(Transfer.RelativeWeights(Paranum)) ;
    end
    YieldDevelopments(yearstodoFut,2) = Value;
    YieldDevelopments(yearstodoFut,3) = ValueSTD;
    YieldDevelopments(1,4) = FutureConstant;
    YieldDevelopmentFuture(Count,1) = FutureYearsArray(Count);
    YieldDevelopmentFuture(Count,2) = Value;
end
% ReScale Development around 1
for i = 1:length(YieldDevelopments(:,2))
    if YieldDevelopments(i,2) < 1
        YieldDevelopments(i,2) = ((YieldDevelopments(i,2)-1).*Transfer.MinFactor)+1;
    else
        YieldDevelopments(i,2) = ((YieldDevelopments(i,2)-1).*Transfer.MaxFactor)+1;
    end
end
YieldDevelopment = dataset(YieldDevelopments(:,1),'Varnames','Years');
YieldDevelopment.Prediction = YieldDevelopments(:,2);
YieldDevelopment.STD = YieldDevelopments(:,3);
YieldDevelopment.MinRange = YieldDevelopments(:,2) - (Zvalue.*YieldDevelopments(:,3));
YieldDevelopment.MaxRange = YieldDevelopments(:,2) + (Zvalue.*YieldDevelopments(:,3));
YieldDevelopment.FutureConstant = YieldDevelopments(:,4);
PLS.(genvarname(['YieldDevelopment_Axis_',mat2str(AxNum)])) = YieldDevelopment;
end

%%
function PLSIn = YieldDevelopmentCombine(PLSIn,~,PerYear,NLMNeed,ActualYield)
NLM = NLMNeed;
AxName1 = 'YieldDevelopment_Axis_1';
AxName2 = 'YieldDevelopment_Axis_2';
Yield = double(ActualYield(:,1));
Overview = PLSIn.Overview(2,:);

Years = double(PerYear.Median.Year);
LengthYears =size(Years,1);
LengthYearsStart = (LengthYears+1);
RangeYears = (max(Years) -  min(Years))+1;
CorrectionCurrent = nanmedian(Yield((end-4):end,1));

NewAx = PLSIn.(genvarname(AxName1))(:,1:5);
Ax1 =  double(PLSIn.(genvarname(AxName1))(:,2:5)).*double(Overview(1,5));
Ax2 =  double(PLSIn.(genvarname(AxName2))(:,2:5)).*double(Overview(1,7));
New = (Ax1 + Ax2)./((double(Overview(1,5)))+(double(Overview(1,7))));
NewAx.Prediction = New(:,1);
NewAx.STD = New(:,2);
NewAx.MinRange = New(:,3);
NewAx.MaxRange = New(:,4);
CorrectionFuture =nanmean(New(LengthYearsStart:(LengthYearsStart+2),1));
Correction = CorrectionFuture-CorrectionCurrent;

NewAx.YieldPlusPrediction = NewAx.Prediction-Correction;
NewAx.YieldPlusPrediction(1:RangeYears,1) = ActualYield;
NewAx.MinRangeYplusP = NewAx.MinRange - Correction;
NewAx.MaxRangeYplusP = NewAx.MaxRange - Correction;
NewAx.MinRangeYplusP(1:RangeYears,1) = NewAx.YieldPlusPrediction(1:RangeYears,1) - NewAx.STD(1:RangeYears,1);
NewAx.MaxRangeYplusP(1:RangeYears,1) = NewAx.YieldPlusPrediction(1:RangeYears,1) + NewAx.STD(1:RangeYears,1);
PLSIn.YieldDevelopmentTwoAxes = NewAx;


% Make combi relative weights
PLSIn.PerEnv.CombiCleanedFull = PLSIn.PerEnv.CleanedAxisFull1;
Correct = 1./(((double(Overview(1,5))) + (double(Overview(1,7)))));

Summer = (double(Overview(1,5)).*(double(PLSIn.PerEnv.CleanedAxisFull1(:,2)))) + (double(Overview(1,7)).*(double(PLSIn.PerEnv.CleanedAxisFull2(:,2))));
Summer = Summer.*Correct;
PLSIn.PerEnv.CombiCleanedFull.Weightfactors = Summer;
Summer =  (double(Overview(1,5)).*(double(PLSIn.PerEnv.CleanedAxisFull1(:,3)))) + (double(Overview(1,7)).*(double(PLSIn.PerEnv.CleanedAxisFull2(:,3))));
Summer = Summer.*Correct;
Summer = Summer./sum(Summer);
PLSIn.PerEnv.CombiCleanedFull.RelativeWeights = Summer;
Summer =  (double(Overview(1,5)).*(double(PLSIn.PerEnv.CleanedAxisFull1(:,4)))) + (double(Overview(1,7)).*(double(PLSIn.PerEnv.CleanedAxisFull2(:,4))));
Summer = Summer.*Correct;
PLSIn.PerEnv.CombiCleanedFull.Betas = Summer;
Len = length(PLSIn.PerEnv.CombiCleanedFull.Parameters)-1;
PLSIn.PerEnv.CombiCleanedFull.NLMCoefficients = NLMNeed.Coefficient(1:Len,1);
end

%%
function CorrelationMatrix = CorrelationMatrixGeneration(VarIn,NameList)
CorrelationMatrix.Rho = dataset({'Dummy'},'VarNames','Parameter');
CorrelationMatrix.PValue = dataset({'Dummy'},'VarNames','Parameter');
for par = 1:length(NameList)
    for par2 = 1:length(NameList)
        CorrelationMatrix.Rho.Parameter(par,1) = NameList(par);
        CorrelationMatrix.PValue.Parameter(par,1) = NameList(par);
        CorrelationMatrix.Rho.N(par,1) = length(VarIn(:,par));
        CorrelationMatrix.PValue.N(par,1) = length(VarIn(:,par));
        [RMat,PMat] = corrcoef(VarIn(:,par),VarIn(:,par2));
        CorrelationMatrix.Rho.(genvarname(char(NameList(par))))(par2,1) = RMat(1,2);
        CorrelationMatrix.PValue.(genvarname(char(NameList(par))))(par2,1) = PMat(1,2);
    end
end
end
%%
function ChangeType = RainCWDCutOffFunc(OnetoFind,InPutArray,NumberofObservations,TypesYield,YieldToRun,RainCWDCutOffs,Input)
if strcmp(Input,'Run it with preset cut-offs') == 1
    RunIt = 1;
end
OnetoFind = char(OnetoFind);
Typea = OnetoFind(1:2);
Typeb = OnetoFind(3:4);
Test1 = strcmp(Typea,'Ra');
Test2 = strcmp(Typeb,'Ra');
Test = Test1 + Test2;
if Test > 0
    InPuts = double(InPutArray.RainFallChange);
    RainName = [{'RainDecrease'},{'RainIncrease'}];
else
    InPuts = double(InPutArray.CWDChange);
    RainName = [{'CWDDecrease'},{'CWDIncrease'}];
end
if RunIt == 1
    CutOff = double(RainCWDCutOffs.(genvarname(char(TypesYield(YieldToRun)))).(genvarname(char(OnetoFind))));
else
    CutOff = 0;
end

RainFallChange = InPuts;
ChangeType(1:NumberofObservations,1) = RainName(1);
ListChange = find(RainFallChange >CutOff);
ChangeType(ListChange,1) = RainName(2);
end


%% 
function StarSignificant = HoghbergFunc(VarIn) %Hoghberg
Nr = find((isnan(VarIn))==1);
[h(:,1)]=Hoghberg(VarIn,0.05,'pdep','no');
[h(:,2)]=Hoghberg(VarIn,0.01,'pdep','no');
[h(:,3)]=Hoghberg(VarIn,0.001,'pdep','no');
StarSignificant = sum(h,2);
StarSignificant(Nr) = NaN;
end
%%
function VarOut = ExtractNLMResults(VarIn)
NameList = [{'DailyTemperature'},'MinTemp','MaxTemp',{'DailyTempRange'},'AnnualTempRange',...
            'Isothermality','CO2Norm','PET','ETOSeasonality','Rainfall','PrecipitationSeasonality',...
            'CWD','Radiation','PopPressure','Ozone','DailyMeanTempChange','Years','YearsNotVarityCorrected']; % 17

RunningPars.Types = fieldnames(VarIn.HeteroCorrected); 
% RunningPars.Types = [{'AllData'};{'C3'};{'C4'};... 
%         'Temperate';'Tropical';...
%         'C3Temperate';'C3Tropical';'C4Temperate';'C4Tropical';... 
%         'CWDDecrease';'CWDIncrease'; ... 
%         'C3CWDDecrease';'C3CWDIncrease';'C4CWDDecrease';'C4CWDIncrease';...
%         'EuropeGlobal';'SubSaharanAfrica';'USMexico'; 'OtherParts';... 
%         'Maize';'Wheat';'Oats';'Barley';'Rice';... 
%         'MaizeTropical'; 'MaizeTemperate';'WheatTropical'; 'WheatTemperate';... 
%         'MaizeAfrica';'MaizeEurope';'MaizeNorthAmericas';'MaizeSouthAmericas';... 
%         'WheatEurope';'WheatNorthAmericas';'WheatSouthAmericas';... 
%        'C3SubSaharanAfrica'; 'C3EuropeGlobal';'C3USMexico';'C3OtherParts';... 
%         'C4SubSaharanAfrica'; 'C4EuropeGlobal';'C4USMexico';'C4OtherParts';... 
%         'RainDecrease';'RainIncrease';'C3RainDecrease';'C3RainIncrease';... 
%         'C4RainDecrease';'C4RainIncrease']; % 4 (49) = 44:49

Testy = find((strcmp(VarIn.HeteroCorrected.(genvarname(char(RunningPars.Types(1)))).Parameters,'YearsNotVarityCorrected'))==1);
if isempty(Testy)== 1
    ToDoList = length(NameList)-2;
else
    ToDoList = length(NameList);
end

for i = 1:ToDoList
    if isempty(Testy)== 1
        Test = VarIn.HeteroCorrectedRunning.(genvarname(char(RunningPars.Types(1)))).Parameters;
        NameListNumber = find((strcmp((Test),(NameList(i))))==1);
    else
        if i == (length(NameList)-1)
            NameListNumber = 20;
        elseif i == length(NameList)
            NameListNumber = 21;
        else
            Test = VarIn.HeteroCorrectedRunning.(genvarname(char(RunningPars.Types(1)))).Parameters;
            NameListNumber = find((strcmp((Test),(NameList(i))))==1);
        end
    end

    List = dataset({'Dummy'},'VarNames','Parameter');
    Tester = fieldnames(VarIn.HeteroCorrectedRunning);
    for j = 1:length(RunningPars.Types)
       if sum(strcmp(Tester,RunningPars.Types(j))) > 0
       List.Parameter(j,1) =  {NameList(i)};
       List.Region(j,1) =  {RunningPars.Types(j)};
       List.Coefficient(j,1) = {VarIn.HeteroCorrectedRunning.(genvarname(char(RunningPars.Types(j)))).Coefficient(NameListNumber)};
       List.HochBergRestricted(j,1) = {VarIn.HeteroCorrected.(genvarname(char(RunningPars.Types(j)))).HochbergRestricted(NameListNumber)};
       List.Constant(j,1) = {VarIn.HeteroCorrectedRunning.(genvarname(char(RunningPars.Types(j)))).Constant(NameListNumber)};
       List.TStat(j,1) = {VarIn.HeteroCorrected.(genvarname(char(RunningPars.Types(j)))).Tstat(NameListNumber)};
       List.PValue(j,1) = {VarIn.HeteroCorrected.(genvarname(char(RunningPars.Types(j)))).PValue(NameListNumber)};
       List.RSquare(j,1) = {VarIn.HeteroCorrectedRunning.(genvarname(char(RunningPars.Types(j)))).Rsquared(NameListNumber)};
       List.N(j,1) = {VarIn.HeteroCorrected.(genvarname(char(RunningPars.Types(j)))).N(NameListNumber)};
       List.HochBerg(j,1) = {VarIn.HeteroCorrected.(genvarname(char(RunningPars.Types(j)))).Hochberg(NameListNumber)};
       end
    end
    VarOut.(genvarname(char(NameList(i)))) = List;
    clear List
end
end
%%
function StructOut = ExtractPLSResults(StructIn)
warning off
RunningPars.Types = fieldnames(StructIn);


for j = 1:length(RunningPars.Types)
    Twentyeleven = find(StructIn.(genvarname(char(RunningPars.Types(j)))).YieldDevelopmentTwoAxes.Years == 2011);
    NowAvg = nanmean(StructIn.(genvarname(char(RunningPars.Types(j)))).YieldDevelopmentTwoAxes.YieldPlusPrediction(Twentyeleven:(Twentyeleven+9),1));
    Twentyninity = find(StructIn.(genvarname(char(RunningPars.Types(j)))).YieldDevelopmentTwoAxes.Years == 2090);
    ThenAvg = nanmean(StructIn.(genvarname(char(RunningPars.Types(j)))).YieldDevelopmentTwoAxes.YieldPlusPrediction(Twentyninity:(Twentyninity+9),1));
    PropChange = ThenAvg;%./NowAvg;

    ThenAvg = nanmean(StructIn.(genvarname(char(RunningPars.Types(j)))).YieldDevelopmentTwoAxes.MinRangeYplusP(Twentyninity:(Twentyninity+9),1));
    PropChangeMin = ThenAvg;%./NowAvg;

    ThenAvg = nanmean(StructIn.(genvarname(char(RunningPars.Types(j)))).YieldDevelopmentTwoAxes.MaxRangeYplusP(Twentyninity:(Twentyninity+9),1));
    PropChangeMax = ThenAvg;%./NowAvg;

    StructIn.(genvarname(char(RunningPars.Types(j)))).Overview.PropChange(2,1) = PropChange;
    StructIn.(genvarname(char(RunningPars.Types(j)))).Overview.PropChangeMinRange(2,1) = PropChangeMin;
    StructIn.(genvarname(char(RunningPars.Types(j)))).Overview.PropChangeMaxRange(2,1) = PropChangeMax;

    StructOut.ResultsOverview(j,:) = StructIn.(genvarname(char(RunningPars.Types(j)))).Overview(2,:);
    StructOut.ResultsOverview.Parameters{j,1} = RunningPars.Types(j);
    StructOut.(genvarname(char(RunningPars.Types(j)))).CombiCleanedFull = StructIn.(genvarname(char(RunningPars.Types(j)))).PerEnv.CombiCleanedFull;
    StructOut.(genvarname(char(RunningPars.Types(j)))).CombiCleanedFull = sortrows(StructOut.(genvarname(char(RunningPars.Types(j)))).CombiCleanedFull,3,'descend');
    StructOut.(genvarname(char(RunningPars.Types(j)))).YieldDevelopmentTwoAxes = StructIn.(genvarname(char(RunningPars.Types(j)))).YieldDevelopmentTwoAxes;
    StructOut.(genvarname(char(RunningPars.Types(j)))).YieldDevelopmentTwoAxes = StructOut.(genvarname(char(RunningPars.Types(j)))).YieldDevelopmentTwoAxes(:,[1:2 4:6 3 7:end]);
end
end
%%
%% 
% NOT MINE, via https://nl.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh

% fdr_bh() - Executes the Benjamini & Hochberg (1995) and the Benjamini &
%            Yekutieli (2001) procedure for controlling the false discovery
%            rate (FDR) of a family of hypothesis tests. FDR is the expected
%            proportion of rejected hypotheses that are mistakenly rejected
%            (i.e., the null hypothesis is actually true for those tests).
%            FDR is a somewhat less conservative/more powerful method for
%            correcting for multiple comparisons than procedures like Bonferroni
%            correction that provide strong control of the family-wise
%            error rate (i.e., the probability that one or more null
%            hypotheses are mistakenly rejected).
%
%            This function also returns the false coverage-statement rate
%            (FCR)-adjusted selected confidence interval coverage (i.e.,
%            the coverage needed to construct multiple comparison corrected
%            confidence intervals that correspond to the FDR-adjusted p-values).
%
%
% Usage:
%  >> [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals,q,method,report);
%
% Required Input:
%   pvals - A vector or matrix (two dimensions or more) containing the
%           p-value of each individual test in a family of tests.
%
% Optional Inputs:
%   q       - The desired false discovery rate. {default: 0.05}
%   method  - ['pdep' or 'dep'] If 'pdep,' the original Bejnamini & Hochberg
%             FDR procedure is used, which is guaranteed to be accurate if
%             the individual tests are independent or positively dependent
%             (e.g., Gaussian variables that are positively correlated or
%             independent).  If 'dep,' the FDR procedure
%             described in Benjamini & Yekutieli (2001) that is guaranteed
%             to be accurate for any test dependency structure (e.g.,
%             Gaussian variables with any covariance matrix) is used. 'dep'
%             is always appropriate to use but is less powerful than 'pdep.'
%             {default: 'pdep'}
%   report  - ['yes' or 'no'] If 'yes', a brief summary of FDR results are
%             output to the MATLAB command line {default: 'no'}
%
%
% Outputs:
%   h       - A binary vector or matrix of the same size as the input "pvals."
%             If the ith element of h is 1, then the test that produced the
%             ith p-value in pvals is significant (i.e., the null hypothesis
%             of the test is rejected).
%   crit_p  - All uncorrected p-values less than or equal to crit_p are
%             significant (i.e., their null hypotheses are rejected).  If
%             no p-values are significant, crit_p=0.
%   adj_ci_cvrg - The FCR-adjusted BH- or BY-selected
%             confidence interval coverage. For any p-values that
%             are significant after FDR adjustment, this gives you the
%             proportion of coverage (e.g., 0.99) you should use when generating
%             confidence intervals for those parameters. In other words,
%             this allows you to correct your confidence intervals for
%             multiple comparisons. You can NOT obtain confidence intervals
%             for non-significant p-values. The adjusted confidence intervals
%             guarantee that the expected FCR is less than or equal to q
%             if using the appropriate FDR control algorithm for the
%             dependency structure of your data (Benjamini & Yekutieli, 2005).
%             FCR (i.e., false coverage-statement rate) is the proportion
%             of confidence intervals you construct
%             that miss the true value of the parameter. adj_ci=NaN if no
%             p-values are significant after adjustment.
%   adj_p   - All adjusted p-values less than or equal to q are significant
%             (i.e., their null hypotheses are rejected). Note, adjusted
%             p-values can be greater than 1.
%
%
% References:
%   Benjamini, Y. & Hochberg, Y. (1995) Controlling the false discovery
%     rate: A practical and powerful approach to multiple testing. Journal
%     of the Royal Statistical Society, Series B (Methodological). 57(1),
%     289-300.
%
%   Benjamini, Y. & Yekutieli, D. (2001) The control of the false discovery
%     rate in multiple testing under dependency. The Annals of Statistics.
%     29(4), 1165-1188.
%
%   Benjamini, Y., & Yekutieli, D. (2005). False discovery rate?adjusted
%     multiple confidence intervals for selected parameters. Journal of the
%     American Statistical Association, 100(469), 71?81. doi:10.1198/016214504000001907
%
%
% Example:
%  nullVars=randn(12,15);
%  [~, p_null]=ttest(nullVars); %15 tests where the null hypothesis
%  %is true
%  effectVars=randn(12,5)+1;
%  [~, p_effect]=ttest(effectVars); %5 tests where the null
%  %hypothesis is false
%  [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh([p_null p_effect],.05,'pdep','yes');
%  data=[nullVars effectVars];
%  fcr_adj_cis=NaN*zeros(2,20); %initialize confidence interval bounds to NaN
%  if ~isnan(adj_ci_cvrg),
%     sigIds=find(h);
%     fcr_adj_cis(:,sigIds)=tCIs(data(:,sigIds),adj_ci_cvrg); % tCIs.m is available on the
%     %Mathworks File Exchagne
%  end
%
%
% For a review of false discovery rate control and other contemporary
% techniques for correcting for multiple comparisons see:
%
%   Groppe, D.M., Urbach, T.P., & Kutas, M. (2011) Mass univariate analysis
% of event-related brain potentials/fields I: A critical tutorial review.
% Psychophysiology, 48(12) pp. 1711-1725, DOI: 10.1111/j.1469-8986.2011.01273.x
% http://www.cogsci.ucsd.edu/~dgroppe/PUBLICATIONS/mass_uni_preprint1.pdf
%
%
% For a review of FCR-adjusted confidence intervals (CIs) and other techniques
% for adjusting CIs for multiple comparisons see:
%
%   Groppe, D.M. (in press) Combating the scientific decline effect with
% confidence (intervals). Psychophysiology.
% http://biorxiv.org/content/biorxiv/early/2015/12/10/034074.full.pdf
%
%
% Author:
% David M. Groppe
% Kutaslab
% Dept. of Cognitive Science
% University of California, San Diego
% March 24, 2010
%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
%
% 5/7/2010-Added FDR adjusted p-values
% 5/14/2013- D.H.J. Poot, Erasmus MC, improved run-time complexity
% 10/2015- Now returns FCR adjusted confidence intervals
function [h, crit_p, adj_ci_cvrg, adj_p]=Hoghberg(pvals,q,method,report)
if nargin<1
    error('You need to provide a vector or matrix of p-values.');
else
    if ~isempty(find(pvals<0,1))
        error('Some p-values are less than 0.');
    elseif ~isempty(find(pvals>1,1))
        error('Some p-values are greater than 1.');
    end
end
if nargin<2
    q=.05;
end
if nargin<3
    method='pdep';
end
if nargin<4
    report='no';
end
s=size(pvals);
if (length(s)>2) || s(1)>1
    [p_sorted, sort_ids]=sort(reshape(pvals,1,prod(s)));
else
    %p-values are already a row vector
    [p_sorted, sort_ids]=sort(pvals);
end
[~, unsort_ids]=sort(sort_ids); %indexes to return p_sorted to pvals order
m=length(p_sorted); %number of tests
if strcmpi(method,'pdep')
    %BH procedure for independence or positive dependence
    % see also https://www.statisticshowto.com/benjamini-hochberg-procedure/
    thresh=(1:m)*q/m;
    wtd_p=m*p_sorted./(1:m);
elseif strcmpi(method,'dep')
    %BH procedure for any dependency structure
    denom = m*sum(1./(1:m)) ;
    thresh=(1:m)*q/denom;
    wtd_p=denom*p_sorted./[1:m]; %#ok<*NBRAK1>
    %Note, it can produce adjusted p-values greater than 1!
    %compute adjusted p-values
else
    error('Argument ''method'' needs to be ''pdep'' or ''dep''.');
end
if nargout>3
    %compute adjusted p-values; This can be a bit computationally intensive
    adj_p=zeros(1,m)*NaN;
    [wtd_p_sorted, wtd_p_sindex] = sort( wtd_p );
    nextfill = 1;
    for k = 1 : m
        if wtd_p_sindex(k)>=nextfill
            adj_p(nextfill:wtd_p_sindex(k)) = wtd_p_sorted(k);
            nextfill = wtd_p_sindex(k)+1;
            if nextfill>m
                break;
            end
        end
    end
    adj_p=reshape(adj_p(unsort_ids),s);
end
rej=p_sorted<=thresh;
max_id=find(rej,1,'last'); %find greatest significant pvalue
if isempty(max_id)
    crit_p=0;
    h=pvals*0;
    adj_ci_cvrg=NaN;
else
    crit_p=p_sorted(max_id);
    h=pvals<=crit_p;
    adj_ci_cvrg=1-thresh(max_id);
end
if strcmpi(report,'yes')
    n_sig=sum(p_sorted<=crit_p);
    if n_sig==1
        fprintf('Out of %d tests, %d is significant using a false discovery rate of %f.\n',m,n_sig,q);
    else
        fprintf('Out of %d tests, %d are significant using a false discovery rate of %f.\n',m,n_sig,q);
    end
    if strcmpi(method,'pdep')
        fprintf('FDR/FCR procedure used is guaranteed valid for independent or positively dependent tests.\n');
    else
        fprintf('FDR/FCR procedure used is guaranteed valid for independent or dependent tests.\n');
    end
end
end
%%
% See April version for all redundant codes