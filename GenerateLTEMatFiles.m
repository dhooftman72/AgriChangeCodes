% Extract Data from preformatted csv files for AgriChange
% Calcuates the Variety corrected yields, as well as the median corrected
% yield (for time series)
% Calcuates/extracts the best performing and restricted treatments
% Includes testing for common errors
% Starts in current and ends in current directory.

% Function include
% Main function: GenerateLTEMatFiles()
% Daughters:
% - ExtractLTEData(File, Current, FileDir, VarMis, MisVarieties, LimitedCount)
% - ExtractYieldData(File, LimitedCount, Current)
% - CalculateYearMeans(InputTreatmentYield, ListerYears, File)
% - SetAgainstMedian(Yields)
% - ReadjustVarieties(Yields, ListerYears)

function [MisVarieties,LTEFullSelected] = GenerateLTEMatFiles
Current = pwd;
FileDir = 'C:\DannyData\Projects\Agrichange\Yield csv files';
cd(FileDir)
List = dir;
length_names = length(List);
MisVarieties = {'Dummy'};
% Renove non csv files
Count = 1;
for x = 1:1:length_names
    name_list_temp = cellstr(List(x,1).name);
    if strcmp(name_list_temp, '.') == 1 || strcmp(name_list_temp, '..') == 1
        ToRemoveList(Count) = x;
        Count = Count + 1;
    else
        NameChar = char(name_list_temp);
        Last = NameChar((end-2):end);
        if strcmp(Last, 'csv') ~= 1
            ToRemoveList(Count) = x; %#ok<*AGROW>
            Count = Count + 1;
        end
    end
end
List(ToRemoveList) = [];

length_names = length(List);
fprintf('Looping through %i files \n', length_names)
% extrac the data from the standard format
VarMis = 1;
LimitedCount.NPKNr = 0;
LimitedCount.MonoNr =  0;
LimitedCount.NPK = {'Dummy'};
LimitedCount.Mono = {'Dummy'};

for x = 1:1:length_names
    if x/25 == round(x/25)
       Message = ['# LTEs done: ', mat2str(x)]; 
       disp(Message)
    end
    %Message = ['Running LTE: ',cellstr(List(x,1).name)];
    %disp(Message)
    [MisVarieties,VarMis,LTETesterSelect,LimitedCount] = ExtractLTEData(char(cellstr(List(x,1).name)),Current,FileDir,VarMis,MisVarieties,LimitedCount); 
    if x == 1
        LTEFullSelected = LTETesterSelect;
    else
        LTEFullSelected = [LTEFullSelected;LTETesterSelect];
    end
    clear LTETesterSelect 
end
 cd(Current)
 save('InputData.mat','LimitedCount', '-append');
end

% The extraction Function
function [MisVarieties,VarMis,LTETesterSelect,LimitedCount] = ExtractLTEData(File,Current,FileDir,VarMis,MisVarieties,LimitedCount)
warning off
CheckVars =  [{'Order'},'TreatType','Replica','Year','Yield','Nutrients_Manure','Residue',...
    'Irrigation','Ctype','Variety','Exclude','Rotation','Crop','SubReplica'];
FileChar = char(File);
FileShort = FileChar(1:3);
FirstLetter = FileChar(1);
FileSecond = FileChar(4:end-4);
Array = readtable(File);
ArrayDataSet = table2dataset(Array);
Tester = ArrayDataSet.Properties.VarNames;

% Set subreplica if not present
SubRepPres = sum(strcmp(Tester,"SubReplica"));
if SubRepPres == 0
    Columns = 13;
else
    Columns = 14;
end
ArrayDataSet(:,(Columns+1):end) = [];
if SubRepPres == 0
    ArrayDataSet.SubReplica(:,1) = 1;
end

% Test for correct variable names
Tester = ArrayDataSet.Properties.VarNames;
VarCheck = sum(strcmp(CheckVars,Tester));
if VarCheck ~= 14
   Message = ['The variable names are not correct for file: ',FileChar];
   disp(Message)
   ccc
end

% Test for correct crop C-type combis' (the few others are manually
% checked)
Crop = char(ArrayDataSet.Crop(1));
Ctype = char(ArrayDataSet.Ctype(1));
if strcmp(Crop, 'Wheat') == 1 || strcmp(Crop, 'Barley') == 1 || strcmp(Crop, 'Oats') == 1 || strcmp(Crop, 'Soybean') == 1
    if strcmp(Ctype, 'C3') ~= 1
        Message = ['C-type is not correct for file: ',FileChar];
        disp(Message)
        ccc
    end
elseif strcmp(Crop, 'Maize') == 1 || strcmp(Crop, 'Millet') == 1 || strcmp(Crop, 'Sorghum') == 1
    if strcmp(Ctype, 'C4') ~= 1
        Message = ['C-type is not correct for file: ',FileChar];
        disp(Message)
        ccc
    end
end

% Test for wrong Yields (kgs instead of tons per ha)
Yield = char(ArrayDataSet.Yield(:));
TestVar = find(Yield>120); %#ok<EFIND> % strange number due to PSP sugarbeet
FirstMention = 1;
if isempty(TestVar)~=1 && FirstMention == 1
   Message = ['Wrong Calcuation of Yields (kg per ha): ',FileChar];
   disp(Message)
end


% Test for unknown Varieties and change to Name
Variety = char(ArrayDataSet.Variety(:));
TestVar = sum(strcmp(Variety,{'Unknown'}));
TestVar2 = sum(strcmp(Variety,{'-'}));
TestVar = TestVar + TestVar2;
if TestVar > 0
   Message = ['Varietes recorded as unknown: ',FileChar];
   if strcmp(FirstLetter,'D')~=1 % Do nor show when DRIVES
   disp(Message)
   end
   MisVarieties(VarMis,1) = {FileChar};
   VarMis = VarMis + 1;
   ListUnknown(:,1) = find((strcmp(Variety,{'Unknown'}))==1);
   ListUnknown2(:,1) = find((strcmp(Variety,{'-'}))==1);
   ListUnknown = [ListUnknown;ListUnknown2];
   for i = 1:length(ListUnknown)
   ArrayDataSet.Variety(ListUnknown(i),1) = {['Variety_', mat2str(ArrayDataSet.Year(ListUnknown(i),1))]};
   end
end

% Write away the file
save('all')
v = genvarname([FileShort, FileSecond]);
eval([v '= ArrayDataSet;']);
OutputName = v;
cd(Current)
cd('MatlabLTEFiles')
save(OutputName,v)
cd(Current)
[LTETesterSelect,AverageYield,MaxYieldTreat,MinimumYieldTreat,LimitedCount] = ExtractYieldData(OutputName,LimitedCount,Current);
cd(Current)
cd('MatlabLTEFiles')
save(OutputName,v,'LTETesterSelect','AverageYield','MaxYieldTreat','MinimumYieldTreat')
cd(FileDir)
end
%%
%%
function [LTETesterSelect,AverageYield,MaxYieldTreat,MinimumYieldTreat,LimitedCount] = ExtractYieldData(File,LimitedCount,Current)
Filename = [File,'.mat'];
cd('MatlabLTEFiles')
load(Filename)
Array = eval(File);
FileShort = File(1:3);
ListExclude = find(strcmp(Array.Exclude,'Y')==1);
Array(ListExclude,:) = [];

LTETester = dataset({FileShort},'Varnames', 'LTE');
if length(File)>3
    SubLTE = File(5:end);
else
    SubLTE = 'One Crop';
end
LTETester.SubLTE = {SubLTE};
LTETester.Treatment = unique(Array.TreatType);
for i = 1: length(LTETester.Treatment)
    Tester = find(strcmp(Array.TreatType,LTETester.Treatment(i))==1);
    LTETester.SubLTE(i,1) = {SubLTE};
    LTETester.FullLTE(i,1) = {File};
    LTETester.Selected(i,1) = {'-'};
    LTETester.MeanYields(i,1) = nanmean(Array.Yield(Tester)); % for checking purposes
    LTETester.StartYear (i,1) = min(unique(Array.Year(Tester)));
    LTETester.EndYear (i,1) = max(unique(Array.Year(Tester)));
    LTETester.Years(i,1) = LTETester.EndYear (i,1) - LTETester.StartYear (i,1); % Allow for gaps because of rotations
    LTETester.Nutrients(i,1) = (Array.Nutrients_Manure(Tester(1,1)));
    LTETester.Residue(i,1) = (Array.Residue(Tester(1,1)));
    LTETester.Ctype(i,1) = (Array.Ctype(Tester(1,1)));
    LTETester.Crop(i,1) = (Array.Crop(Tester(1,1)));
    LTETester.LTE(i,1) = {FileShort};
end

MaximumYears = max(LTETester.Years);
LTETester.PercYears =  LTETester.Years./MaximumYears;
SeventyFive = find(LTETester.PercYears >= 0.75); % so present for least 75% of the LTE time
LTETesterSelect =  LTETester (SeventyFive,:);

clear LTETester SeventyFive MaximumYears i Tester ListExclude
SelectedTreatments = LTETesterSelect.Treatment;
for i = 1:length(SelectedTreatments)
    List = find(strcmp(Array.TreatType,SelectedTreatments(i))==1);
    if i == 1
        SelectedTreatmentsperLTE = Array(List,:);
    else
        SelectedTreatmentsperLTE = [SelectedTreatmentsperLTE;Array(List,:)];  %#ok<*AGROW>
    end
end
clear i List


% Stats
Yield = SelectedTreatmentsperLTE.Yield;
Year = SelectedTreatmentsperLTE.Year;
TreatType = SelectedTreatmentsperLTE.TreatType;
Replica1 = SelectedTreatmentsperLTE.Replica;
Replica2 = SelectedTreatmentsperLTE.SubReplica;

% Selection for coefficients per Treatment type
[~,StatsOutput,stats]= anovan((Yield),{Year, TreatType, Replica1, Replica2},'sstype',1,...
    'model',[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1],'nested',[0 0 0 0; 0 0 0 0; 0 1 0 0; 0 0 1 0],'display', 'off',...
    'continuous',1, 'random', [2 3 4],'varnames', {'Year', 'TreatType' , 'Replica', 'SubReplica'}); 
Endvalue = (3 + size(LTETesterSelect,1)) -1;
LTETesterSelect.NameCheck = stats.coeffnames(3:Endvalue);
LTETesterSelect.Coeffs = stats.coeffs(3:Endvalue);

cd(Current)
load('InputData.mat', 'LimitedTreatments');
% Per year statistics 
for i = 1:length(SelectedTreatments)
    List = find(strcmp(SelectedTreatmentsperLTE.TreatType,SelectedTreatments(i))==1);
    PerTreat = SelectedTreatmentsperLTE(List,:);
    Year = PerTreat.Year-1970;
    Replica1 = PerTreat.Replica;
    Replica2 = PerTreat.SubReplica;
    ReferenceYield = nanmedian(PerTreat.Yield);
    YieldMedian = PerTreat.Yield./ReferenceYield;
    [~,StatsOutput2,stats2]= anovan((YieldMedian),{Year, Replica1, Replica2},'sstype',1,...
        'model',[1 0 0; 0 1 0; 0 0 1],'nested',[0 0 0; 0 0 0; 0 1 0],'display', 'off',...
        'continuous',1, 'random', [2 3],'varnames', {'Year', 'Replica', 'SubReplica'});
    LTETesterSelect.YearConstant(i,1) =  stats2.coeffs(1);
    LTETesterSelect.YearCoef(i,1) =  stats2.coeffs(2);
    LTETesterSelect.Year_FValue(i,1) =  StatsOutput2(2,6);
    LTETesterSelect.Year_PValue(i,1) =  StatsOutput2(2,7);
end

LTETesterSelect(i+1,[1,2,4,11]) = LTETesterSelect(i,[1,2,4,11]);
LTETesterSelect.Years(i+1,1) = round(nanmedian(LTETesterSelect.Years));
LTETesterSelect.StartYear(i+1,1) = round(nanmedian(LTETesterSelect.StartYear));
LTETesterSelect.EndYear(i+1,1) = round(nanmedian(LTETesterSelect.EndYear));
LTETesterSelect.Treatment(i+1,1) = {'Combined'};
LTETesterSelect.NameCheck(i+1,1) = {'Combined'};
LTETesterSelect.YearConstant(i+1,1) =  stats.coeffs(1);
LTETesterSelect.YearCoef(i+1,1) =  stats.coeffs(2);
LTETesterSelect.Year_FValue(i+1,1) =  StatsOutput(2,6);
LTETesterSelect.Year_PValue(i+1,1) =  StatsOutput(2,7);

% Extract Limited treatments
Count = 0;
for is = 1:length(SelectedTreatments)
    Count = Count +1;
    ListLTE(Count) = find(strcmp(LimitedTreatments.FullLTE,File) & strcmp(LimitedTreatments.Treatment,SelectedTreatments(is)));
end
LTETesterSelect.Prob_NPKLimit(1:Count,1) = LimitedTreatments.Prob_NPKLimit(ListLTE);
LTETesterSelect.Maybe_NPKLimit(1:Count,1) = LimitedTreatments.Maybe_NPKLimit(ListLTE);
LTETesterSelect.Prob_MonoLimit(1:Count,1) = LimitedTreatments.Prob_MonoLimit(ListLTE);
LTETesterSelect.Maybe_MonoLimit(1:Count,1) = LimitedTreatments.Maybe_MonoLimit(ListLTE);

%% MEDIAN YIELD
% Refernce Yield Median
[AverageYield,FullYears] = CalculateYearMeans(SelectedTreatmentsperLTE,-999,File);
AverageYield = ReadjustVarieties(AverageYield,FullYears);
AverageYield = SetAgainstMedian(AverageYield);
%% MAXIMUM YIELD TREATMENT
NoCombi = (length(LTETesterSelect.Coeffs)-1);
MaxTreat = find(LTETesterSelect.Coeffs(1:NoCombi) == max(LTETesterSelect.Coeffs(1:NoCombi)));
ListMax = find(strcmp(SelectedTreatmentsperLTE.TreatType,SelectedTreatments(MaxTreat))==1); %#ok<*FNDSB>
SelectedforMax = SelectedTreatmentsperLTE(ListMax,:);
[MaxYieldTreat,~] = CalculateYearMeans(SelectedforMax,FullYears,File);
MaxYieldTreat = ReadjustVarieties(MaxYieldTreat,FullYears);
MaxYieldTreat = SetAgainstMedian(MaxYieldTreat);
LTETesterSelect.Selected(MaxTreat,1) = {'Maximum Treatment'} ;
%% NO NUTRIENT OR MINIMUM YIELD]
NoCombi = (length(LTETesterSelect.Coeffs)-1);
Text = {'Min Limited Treatment'};
NoNutMin = find(strcmp(LTETesterSelect.Prob_NPKLimit(1:NoCombi),'limited')); 
if isempty(NoNutMin) == 1
    NoNutMin = find(strcmp(LTETesterSelect.Maybe_NPKLimit(1:NoCombi),'limited'));
end
if isempty(NoNutMin) == 1
    NoNutMin = 1:NoCombi;
    Text = {'Min coefficient Treatment'};
else
    LimitedCount.NPKNr = LimitedCount.NPKNr + 1;
    LimitedCount.NPK(LimitedCount.NPKNr,1) = {File};
end
NoNutListMin = find(LTETesterSelect.Coeffs(1:NoCombi) == min(LTETesterSelect.Coeffs(NoNutMin)));
ListMin = find(strcmp(SelectedTreatmentsperLTE.TreatType,LTETesterSelect.Treatment(NoNutListMin))==1);
SelectedforMin= SelectedTreatmentsperLTE(ListMin,:);
[MinimumYieldTreat,~] = CalculateYearMeans(SelectedforMin,FullYears,File);
MinimumYieldTreat = ReadjustVarieties(MinimumYieldTreat,FullYears);
MinimumYieldTreat = SetAgainstMedian(MinimumYieldTreat);
LTETesterSelect.Selected(NoNutListMin,1) = Text;
if NoCombi == 1
    LTETesterSelect.Selected(NoNutListMin,1) = {'Min and Max Treatment'};
end

% Make a list for Monoculture Limited
NoNutList = find(strcmp(LTETesterSelect.Prob_MonoLimit(1:NoCombi),'mono_not_limited')); 
% if isempty(NoNutList) == 1
%     NoNutList = find(strcmp(LTETesterSelect.Maybe_MonoLimit(1:NoCombi),'mono_not_limited'));
% end
if isempty(NoNutList) ~= 1
    LimitedCount.MonoNr = LimitedCount.MonoNr + 1;
    LimitedCount.Mono(LimitedCount.MonoNr,1) = {File};
end
end % of function

%%
function [Yields,ListerYears] = CalculateYearMeans(InputTreatmentYield,ListerYears,File)
ReferenceYield = nanmedian(InputTreatmentYield.Yield);

if ListerYears == -999
    ListerYearMin = min(InputTreatmentYield.Year);
    if ListerYearMin < 1965
        ListerYearMin = 1965;
    end
    ListerYearMax = max(InputTreatmentYield.Year);
    if ListerYearMax > 2022
        ListerYearMax = 2022;
    end
    ListerYears = ListerYearMin:1:ListerYearMax;
end

TreatTyper = unique(InputTreatmentYield.TreatType);
if length(TreatTyper)>1
    TreatTyper = {'Combined'};
end
FileShort = File(1:3);
Yields = dataset({FileShort},'Varnames', 'LTE'); %#ok<*DTSET>
if length(File)>3
    SubLTE = File(5:end);
else
    SubLTE = 'One Crop';
end
Yields.SubLTE = {SubLTE};
warning off 
for i = 1: length(ListerYears)
    ValuesPerYear = find(InputTreatmentYield.Year ==ListerYears(i));
    Yields.LTE(i,1) = {FileShort};
    Yields.SubLTE(i,1) = {SubLTE};
    Yields.FullLTE(i,1) = {File};
    if isempty(ValuesPerYear) ~= 1
        PerYearPerLTE = InputTreatmentYield(ValuesPerYear,:);
        Tester = unique(PerYearPerLTE.TreatType);
        PerTreat = [];
        PerTreatNorm = [];
        for x = 1:length(Tester)
            ListAvg = find(strcmp(PerYearPerLTE.TreatType,Tester(x))==1);
            PerTreat(x) = nanmedian(PerYearPerLTE.Yield(ListAvg));
            % Normalised
            PerTreatNorm(x) = nanmedian(PerYearPerLTE.Yield(ListAvg)./ReferenceYield);
        end
        Yields.TreatType(i,1) = TreatTyper;
        Yields.Year(i,1) = ListerYears(i);
        Variety = unique(PerYearPerLTE.Variety);
        if length(Variety) >1 && (strcmp(FileShort,{'IRR'}) ~= 1) && (strcmp(FileShort,{'C29'}) ~= 1) &&...
                (strcmp(FileShort,{'KBM'}) ~= 1) && (strcmp(FileShort,{'KRC'}) ~= 1) &&...
                (strcmp(FileShort,{'URI'}) ~= 1) && (strcmp(FileShort,{'DCA'}) ~= 1) &&...
                (strcmp(FileShort,{'DIA'}) ~= 1) && (strcmp(FileShort,{'DMI'}) ~= 1) &&...
                (strcmp(FileShort,{'DWA'}) ~= 1) && (strcmp(FileShort,{'DMD'}) ~= 1)
            disp(['Multiple varieties in a year in ',File])
            disp(ListerYears(i))
        end
        Yields.Variety(i,1) = Variety(1);
        Yields.Ctype(i,1) = PerYearPerLTE.Ctype(1);
        Yields.Crop(i,1) = PerYearPerLTE.Crop(1);
        Yields.NrValuesYear(i,1) = length(PerTreat);
        Yields.MedianYield(i,1) = nanmedian(PerTreat); %#ok<*NANMEDIAN>
        Yields.MedianYieldCV(i,1) = nanstd(PerTreat)./nanmean(PerTreat); %#ok<*NANMEAN,*NANSTD>
        Yields.MedianYieldNorm(i,1) = nanmedian(PerTreatNorm);
        clear PerTreat PerTreatNorm
    else
        Yields.TreatType(i,1) = TreatTyper;
        Yields.Year(i,1) = ListerYears(i);
        Yields.Ctype(i,1) = InputTreatmentYield.Ctype(1);
        Yields.Crop(i,1) = InputTreatmentYield.Crop(1);
        Yields.Variety(i,1) = {'No crop'};
        Yields.NrValuesYear(i,1) = NaN;
        Yields.MedianYield(i,1) = NaN;
        Yields.MedianYieldCV(i,1) = NaN;
        Yields.MedianYieldNorm(i,1) = NaN;
    end
end
end % of function


%% IMPORTANT STEP Set against the median per variety
function Yields = SetAgainstMedian(Yields)
% Set against the Median per variety
VarList = unique(Yields.AdjustedVariety);
for Var = 1:1:length(VarList)
    list = find(strcmp(Yields.AdjustedVariety,VarList(Var))==1);
    if strcmp(VarList(Var),'-') ~= 1 || strcmp(VarList(Var),'No crop') ~= 1
        Yield = double(Yields.MedianYieldNorm(list));
        MedianYield = nanmedian(Yield);
        YieldMedian = Yield/MedianYield;
        for Ind = 1:length(YieldMedian)
        if YieldMedian(Ind) > 2
           YieldMedian(Ind) = 2;
        end
        end
        Yields.MedianYieldVariety(list,1) = YieldMedian;
    else
        Yields.MedianYieldVariety(list,1) = NaN;
    end
end
end % of function

%% RE_ADJUST VARIETIES
% This re-adjust varieties that don't occur at least 3 years.
% if for less than three years, they are joined into parts of three years
% by adding to the previous. This reduces the amount of single varieties.
% Such are either added to more occuring other ones or when all single, put
% in blocks of three years. This includes rotations, so three years can be
% spread in time (e.g. 2005; 2008 & 2011, when occuring every three years)

function Yields = ReadjustVarieties(Yields,ListerYears)
TotVar = 1;
CountExist = 1;
countpresence = 1;
for i = 1: length(ListerYears)
    Variety = Yields.Variety(i);
    if strcmp(Variety,{'No crop'})==1
        AdjustedVariety(i,1) = {'No crop'};
    else
        if countpresence == 1 % Intitate
            AdjustedVariety(i,1) = Variety;
            VarExist(CountExist) = Variety;
            VarExistCount(CountExist) = 1;
            CountExist = CountExist + 1;
            countpresence = countpresence + 1;
        else
            Test = find(strcmp(VarExist,Variety)==1);
            if isempty(Test) ~= 1 % this is an already existing variety
                AdjustedVariety(i,1) = VarExist(Test);
                VarExistCount(Test) = VarExistCount(Test) +1;
                TotVar = TotVar + 1;
            else
                if TotVar <3 % this one is new but will be added to the previous
                    % equals in case this variety has multiple, the previous 
                    % incomplete one (less than 3) is added to this one.
                    AdjustedVariety(i,1) = Previous;
                    TotVar = TotVar + 1;
                    % if this one is longer run (3 or more), give this name instead
                    Tester = sum(strcmp(Yields.Variety,Variety));
                    if Tester >=3
                        AdjustedVariety(i,1) = Variety;
                        Lister = find(strcmp(AdjustedVariety,Previous)==1);
                        AdjustedVariety(Lister,1) = Variety;
                        VarExist(CountExist-1) = Variety;
                    end

                else % This is a new variety
                    AdjustedVariety(i,1) = Variety;
                    VarExist(CountExist) = AdjustedVariety(i,1);
                    VarExistCount(CountExist) = 1;
                    CountExist = CountExist + 1;
                    TotVar = 1;
                end
            end
        end % countpresence == 1
        Previous = AdjustedVariety(i,1);
    end % no crop
end % length listeryears

% Correct for glitch single varieties among majorities by setting to
% previous (e.g. one time occuring among before and after the same one)
UniList = unique(AdjustedVariety);
Test = find(strcmp(UniList,{'No crop'})==1);
if isempty(Test)~=1
UniList(Test) = [];
end
for i = 1:length(UniList)
CountPerVar(i) = sum(strcmp(AdjustedVariety,UniList(i)));
end
Misslist = find(CountPerVar<3);
if isempty(Misslist)~=1
    for i = 1:length(Misslist)
        SingeList = find(strcmp(AdjustedVariety,UniList(Misslist(i)))==1);
        for x = 2:length(SingeList)
            AdjustedVariety(SingeList(x)) = AdjustedVariety(SingeList(x-1));
        end
    end
end
for i = 1: length(ListerYears)
    if strcmp(AdjustedVariety(i),{'No crop'})==1
        AdjustedVariety(i) = Yields.Variety(i);
    end
end
Yields.AdjustedVariety = AdjustedVariety;
end % of function







