function [NumbersFixedStore, IterationListStore,PLSOutStore] =  ForwardRegressionAlgorithm(XVar, InVar, EnvironmentNLM,NameList,Weights,AxNum) 
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