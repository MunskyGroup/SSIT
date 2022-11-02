function app = updateModelSolveAndCompareToData(app,computeSensitivity)
% This function Projects the FSP results onto species of interest.
% This compares the data to the FSP results.
arguments
    app
    computeSensitivity = false;
end
if isempty(app.FspConstraintTable.Data)
    makeDefaultConstraints(app);
    readConstraintsForAdaptiveFsp(app);
end

if ~isfield(app.DataLoadingAndFittingTabOutputs,'ModelMerge')
    app.ModelParameterTable.Data=app.fit_parameters_table.Data(:,1:2);
    app.ReactionsTabOutputs.parameters=app.fit_parameters_table.Data(:,1:2);
end
app.ReactionsTabOutputs.parameters=app.ModelParameterTable.Data(:,1:2);

if computeSensitivity
    app.SensPrintTimesEditField.Value = app.FspPrintTimesField.Value;
    app = runSensitivity(app,false);
    solutions = app.SensFspTabOutputs.solutions.data;
else
    app = runFsp(app);
    solutions = app.FspTabOutputs.solutions;
end
app = ssit.pdo.generatePDO(app);

%% Project FSP result onto species of interest.
T_array = unique(eval(app.FspPrintTimesField.Value));
Nd = solutions{1}.p.dim;
for i=1:Nd
    indsPlots(i) = max(contains(app.SpeciesForFitPlot.Value,app.SpeciesForFitPlot.Items{i}));
end

indsUnobserved = setdiff([1:Nd],find(indsPlots));

szP = zeros(1,Nd);
for it = length(T_array):-1:1
    szP = max(szP,size(solutions{it}.p.data));
end

P = zeros([length(T_array),szP(indsPlots)]);
for it = length(T_array):-1:1
    if ~isempty(solutions{it})

        px = solutions{it}.p;

        if computeSensitivity
            Sx = solutions{it}.S;
            parCount = length(Sx);
            % Add effect of PDO.
            if ~isempty(app.FIMTabOutputs.distortionOperator)
                for iPar = 1:parCount
                    Sx(iPar) = app.FIMTabOutputs.distortionOperator.computeObservationDistDiff(px, Sx(iPar), iPar);
                end
            end
        end

        % Add effect of PDO.
        if ~isempty(app.FIMTabOutputs.distortionOperator)
            px = app.FIMTabOutputs.distortionOperator.computeObservationDist(px);
        end

        if ~isempty(indsUnobserved)
            d = double(px.sumOver(indsUnobserved).data);
        else
            d = double(px.data);
        end

        P(it,d~=0) = d(d~=0);

        if computeSensitivity
            for iPar = parCount:-1:1
                d = double(Sx(iPar).sumOver(indsUnobserved).data);
                S{iPar}(it,1:length(d)) = d;
            end
        end
    end
end

%% Padd P or Data to match sizes of tensors.
NP = size(P);
NDat = size(app.DataLoadingAndFittingTabOutputs.dataTensor);
if length(NP)<Nd; NP(end+1:Nd)=1; end
if max(NDat(2:end)-NP(2:length(NDat)))>0   % Pad if data longer than model
    NP(2:length(NDat)) = max(NP(2:length(NDat)),NDat(2:end));
    tmp = 'P(end';
    for j = 2:length(NDat)
        tmp = [tmp,',NP(',num2str(j),')'];
    end
    tmp = [tmp,')=0;'];
    eval(tmp)
end
if max(NP(2:length(NDat))-NDat(2:end))>0   % truncate if model longer than data
    tmp = 'P = P(:';
    for j = 2:length(NDat)
        tmp = [tmp,',1:',num2str(NDat(j))];
    end
    for j = (length(NDat)+1):4
        tmp = [tmp,',1'];
    end
    tmp = [tmp,');'];
    eval(tmp)
end
P = max(P,1e-10);

if computeSensitivity
    for iPar = parCount:-1:1
        NS = size(S{iPar});
        if length(NS)<4; NS(end+1:4)=1; end
        if max(NDat(2:end)-NS(2:length(NDat)))>0   % Pad if data longer than model
            NS(2:length(NDat)) = max(NS(2:length(NDat)),NDat(2:end));
            S{iPar}(end,NS(2),NS(3),NS(4)) = 0;
        end
        if max(NS(2:length(NDat))-NDat(2:end))>0   % truncate if model longer than data
            tmp = 'S{iPar} = S{iPar}(:';
            for j = 2:length(NDat)
                tmp = [tmp,',1:',num2str(NDat(j))];
            end
            for j = (length(NDat)+1):4
                tmp = [tmp,',1'];
            end
            tmp = [tmp,');'];
            eval(tmp)
        end
    end
end

%%
times = app.DataLoadingAndFittingTabOutputs.fittingOptions.fit_times;

%% Compute log likelihood using equal sized P and Data tensors.
sz = size(P);
app.DataLoadingAndFittingTabOutputs.fitResults.current = zeros([length(times),sz(2:end)]);
app.DataLoadingAndFittingTabOutputs.fitResults.currentData = zeros([length(times),sz(2:end)]);
LogLk = zeros(1,length(times));
numCells = zeros(1,length(times));
perfectMod = zeros(1,length(times));
perfectModSmoothed = zeros(1,length(times));
for i=1:length(times)
    [~,j] = min(abs(T_array-times(i)));
    Jind = app.DataLoadingAndFittingTabOutputs.dataTensor.subs(:,1) == i;
    SpInds = app.DataLoadingAndFittingTabOutputs.dataTensor.subs(Jind,:);
    SpVals = app.DataLoadingAndFittingTabOutputs.dataTensor.vals(Jind);
    H = sptensor([ones(length(SpVals),1),SpInds(:,2:end)],SpVals,[1,NDat(2:end)]);
    H = double(H);
    Pt = P(j,:,:,:,:,:,:);
    LogLk(i) = sum(H(:).*log(Pt(:)));
    Q = H(:)/sum(H(:));
    smQ = smooth(Q);
    logQ = log(Q); logQ(H==0)=1;
    logSmQ = log(smQ); logSmQ(H==0)=1;
    perfectMod(i) = sum(H(:).*logQ);
    perfectModSmoothed(i) = sum(H(:).*logSmQ);
    numCells(i) = sum(H(:));
    app.DataLoadingAndFittingTabOutputs.fitResults.current(i,:,:,:) = Pt;
    app.DataLoadingAndFittingTabOutputs.fitResults.currentData(i,:,:,:) = ...
         reshape(Q,size(app.DataLoadingAndFittingTabOutputs.fitResults.currentData(i,:,:,:)));
    if computeSensitivity
        for iPar = parCount:-1:1
            St = S{iPar}(j,:,:,:);
            dlogL_dPar(iPar,i) = sum(H(:).*St(:)./Pt(:));
        end
    end
end
app.DataLoadingAndFittingTabOutputs.J_LogLk = sum(LogLk);
if computeSensitivity
    app.DataLoadingAndFittingTabOutputs.dLogLk_dpar = sum(dlogL_dPar,2);
end
app.DataLoadingAndFittingTabOutputs.V_LogLk = LogLk;
app.DataLoadingAndFittingTabOutputs.numCells = numCells;
app.DataLoadingAndFittingTabOutputs.perfectMod = perfectMod;
app.DataLoadingAndFittingTabOutputs.perfectModSmoothed = perfectModSmoothed;

