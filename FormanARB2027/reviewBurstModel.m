%% Plots of means and fano factor versus parameters.

kr = 100;
gr = 1;
N = 100;
kArr = logspace(-2,3,N);
MEAN = zeros(N,N);
FANO = zeros(N,N);
for i = 1:N
    for j = 1:N
        F(i,j) = kArr(i)/(kArr(i)+kArr(j));
        MEAN(i,j) = kr/gr*F(i,j);
        FANO(i,j) = 1+((1-F(i,j))*kr)/(kArr(i)+kArr(j)+gr);
    end
end
figure(1); clf;
ax1 = axes;
[~,cF] = contourf(ax1,log10(kArr),log10(kArr),log10(FANO),linspace(0,log10(64),30)); hold on;

colormap(ax1,jet)
cb = colorbar(ax1);
cb.Ticks = log10([1.0001,4,16,64]);
cb.TickLabels = {'1','4','16','64'};

set(gca,'xtick',[-2:3],'ytick',[-2:3],...
    'XTickLabel',{'10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}'},...
    'YTickLabel',{'10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}'},...
    'FontSize',16);

ax2 = axes;
[~,cM] = contour(ax2,log10(kArr),log10(kArr),MEAN,[2,25,75],'k--','LineWidth',2);
ax2.Color = 'none';
ax2.XColor = 'none';
ax2.YColor = 'none';
linkaxes([ax1 ax2]);
ax2.Position = ax1.Position;

%% Burst Frequency Model

Model = SSIT('Empty');

Model.species = {'gON','mRNA'};

Model.initialCondition = [0;0];

%                                     PRIOR
Model.parameters = {'kon0',0.01;...  % logn(-1,2)
    'koff0',0.1;...                   % logn(0,2)
    'kr',10;...                     % logn(1,2)
    'g',0.01;...                    % logn(-2,2)
    'kD',10;...                     % logn(1,2)
    'S0',1;...                      % NA (initial input concentration)
    'S1',5};                        % NA (final input concentration)

Model.inputExpressions = {'Iu','S0+(S1-S0)*(t>0)'};

Model = Model.addReaction(struct(...
    'propensity',{'kon0*(Iu/(kD+Iu))*(1-gON)'},...
    'stoichiometry',{{'gON',1}}));

Model = Model.addReaction(struct(...
    'propensity',{'koff0*gON'},...
    'stoichiometry',{{'gON',-1}}));

Model = Model.addReaction(struct(...
    'propensity',{'kr*gON'},...
    'stoichiometry',{{'mRNA',1}}));

Model = Model.addReaction(struct(...
    'propensity',{'g*mRNA'},...
    'stoichiometry',{{'mRNA',-1}}));

Model.fspOptions.initApproxSS = true;
Model.tSpan = linspace(0,300,31);

Model = Model.solve;

Model.plotFSP

%% Verification of FIM using CRLB (spread of MLE)
nCellsInExperiment = 0*Model.tSpan;
nCellsInExperiment([1,11,31]) = 200;
nMLE = 10;
Model.fittingOptions.modelVarsToFit = [1:5];
MLE = Model.estimateMLEspread(nCells=nCellsInExperiment,observableSpecies={'mRNA'},nMLE=nMLE,simsSaveFile='BurstFIMSims.csv',freePars=[1:5],restart=true);
MLE = Model.estimateMLEspread(nCells=nCellsInExperiment,observableSpecies={'mRNA'},nMLE=nMLE,simsSaveFile='BurstFIMSims.csv',freePars=[1:5],startPars=exp(MLE.mhSamples),restart=false);

FIMs = Model.computeFIM(scale='log',freePars=[1:5],...
    observed={'mRNA'});
FIMTotal = Model.totalFim(FIMs,nCellsInExperiment);

Model.plotMHResults(MLE,FIM=FIMTotal,fimScale='log',truncateChain=false);

return

%% FIM Calculations
Sarray = [1:5];

Model.solutionScheme = 'fspsens';
FIM = cell(length(Sarray),length(Sarray),length(Model.tSpan));
for iS0 = 1:length(Sarray)
    for iS1 = 1:length(Sarray)
        Model = Model.changeParameter({'S0',Sarray(iS0);'S1',Sarray(iS1)-Sarray(iS0)});
        Model = Model.solve;
        FIM(iS0,iS1,:) = Model.computeFIM(freePars=(1:4),scale='log');
    end
end

%% FIM for different experiment designs.
Ncells = 600;
% Measurment at one steady state values.
iS = 3;
FIM_One_SS = Ncells*FIM{iS,iS,1};
disp(['Determinant of FIM for one SS measurement: ',num2str(det(FIM_One_SS))])

% Measurment at two steady state values.
iS1 = 1;
iS2 = 5;
FIM_Two_SS = Ncells/2*FIM{iS1,iS1,1}+Ncells/2*FIM{iS2,iS2,1};
disp(['Determinant of FIM for two SS measurements: ',num2str(det(FIM_Two_SS))])

% Measurement at change from one SS to another at three time points.
iS1 = 1;
iS2 = 5;
itimes = [1,11,31];
FIM_Dynamic = 0;
for it = 1:length(itimes)
    FIM_Dynamic = FIM_Dynamic + Ncells/length(itimes)*FIM{iS1,iS2,itimes(it)};
end
disp(['Determinant of FIM for dynamic measurements: ',num2str(det(FIM_Dynamic))])

%% Optimized Experiment Design
% AllFims = reshape(FIM,numel(FIM),1);
allFims = {};%cell(numel(FIM),1);
indsFims = [];zeros(numel(FIM),3);
k = 0;
for iS0 = 1:length(Sarray)
    for iS1 = 1:length(Sarray)
        for iT = 1:length(Model.tSpan)
            k = k+1;
            allFims(k,1) = FIM(iS0,iS1,iT);
            indsFims(k,:) = [iS0,iS1,iT];
        end
    end
end
OptExperiment = Model.optimizeCellCounts(allFims,600,'D-opt');
J = find(OptExperiment);
disp(['Optimized Experiment Design:'])
for j = 1:length(J)
    % Store optimized parameters and their corresponding indices
    optimizedParams = OptExperiment(J(j));
    paramIndices = indsFims(J(j), :);
    if paramIndices(3)==1||paramIndices(1)==paramIndices(2) % SS experiment
        disp(['   ',num2str(optimizedParams),' cells at steady state for S0 = ',num2str(Sarray(paramIndices(1)))])
    else
        disp(['   ',num2str(optimizedParams),' cells at time ',num2str(Model.tSpan(paramIndices(3))),' for S0 = ',num2str(Sarray(paramIndices(1))),' and S1 = ',num2str(Sarray(paramIndices(2)))])
    end
end

FIM_Opt = 0;
for i = 1:length(OptExperiment)
    FIM_Opt = FIM_Opt + OptExperiment(i)*allFims{i};
end
disp(['Determinant of FIM for optimized measurements: ',num2str(det(FIM_Opt))])

return
%% Generalized model
% 
% In the most general form of this model, the parameter 'alpha' determines
% if the model is burst size regulate (koff, alpha=0) or burst frequency
% regulated (kon, alpha=1) or something in between.
ModelGen = SSIT('Empty');

ModelGen.species = {'gON','mRNA'};

ModelGen.initialCondition = [0;0];

%                                     PRIOR
ModelGen.parameters = {'kon0',0.1;...  % logn(-1,2)
    'koff0',1;...                   % logn(0,2)
    'kr',10;...                     % logn(1,2)
    'g',0.01;...                    % logn(-2,2)
    'kD',10;...                     % logn(1,2)
    'alpha',1e-6;...                   % n(0,2)
    'S0',1;...                      % NA (initial input concentration)
    'S1',5};                        % NA (final input concentration)

ModelGen.inputExpressions = {'Iu','S0+S1*(t>0)'};

ModelGen = ModelGen.addReaction(struct(...
    'propensity',{'kon0*((1-alpha)+(alpha)*Iu/(kD+Iu))*(1-gON)'},...
    'stoichiometry',{{'gON',1}}));

ModelGen = ModelGen.addReaction(struct(...
    'propensity',{'koff0*((alpha)+(1-alpha)*(kD+Iu)/Iu)*gON'},...
    'stoichiometry',{{'gON',-1}}));

ModelGen = ModelGen.addReaction(struct(...
    'propensity',{'kr*gON'},...
    'stoichiometry',{{'mRNA',1}}));

ModelGen = ModelGen.addReaction(struct(...
    'propensity',{'g*mRNA'},...
    'stoichiometry',{{'mRNA',-1}}));

ModelGen.tSpan = linspace(0,300,31);

ModelGen = ModelGen.solve;

ModelGen.plotFSP

%% Parameters for key points.
kr = 100;
gr = 1;
KD = 10;
kon0a = 2.5;
koff0a = 10;
S0 = KD/24;
S1 = 23*KD/24;

%                                     PRIOR
ModelGen.parameters = {'kon0',kon0;...  % logn(-1,2)
    'koff0',koff0;...                   % logn(0,2)
    'kr',kr;...                     % logn(1,2)
    'g',gr;...                    % logn(-2,2)
    'kD',KD;...                     % logn(1,2)
    'alpha',-6;...                   % n(0,2)
    'S0',S0;...                      % NA (initial input concentration)
    'S1',S1};                        % NA (final input concentration)

ModelGen.tSpan = linspace(0,6,61);

ModelGen.fspOptions.initApproxSS = true;
ModelGen = ModelGen.solve;

ModelGen.plotFSP
