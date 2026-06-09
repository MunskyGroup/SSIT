%% Install SSIT Tools
% Make sure that you are in the main SSIT directory before running.
% We recommend running this one cell at a time so that results are not
% overwritten.
install

%% Example 1 - Toggle Switch
% Create Toggle Model
clear
Toggle = SSIT();
Toggle.species = {'Gx_on','Gx_off','Gy_on','Gy_off','Px','Py'};
Toggle.initialCondition = [1;0;1;0;0;0];
Toggle.parameters = {'sxy',50;
    'dxy',1;
    'bxy',1e-4;
    'uxy',0.1};
Toggle.propensityFunctions = {};
Toggle.stoichiometry = [];
Toggle = Toggle.addReaction(struct('propensity',{'sxy*Gx_on'},'stoichiometry',{{'Px',1}}));
Toggle = Toggle.addReaction(struct('propensity',{'sxy*Gy_on'},'stoichiometry',{{'Py',1}}));
Toggle = Toggle.addReaction(struct('propensity',{'dxy*Px'},'stoichiometry',{{'Px',-1}}));
Toggle = Toggle.addReaction(struct('propensity',{'dxy*Py'},'stoichiometry',{{'Py',-1}}));
Toggle = Toggle.addReaction(struct('propensity',{'bxy*Px*(Px-1)*Gy_on'},'stoichiometry',{{'Px',-2;'Gy_on',-1;'Gy_off',1}}));
Toggle = Toggle.addReaction(struct('propensity',{'bxy*Py*(Py-1)*Gx_on'},'stoichiometry',{{'Py',-2;'Gx_on',-1;'Gx_off',1}}));
Toggle = Toggle.addReaction(struct('propensity',{'uxy*Gy_off'},'stoichiometry',{{'Px',2;'Gy_on',1;'Gy_off',-1}}));
Toggle = Toggle.addReaction(struct('propensity',{'uxy*Gx_off'},'stoichiometry',{{'Py',2;'Gx_on',1;'Gx_off',-1}}));
Toggle.tSpan = linspace(0,40,81);
Toggle.fspOptions.verbose = true;
Toggle = Toggle.formPropensitiesGeneral('Toggle');

%% Solve Toggle Model
Toggle.fspOptions.fspTol = 0.001;
tic
[~,~,Toggle] = Toggle.solve;
toc
Toggle.plotFSP(speciesNames={'Px','Py'});

%% Example 2 - MAPK Model
% Set up MAPK Model
clear
MAPK = SSIT();

% Species
MAPK.species =   {'MKP3'      
    'Kpp_MKP3'  
    'KpT_MKP3_Y'
    'KpT_MKP3_T'
    'K_MKP3_T'  
    'KpY_MKP3'  
    'K_MKP3_Y'  
    'KpT_MEK'   
    'KpT'       
    'KpY_MEK'   
    'K_MEK_T'   
    'KpY'       
    'K_MEK_Y'   
    'K'         
    'MEK'       
    'Kpp'       };

MAPK.initialCondition = zeros(length(MAPK.species),1);
MAPK.initialCondition(find(strcmp(MAPK.species,'MKP3'))) = 1; % MKP3 IC is 1
MAPK.initialCondition(find(strcmp(MAPK.species,'K'))) = 3; % K IC is 3

% Parameters
MAPK.parameters = {
    's1',0.00024; 
    'd1',0.0001;
    's2',0.001; 
    'd2',0.15;
    's3',0.005; 
    'k1',0.375; 
    'k_1',1.0;
    'k2',0.06;
    'k3',0.375; 
    'k_3',1.0;
    'k4',4.5;
    'k5',0.375; 
    'k_5',1.0;
    'k6',0.06;
    'k7',0.375; 
    'k_7',1.0;
    'k8',4.5;
    'h1',0.015; 
    'h_1',1.0;
    'h2',0.032;
    'h3',0.31; 
    'h_3',0.01;
    'h4',0.01; 
    'h_4',1.0;
    'h5',0.5;
    'h6',0.086; 
    'h_6',0.0011;
    'h7',0.01; 
    'h_7',1.0;
    'h8',0.47;
    'h9',0.14; 
    'h_9',0.0018;
};

% Reactions
MAPK.propensityFunctions = {};
MAPK.stoichiometry = [];

% R1: synthesis/degradation of MEK
MAPK = MAPK.addReaction(struct('propensity',{'s2'},'stoichiometry',{{'MEK',1}}));
MAPK = MAPK.addReaction(struct('propensity',{'d2*MEK'},'stoichiometry',{{'MEK',-1}}));

% R2: synthesis/degradation of MEK
MAPK = MAPK.addReaction(struct('propensity',{'s3*Kpp'},'stoichiometry',{{'MEK',1}}));

% R3: synthesis/degradation of K
MAPK = MAPK.addReaction(struct('propensity',{'s1'},'stoichiometry',{{'K',1}}));
MAPK = MAPK.addReaction(struct('propensity',{'d1*K'},'stoichiometry',{{'K',-1}}));

% R4:
MAPK = MAPK.addReaction(struct('propensity',{'d1*KpY'},'stoichiometry',{{'KpY',-1}}));
MAPK = MAPK.addReaction(struct('propensity',{'d1*KpT'},'stoichiometry',{{'KpT',-1}}));
MAPK = MAPK.addReaction(struct('propensity',{'d1*Kpp'},'stoichiometry',{{'Kpp',-1}}));

% R5: K + MEK -> K_MEK_Y
MAPK = MAPK.addReaction(struct('propensity',{'k1*K*MEK'},'stoichiometry',{{'K',-1;'MEK',-1;'K_MEK_Y',1}}));

% R6: K_MEK_Y -> KpY + MEK, k2 = 0.06
MAPK = MAPK.addReaction(struct('propensity',{'k2*K_MEK_Y'},'stoichiometry',{{'K_MEK_Y',-1;'KpY',1;'MEK',1}}));

% R7: KpY + MEK <=> KpY_MEK, k3 = 0.375, k-3 = 1.0
MAPK = MAPK.addReaction(struct('propensity',{'k3*KpY*MEK'},'stoichiometry',{{'KpY',-1;'MEK',-1;'KpY_MEK',1}}));
MAPK = MAPK.addReaction(struct('propensity',{'k_3*KpY_MEK'},'stoichiometry',{{'KpY_MEK',-1;'KpY',1;'MEK',1}}));

% R8: KpY_MEK -> Kpp + MEK, k4 = 4.5
MAPK = MAPK.addReaction(struct('propensity',{'k4*KpY_MEK'},'stoichiometry',{{'KpY_MEK',-1;'Kpp',1;'MEK',1}}));

% R9: K + MEK <=> K_MEK_T, k5 = 0.375, k-5 = 1.0
MAPK = MAPK.addReaction(struct('propensity',{'k5*K*MEK'},'stoichiometry',{{'K',-1;'MEK',-1;'K_MEK_T',1}}));
MAPK = MAPK.addReaction(struct('propensity',{'k_5*K_MEK_T'},'stoichiometry',{{'K_MEK_T',-1;'K',1;'MEK',1}}));

% R10: K_MEK_T -> KpT + MEK, k6 = 0.06
MAPK = MAPK.addReaction(struct('propensity',{'k6*K_MEK_T'},'stoichiometry',{{'K_MEK_T',-1;'KpT',1;'MEK',1}}));

% R11: KpT + MEK <=> KpT_MEK, k7 = 0.375, k-7 = 1.0
MAPK = MAPK.addReaction(struct('propensity',{'k7*KpT*MEK'},'stoichiometry',{{'KpT',-1;'MEK',-1;'KpT_MEK',1}}));
MAPK = MAPK.addReaction(struct('propensity',{'k_7*KpT_MEK'},'stoichiometry',{{'KpT_MEK',-1;'KpT',1;'MEK',1}}));

% R12: KpT_MEK -> Kpp + MEK, k8 = 4.5
MAPK = MAPK.addReaction(struct('propensity',{'k8*KpT_MEK'},'stoichiometry',{{'KpT_MEK',-1;'Kpp',1;'MEK',1}}));

% R13: Kpp + MKP3 <=> Kpp_MKP3, h1 = 0.015, h-1 = 1.0
MAPK = MAPK.addReaction(struct('propensity',{'h1*Kpp*MKP3'},'stoichiometry',{{'Kpp',-1;'MKP3',-1;'Kpp_MKP3',1}}));
MAPK = MAPK.addReaction(struct('propensity',{'h_1*Kpp_MKP3'},'stoichiometry',{{'Kpp_MKP3',-1;'Kpp',1;'MKP3',1}}));

% R14: Kpp_MKP3 -> KpT_MKP3_Y, h2 = 0.032
MAPK = MAPK.addReaction(struct('propensity',{'h2*Kpp_MKP3'},'stoichiometry',{{'Kpp_MKP3',-1;'KpT_MKP3_Y',1}}));

% R15: KpT_MKP3_Y <=> KpT + MKP3, h3 = 0.31, h-3 = 0.01
MAPK = MAPK.addReaction(struct('propensity',{'h3*KpT_MKP3_Y'},'stoichiometry',{{'KpT_MKP3_Y',-1;'KpT',1;'MKP3',1}}));
MAPK = MAPK.addReaction(struct('propensity',{'h_3*KpT*MKP3'},'stoichiometry',{{'KpT',-1;'MKP3',-1;'KpT_MKP3_Y',1}}));

% R16: KpT + MKP3 <=> KpT_MKP3_T, h4 = 0.01, h-4 = 1.0
MAPK = MAPK.addReaction(struct('propensity',{'h4*KpT*MKP3'},'stoichiometry',{{'KpT',-1;'MKP3',-1;'KpT_MKP3_T',1}}));
MAPK = MAPK.addReaction(struct('propensity',{'h_4*KpT_MKP3_T'},'stoichiometry',{{'KpT_MKP3_T',-1;'KpT',1;'MKP3',1}}));

% R17: KpT_MKP3_T -> K_MKP3_T, h5 = 0.5
MAPK = MAPK.addReaction(struct('propensity',{'h5*KpT_MKP3_T'},'stoichiometry',{{'KpT_MKP3_T',-1;'K_MKP3_T',1}}));

% R18: K_MKP3_T <=> K + MKP3, h6 = 0.086, h-6 = 0.0011
MAPK = MAPK.addReaction(struct('propensity',{'h6*K_MKP3_T'},'stoichiometry',{{'K_MKP3_T',-1;'K',1;'MKP3',1}}));
MAPK = MAPK.addReaction(struct('propensity',{'h_6*K*MKP3'},'stoichiometry',{{'K',-1;'MKP3',-1;'K_MKP3_T',1}}));

% R19: KpY + MKP3 <=> KpY_MKP3, h7 = 0.01, h-7 = 1.0
MAPK = MAPK.addReaction(struct('propensity',{'h7*KpY*MKP3'},'stoichiometry',{{'KpY',-1;'MKP3',-1;'KpY_MKP3',1}}));
MAPK = MAPK.addReaction(struct('propensity',{'h_7*KpY_MKP3'},'stoichiometry',{{'KpY_MKP3',-1;'KpY',1;'MKP3',1}}));

% R20: KpY_MKP3 -> K_MKP3_Y, h8 = 0.47
MAPK = MAPK.addReaction(struct('propensity',{'h8*KpY_MKP3'},'stoichiometry',{{'KpY_MKP3',-1;'K_MKP3_Y',1}}));

% R21: K_MKP3_Y <=> K + MKP3, h9 = 0.14, h-9 = 0.0018
MAPK = MAPK.addReaction(struct('propensity',{'h9*K_MKP3_Y'},'stoichiometry',{{'K_MKP3_Y',-1;'K',1;'MKP3',1}}));
MAPK = MAPK.addReaction(struct('propensity',{'h_9*K*MKP3'},'stoichiometry',{{'K',-1;'MKP3',-1;'K_MKP3_Y',1}}));

MAPK.tSpan = linspace(0,10000,101);

MAPK.summarizeModel;

%% Form FSP Codes (Can Skip this step for SSA Only)
MAPK = MAPK.formPropensitiesGeneral('MAPK');
%% Solve MAPK Model using FSP Approach
MAPK.fspOptions.verbose = true;
MAPK.customConstraintFuns = ...
{'MKP3+Kpp_MKP3+KpT_MKP3_Y+KpT_MKP3_T+K_MKP3_T+KpY_MKP3+K_MKP3_Y';  % Sum of MKP3
    'K+K_MEK_Y+K_MEK_T+K_MKP3_T+K_MKP3_Y';  % Sum of K
    'Kpp+Kpp_MKP3'; % Sum of Kpp
    'KpY+KpY_MEK+KpY_MKP3'; % Sum of KpY
    'KpT+KpT_MEK+KpT_MKP3_Y+KpT_MKP3_T'; % Sum of KpT
    'K+K_MEK_Y+K_MEK_T+K_MKP3_T+K_MKP3_Y+KpY+KpY_MEK+KpY_MKP3+KpT+KpT_MEK+KpT_MKP3_Y+KpT_MKP3_T'}; % Sum of all

MAPK.fspOptions.fspTol = 0.001;
tic
[~,~,MAPK] = MAPK.solve;
toc
MAPK.plotFSP(speciesNames={'MKP3','K','Kpp'})

%% Solve MAPK Model using SSA Approach
MAPKSSA = MAPK; 
MAPKSSA.Solutions = [];
% Run once to write C Codes.
MAPKSSA.solutionScheme = 'ssa';
MAPKSSA.ssaOptions.Nsims = 10;
[~,~,MAPKSSA] = MAPKSSA.solve;

%% Now run for ful set of simulations.
MAPKSSA.ssaOptions.Nsims = 1e4;
MAPKSSA.ssaOptions.useParallel = false;
tic
[~,~,MAPKSSA] = MAPKSSA.solve;
toc
hold on
plot(MAPKSSA.tSpan,mean(MAPKSSA.Solutions.trajs(1,:,:),3),'c--','linewidth',1)
plot(MAPKSSA.tSpan,mean(MAPKSSA.Solutions.trajs(14,:,:),3),'r--','linewidth',1)
plot(MAPKSSA.tSpan,mean(MAPKSSA.Solutions.trajs(16,:,:),3),'b--','linewidth',1)

legend({'MKP3 - FSP','K - FSP','Kpp - FSP','MKP3 - SSA','K - SSA','Kpp - SSA',})

%% 8 Compartment Spatial Schlogel Model
clear

Schlogle = SSIT();
N = 6;

Schlogle.parameters = {
'c1A',2.676;
'c2',0.040;
'c3B', 108.102;
'c4',37.881;
'd',8.2207};

Schlogle.species =   {};
Schlogle.propensityFunctions = {};
Schlogle.stoichiometry = [];
for i = 1:N
    Schlogle.species =   [Schlogle.species,['X',num2str(i)]];
end

Schlogle.initialCondition = 50*ones(N,1);

% Create Reactions
for i = 1:N
    Schlogle = Schlogle.addReaction(struct(...
        'propensity',{['c1A*X',num2str(i),'*(X',num2str(i),'-1)']},...
        'stoichiometry',{{['X',num2str(i)],1}}));
    Schlogle = Schlogle.addReaction(struct(...
        'propensity',{['c2*X',num2str(i),'*(X',num2str(i),'-1)*(X',num2str(i),'-2)']},...
        'stoichiometry',{{['X',num2str(i)],-1}}));
    Schlogle = Schlogle.addReaction(struct(...
        'propensity',{'c3B'},...
        'stoichiometry',{{['X',num2str(i)],1}}));
    Schlogle = Schlogle.addReaction(struct(...
        'propensity',{['c4*X',num2str(i)]},...
        'stoichiometry',{{['X',num2str(i)],-1}}));
    if i>1
        Schlogle = Schlogle.addReaction(struct(...
            'propensity',{['d*X',num2str(i)]},...
            'stoichiometry',{{['X',num2str(i)],-1;['X',num2str(i-1)],1}}));
    end
    if i<N
        Schlogle = Schlogle.addReaction(struct(...
            'propensity',{['d*X',num2str(i)]},...
            'stoichiometry',{{['X',num2str(i)],-1;['X',num2str(i+1)],1}}));
    end
end

%% Solve using SSA
Schlogle.solutionScheme = 'SSA';
Schlogle.ssaOptions.useParallel = false;
Schlogle.tSpan = linspace(0, 1, 101);
Schlogle.ssaOptions.Nsims = 1000000;
tic
[~, ~, Schlogle] = Schlogle.solve;
toc

close all
Schlogle.plotSSA(speciesNames=Schlogle.species,Ylim=[40 50]);

