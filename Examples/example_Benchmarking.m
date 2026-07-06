% Benchmark Examples

for modset = 1:7
    switch modset
        case 1
            Models = {'Pap','Goutsias','BirthDeath','Toggle2','Phage','MichaelisMenten'};
        case 2
            Models = {'Toggle'};
        case 3
            Models = {'TripleRepressor_SS'};
        case 4
            Models = {'GeneExpression_SS'};
        case 5
            Models = {'Goutsias_1000'};
        case 6
            Models = {'MAPK'};
        case 7
            Models = {'Wang2StateNC'};
    end

    %% Expected Time (MacBook Air, Apple M4, 24Gb, MATLAB R2025b)
    %%  * (MackBook Pro, Apple M3 Pro, 18Gb, MATLAB R2025b)
    % NAME      Model Time      Computation Time (Repeat Time)

    % FAST MODELS (<1 m)
    % Pap:      (tf=10)         ~ 1.27 (0.027) s
    %           (tf=100)        ~ 0.46 (0.10) s
    % Goutsias: (tf=100)        ~ 2.04 (0.40) s
    %           (tf=300)        ~ 31.5 (36.9) s
    % BirthDeath: (tf=10)       ~ 0.79  (0.23) s
    % Toggle2:  (tf=30)         ~ 1.9  (0.38) s
    % Phage:    (tf=30)         ~ 4.9  (0.21) s
    % MichaelisMenten: (tf=70)  ~ 5.1  (1.2) s
    % * MAPK:     (tf=10)         ~ 43.3 (1.03) s

    % MODERATE MODELS (1 min to 30 min)
    % Toggle:   (tf=30)         ~ 15.3 (13.9) s
    %           (tf=100)        ~ 43.4 (39.1) s
    % GeneExpression: (tf=SS)   ~ 68.1 (72.4) s

    % SLOW MODELS (> 5 min)
    % NAME                      Model Time  Computation Time (Repeat Time)
    % TripleRepressor_SS:       (tf=SS)     ~ 598 (148) s
    % Goutsias_1000:            (tf=1000)   ~ 6107 ()

    %%
    clear benchmarks
    close all
    for iM = 1:length(Models)
        clear Model verificationCode timeSets followUp
        close all
        [Model,verificationCode,timeSets,followUp] = Generate_Model_from_Benchmark_Library(Models{iM});
        Model.propensityFilePrefix = Models{iM};
        for iT = 1:length(timeSets)
            if timeSets(iT)>0
                Model.tSpan = linspace(0,timeSets(iT),length(Model.tSpan));
            else
                Model.tSpan = 0;
            end

            disp('*******************************************************')
            disp(['Running benchmarks for model: ', Models{iM}, ' with time set: ', num2str(timeSets(iT))]);
            benchmarks.([Models{iM},'_',num2str(timeSets(iT))]) = ...
                run_benchmarks(Model,verbose=false, ...
                ssaInitialize=false, ...
                addCustomConstraints=true, ...
                followUp = followUp);
            disp('Benchmark Complete')
            benchmarks.([Models{iM},'_',num2str(timeSets(iT))])
        end
        save(['BenchmarkResuts',num2str(modset)],"benchmarks");
    end
end

function [Model,verificationCode,timeSets,followUp] = Generate_Model_from_Benchmark_Library(Name)
arguments
    Name
end

verificationCode =[];
followUp =[];

switch Name
    case 'Toggle'
        Model = SSIT('Empty');
        Model.parameters = {'d1',1;'a1',5000;'b',2.5;'d2',1;'a2',1600;'g',1.5};
        Model.species = {'U';'V'};
        Model.stoichiometry = [-1,1,0, 0;...
            0, 0,-1,1];
        Model.propensityFunctions = {'d1*U';'a1/(1+V^b)';'d2*V';'a2/(1+U^g)'};
        Model.initialCondition = [100;0];
        timeSets = [30,100];
        verificationCode = "Model.plotFSP(plotType = 'marginals',indTimes=length(Model.tSpan))";

    case 'Pap'
        Model = SSIT('Empty');
        Model.parameters = {'LRPtot',100;'k1',1;'k2a',0.25;'k2b',2.25;...
            'k3',1;'k4a',1;'k4b',0.2;'k5',0.01;'k6a',1;'k6b',0.2;...
            'k7',0.01;'k8a',0.25;'k8b',2.25;'k9',10;'k10',1};
        Model.species = {'G1';'G2';'G3';'G4';'LRP';'PapI'};
        % Species order:
        % G1   = DNA
        % G2   = LRP.DNA
        % G3   = DNA.LRP
        % G4   = LRP.DNA.LRP
        % LRP  = free LRP
        % PapI = r
        Model.stoichiometry = [
            -1,  1, -1,  1,  0,  0,  0,  0,  0,  0
            1, -1,  0,  0, -1,  1,  0,  0,  0,  0
            0,  0,  1, -1,  0,  0, -1,  1,  0,  0
            0,  0,  0,  0,  1, -1,  1, -1,  0,  0
            -1,  1, -1,  1, -1,  1, -1,  1,  0,  0
            0,  0,  0,  0,  0,  0,  0,  0,  1, -1
            ];
        Model.propensityFunctions = {'k1*LRPtot*G1';...
            '(k2a + k2b/(1 + PapI))*G2';...
            'k3*LRPtot*G1';...
            '(k4a + k4b/(1 + PapI))*G3';...
            'k5*(LRPtot - 1)*G2';...
            '(k6a + k6b/(1 + PapI))*G4';...
            'k7*(LRPtot - 1)*G3';...
            '(k8a + k8b/(1 + PapI))*G4';...
            'k9*G2';...
            'k10*PapI'};
        % Physically valid starting point:
        % one free DNA template, two free LRP molecules, zero PapI.
        Model.initialCondition = [1;0;0;0;2;0];
        Model.tSpan = linspace(0,10,11);
        timeSets = [10,100];
        %Model.fspTol = 1e-5;
        verificationCode = "Model.plotFSP(plotType = 'marginals',indTimes=length(Model.tSpan),speciesNames='PapI')";

    case 'Goutsias'
        Model = SSIT('Empty');
        Model.parameters = {'k1',0.043;'k2',0.0007;'k3',0.0715;'k4',0.0039;'k5',0.012e9/(6.0221415e23*1e-15);...
            'k6',0.4791;'k7',0.00012e9/(6.0221415e23*1e-15);'k8',0.8765e-11;'k9',0.05e9/(6.0221415e23*1e-15);...
            'k10',0.5};
        Model.species = {'D','DNA','M','RNA','DNAD','DNA2D'};
        Model.stoichiometry = [0, 0,0, 0,-1, 1,-1, 1, 1, -1;...
            0, 0,0, 0,-1, 1, 0, 0, 0, 0;...
            1,-1,0, 0, 0, 0, 0, 0,-2, 2;...
            0, 0,1,-1, 0, 0, 0, 0, 0, 0;...
            0, 0,0, 0, 1,-1,-1, 1, 0, 0;...
            0, 0,0, 0, 0, 0, 1,-1, 0, 0];
        Model.initialCondition = [6;2;2;0;0;0];

        % Reorder species DNA-RNA-Protein
        Model.species = Model.species([2,5,6,4,3,1]);
        Model.initialCondition = Model.initialCondition([2,5,6,4,3,1]);
        Model.stoichiometry = Model.stoichiometry([2,5,6,4,3,1],:);

        Model.propensityFunctions = {'k1*RNA';'k2*M';'k3*DNAD';'k4*RNA';'k5*DNA*D';...
            'k6*DNAD';'k7*DNAD*D';'k8*DNA2D';'(k9/2)*M*(M-1)';'k10*D'};
        Model.tSpan = linspace(0,10,11);
        verificationCode = "Model.plotFSP(plotType = 'marginals',indTimes=length(Model.tSpan),speciesNames='D')";
        timeSets = [100,300];
    case 'BirthDeath'
        Model = SSIT('Empty');
        Model.parameters = {'bk',1000;'dk',1};
        Model.species = {'X'};
        Model.stoichiometry = [1,-1];
        Model.propensityFunctions = {'bk';'dk*X'};
        Model.initialCondition = 0;
        timeSets = 10;
        verificationCode = "Model.plotFSP(plotType = 'marginals',indTimes=length(Model.tSpan),speciesNames='X')";
        %Model.fspTol = 1e-6;
    case 'MichaelisMenten'
        Model = SSIT('Empty');
        Model.parameters = {'k1',1;'k2',1;'k3',0.1};
        Model.species = {'S';'E';'ES';'P'};
        Model.stoichiometry = [
            -1,  1,  0
            -1,  1,  1
            1, -1, -1
            0,  0,  1
            ];
        Model.propensityFunctions = {'k1*S*E';'k2*ES';'k3*ES'};
        Model.initialCondition = [100;1000;0;0];
        timeSets = 70;
        verificationCode = "Model.plotFSP(plotType = 'marginals',indTimes=length(Model.tSpan))";
    case 'Toggle2'
        Model = SSIT('Empty');
        Model.parameters = {'k1',40;'k2',20;'k3',1;'k4',1;'k5',1e-5;...
            'k6',3.5e-5;'k7',1;'k8',1};
        Model.species = {'GeneA';'GeneB';'A';'B';'bGeneB';'bGeneA'};
        Model.stoichiometry = [
            0,  0,  0,  0,  0, -1,  0,  1
            0,  0,  0,  0, -1,  0,  1,  0
            1,  0, -1,  0, -2,  0,  2,  0
            0,  1,  0, -1,  0, -2,  0,  2
            0,  0,  0,  0,  1,  0, -1,  0
            0,  0,  0,  0,  0,  1,  0, -1
            ];
        Model.propensityFunctions = {'k1*GeneA';'k2*GeneB';'k3*A';'k4*B';...
            '(k5/2)*A*(A-1)*GeneB';'(k6/2)*B*(B-1)*GeneA';'k7*bGeneB';'k8*bGeneA'};
        Model.initialCondition = [1;1;0;0;0;0];
        timeSets = 30;
        verificationCode = "Model.plotFSP(plotType = 'marginals',indTimes=length(Model.tSpan),speciesNames={'A','B'})";
    case 'Phage'
        Model = SSIT('Empty');
        Model.parameters = {'k1',0.0069;'k2',0.0069;'k3',0.069;'k4',0.0929;'k5',0.0026;...
            'k6',0.0025;'k7',0.021;'k8',0.021;'k9',0.021;'k10',0.021;'k11',0.021;...
            'k12',0.021;'k13',0.00898;'k14',0.00011;'k15',0.01242;'k16',0.00011;...
            'k17',0.00898;'k18',0.2297;'k19',0.0029;'k20',0.0021;'k21',0.0029;...
            'k22',0.2297;'k23',0.2297;'k24',0.2297;'k25',0.0029;'k26',0.0021;...
            'k27',1.13;'k28',0.0106;'k29',0.0106;'k30',0.0106;'k31',1.13;'k32',0.0202;...
            'k33',0.0202;'k34',0.0040;'k35',0.0040;'k36',0.0040;'k37',0.1413;...
            'k38',0.1413;'k39',0.1413;'k40',0.1413;'k41',0.0279;'k42',0.053;...
            'k43',0.0328;'k44',0.053;'k45',0.0279;'k46',0.0022;'k47',0.0022;...
            'k48',0.0008;'k49',0.0008;'k50',0.003};
        Model.species = {'OR1';'OR2';'OR3';'COR1';'COR2';'COR3';'ROR1';'ROR2';'ROR3';...
            'CI2';'Cro2'};
        propensities = {
            'k1*OR3*OR2'
            'k2*OR3*COR2'
            'k3*OR3*ROR2'
            'k4*OR1*OR2'
            'k5*CI2'
            'k6*Cro2'

            'k7*CI2*OR1'
            'k8*CI2*OR2'
            'k9*CI2*OR3'
            'k10*Cro2*OR1'
            'k11*Cro2*OR2'
            'k12*Cro2*OR3'

            'k13*ROR1*OR2'
            'k14*ROR1*ROR2*OR3'
            'k15*ROR1*ROR2*ROR3'
            'k16*ROR1*ROR2*COR3'
            'k17*ROR1*COR2'

            'k18*ROR2*OR1*OR3'
            'k19*ROR2*ROR1*OR3'
            'k20*ROR2*OR1*ROR3'
            'k21*ROR2*ROR1*ROR3'
            'k22*ROR2*COR1*OR3'
            'k23*ROR2*OR1*COR3'
            'k24*ROR2*COR1*COR3'
            'k25*ROR2*ROR1*COR3'
            'k26*ROR2*COR1*ROR3'

            'k27*ROR3*OR2'
            'k28*ROR3*ROR2*OR1'
            'k29*ROR3*ROR2*ROR1'
            'k30*ROR3*ROR2*COR1'
            'k31*ROR3*COR2'

            'k32*COR1*OR2'
            'k33*COR1*ROR2'
            'k34*COR1*COR2*OR3'
            'k35*COR1*COR2*ROR3'
            'k36*COR1*COR2*COR3'

            'k37*COR2*OR1*OR3'
            'k38*COR2*ROR1*OR3'
            'k39*COR2*OR1*ROR3'
            'k40*COR2*ROR1*ROR3'
            'k41*COR2*COR1*OR3'
            'k42*COR2*OR1*COR3'
            'k43*COR2*COR1*COR3'
            'k44*COR2*ROR1*COR3'
            'k45*COR2*COR1*ROR3'

            'k46*COR3*OR2'
            'k47*COR3*ROR2'
            'k48*COR3*COR2*OR1'
            'k49*COR3*COR2*ROR1'
            'k50*COR3*COR2*COR1'
            };
        stoichiometries = {
            {'CI2',1}
            {'CI2',1}
            {'CI2',1}
            {'Cro2',1}
            {'CI2',-1}
            {'Cro2',-1}
            {'CI2',-1;'OR1',-1;'ROR1',1}
            {'CI2',-1;'OR2',-1;'ROR2',1}
            {'CI2',-1;'OR3',-1;'ROR3',1}
            {'Cro2',-1;'OR1',-1;'COR1',1}
            {'Cro2',-1;'OR2',-1;'COR2',1}
            {'Cro2',-1;'OR3',-1;'COR3',1}
            {'ROR1',-1;'CI2',1;'OR1',1}
            {'ROR1',-1;'CI2',1;'OR1',1}
            {'ROR1',-1;'CI2',1;'OR1',1}
            {'ROR1',-1;'CI2',1;'OR1',1}
            {'ROR1',-1;'CI2',1;'OR1',1}
            {'ROR2',-1;'CI2',1;'OR2',1}
            {'ROR2',-1;'CI2',1;'OR2',1}
            {'ROR2',-1;'CI2',1;'OR2',1}
            {'ROR2',-1;'CI2',1;'OR2',1}
            {'ROR2',-1;'CI2',1;'OR2',1}
            {'ROR2',-1;'CI2',1;'OR2',1}
            {'ROR2',-1;'CI2',1;'OR2',1}
            {'ROR2',-1;'CI2',1;'OR2',1}
            {'ROR2',-1;'CI2',1;'OR2',1}
            {'ROR3',-1;'CI2',1;'OR3',1}
            {'ROR3',-1;'CI2',1;'OR3',1}
            {'ROR3',-1;'CI2',1;'OR3',1}
            {'ROR3',-1;'CI2',1;'OR3',1}
            {'ROR3',-1;'CI2',1;'OR3',1}
            {'COR1',-1;'Cro2',1;'OR1',1}
            {'COR1',-1;'Cro2',1;'OR1',1}
            {'COR1',-1;'Cro2',1;'OR1',1}
            {'COR1',-1;'Cro2',1;'OR1',1}
            {'COR1',-1;'Cro2',1;'OR1',1}
            {'COR2',-1;'Cro2',1;'OR2',1}
            {'COR2',-1;'Cro2',1;'OR2',1}
            {'COR2',-1;'Cro2',1;'OR2',1}
            {'COR2',-1;'Cro2',1;'OR2',1}
            {'COR2',-1;'Cro2',1;'OR2',1}
            {'COR2',-1;'Cro2',1;'OR2',1}
            {'COR2',-1;'Cro2',1;'OR2',1}
            {'COR2',-1;'Cro2',1;'OR2',1}
            {'COR2',-1;'Cro2',1;'OR2',1}
            {'COR3',-1;'Cro2',1;'OR3',1}
            {'COR3',-1;'Cro2',1;'OR3',1}
            {'COR3',-1;'Cro2',1;'OR3',1}
            {'COR3',-1;'Cro2',1;'OR3',1}
            {'COR3',-1;'Cro2',1;'OR3',1}
            };

        % Add reactions to model
        for i = 1:numel(propensities)
            newReaction = struct();
            newReaction.propensity = propensities{i};
            newReaction.stoichiometry = stoichiometries{i};

            Model = Model.addReaction(newReaction);
        end
        Model.initialCondition = [1;1;1;0;0;0;0;0;0;0;0];
        timeSets = 30;
        verificationCode = "Model.plotFSP(plotType = 'marginals',indTimes=length(Model.tSpan),speciesNames={'Cro2','CI2'})";

    case 'Goutsias_1000'
        [Model,verificationCode] = Generate_Model_from_Benchmark_Library('Goutsias');
        Model.tSpan = linspace(0,1000,101);

    case 'TripleRepressor_SS'
        Model = SSIT('Empty');
        Model.parameters = {'kOnA0',10;'kOnA_PA',1.5;'kOffA0',7;...
            'kOffA_PC',2;'kOnB0',9;'kOnB_PB',4;'kOffB0',10;'kOffB_PA',4;...
            'kOnC0',11;'kOnC_PC',1.5;'kOffC0',9;'kOffC_PB',2;'kTxA',1.5;...
            'kTxB',1;'kTxC',1.1;'kDegMA',0.5;'kDegMB',0.3;...
            'kDegMC',0.425;'kTlA',9.5;'kTlB',11;'kTlC',10;'kDegPA',14.5;...
            'kDegPB',15;'kDegPC',11};
        Model.species = {'G1A';'G1B';'G1C';'MA';'MB';'MC';'PA';'PB';'PC'};
        Model.stoichiometry = [
            1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
            0,  0,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
            0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
            0,  0,  0,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0
            0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0
            0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0
            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -1,  0,  0
            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -1,  0
            0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -1
            ];
        Model.propensityFunctions = {'(kOnA0 + kOnA_PA*PA)*(1 - G1A)';...
            '(kOffA0 + kOffA_PC*PC)*G1A';'(kOnB0 + kOnB_PB*PB)*(1 - G1B)';...
            '(kOffB0 + kOffB_PA*PA)*G1B';'(kOnC0 + kOnC_PC*PC)*(1 - G1C)';...
            '(kOffC0 + kOffC_PB*PB)*G1C';'kTxA*G1A';'kTxB*G1B';'kTxC*G1C';...
            'kDegMA*MA';'kDegMB*MB';'kDegMC*MC';'kTlA*MA';'kTlB*MB';...
            'kTlC*MC';'kDegPA*PA';'kDegPB*PB';'kDegPC*PC'};
        Model.initialCondition = [0;0;0;0;0;0;0;0;0];
        Model.tSpan = 0;
        Model.fspOptions.initApproxSS = true;
        Model.fspOptions.minSSEscapeRate = 3.6e-3;

        verificationCode = "Model.plotFSP(plotType='marginals')";

    case 'GeneExpression_SS'
        Model = SSIT('Empty');
        Model.parameters = {'t1',50;'t2',4;'t3',0.5;'t4',0.2};
        Model.species = {'M';'P'};
        Model.stoichiometry = [1,0,-1,0;...
            0,1,0,-1];
        Model.propensityFunctions = {'t1';'t2*M';'t3*M';'t4*P'};
        Model.initialCondition = [100;2000];
        Model.tSpan = 0;
        Model.fspOptions.initApproxSS = true;
        Model.fspOptions.minSSEscapeRate = 1e-4;
        verificationCode = "Model.plotFSP(plotType = 'marginals',indTimes=length(Model.tSpan))";
    case 'MAPK'
        Model = SSIT();

        % Species
        Model.species =   {'MKP3'
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

        Model.initialCondition = zeros(length(Model.species),1);
        Model.initialCondition(find(strcmp(Model.species,'MKP3'))) = 1; % MKP3 IC is 1
        Model.initialCondition(find(strcmp(Model.species,'K'))) = 3; % K IC is 3

        % Parameters
        Model.parameters = {
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
        Model.propensityFunctions = {};
        Model.stoichiometry = [];

        % R1: synthesis/degradation of MEK
        Model = Model.addReaction(struct('propensity',{'s2'},'stoichiometry',{{'MEK',1}}));
        Model = Model.addReaction(struct('propensity',{'d2*MEK'},'stoichiometry',{{'MEK',-1}}));

        % R2: synthesis/degradation of MEK
        Model = Model.addReaction(struct('propensity',{'s3*Kpp'},'stoichiometry',{{'MEK',1}}));

        % R3: synthesis/degradation of K
        Model = Model.addReaction(struct('propensity',{'s1'},'stoichiometry',{{'K',1}}));
        Model = Model.addReaction(struct('propensity',{'d1*K'},'stoichiometry',{{'K',-1}}));

        % R4:
        Model = Model.addReaction(struct('propensity',{'d1*KpY'},'stoichiometry',{{'KpY',-1}}));
        Model = Model.addReaction(struct('propensity',{'d1*KpT'},'stoichiometry',{{'KpT',-1}}));
        Model = Model.addReaction(struct('propensity',{'d1*Kpp'},'stoichiometry',{{'Kpp',-1}}));

        % R5: K + MEK -> K_MEK_Y
        Model = Model.addReaction(struct('propensity',{'k1*K*MEK'},'stoichiometry',{{'K',-1;'MEK',-1;'K_MEK_Y',1}}));

        % R6: K_MEK_Y -> KpY + MEK, k2 = 0.06
        Model = Model.addReaction(struct('propensity',{'k2*K_MEK_Y'},'stoichiometry',{{'K_MEK_Y',-1;'KpY',1;'MEK',1}}));

        % R7: KpY + MEK <=> KpY_MEK, k3 = 0.375, k-3 = 1.0
        Model = Model.addReaction(struct('propensity',{'k3*KpY*MEK'},'stoichiometry',{{'KpY',-1;'MEK',-1;'KpY_MEK',1}}));
        Model = Model.addReaction(struct('propensity',{'k_3*KpY_MEK'},'stoichiometry',{{'KpY_MEK',-1;'KpY',1;'MEK',1}}));

        % R8: KpY_MEK -> Kpp + MEK, k4 = 4.5
        Model = Model.addReaction(struct('propensity',{'k4*KpY_MEK'},'stoichiometry',{{'KpY_MEK',-1;'Kpp',1;'MEK',1}}));

        % R9: K + MEK <=> K_MEK_T, k5 = 0.375, k-5 = 1.0
        Model = Model.addReaction(struct('propensity',{'k5*K*MEK'},'stoichiometry',{{'K',-1;'MEK',-1;'K_MEK_T',1}}));
        Model = Model.addReaction(struct('propensity',{'k_5*K_MEK_T'},'stoichiometry',{{'K_MEK_T',-1;'K',1;'MEK',1}}));

        % R10: K_MEK_T -> KpT + MEK, k6 = 0.06
        Model = Model.addReaction(struct('propensity',{'k6*K_MEK_T'},'stoichiometry',{{'K_MEK_T',-1;'KpT',1;'MEK',1}}));

        % R11: KpT + MEK <=> KpT_MEK, k7 = 0.375, k-7 = 1.0
        Model = Model.addReaction(struct('propensity',{'k7*KpT*MEK'},'stoichiometry',{{'KpT',-1;'MEK',-1;'KpT_MEK',1}}));
        Model = Model.addReaction(struct('propensity',{'k_7*KpT_MEK'},'stoichiometry',{{'KpT_MEK',-1;'KpT',1;'MEK',1}}));

        % R12: KpT_MEK -> Kpp + MEK, k8 = 4.5
        Model = Model.addReaction(struct('propensity',{'k8*KpT_MEK'},'stoichiometry',{{'KpT_MEK',-1;'Kpp',1;'MEK',1}}));

        % R13: Kpp + MKP3 <=> Kpp_MKP3, h1 = 0.015, h-1 = 1.0
        Model = Model.addReaction(struct('propensity',{'h1*Kpp*MKP3'},'stoichiometry',{{'Kpp',-1;'MKP3',-1;'Kpp_MKP3',1}}));
        Model = Model.addReaction(struct('propensity',{'h_1*Kpp_MKP3'},'stoichiometry',{{'Kpp_MKP3',-1;'Kpp',1;'MKP3',1}}));

        % R14: Kpp_MKP3 -> KpT_MKP3_Y, h2 = 0.032
        Model = Model.addReaction(struct('propensity',{'h2*Kpp_MKP3'},'stoichiometry',{{'Kpp_MKP3',-1;'KpT_MKP3_Y',1}}));

        % R15: KpT_MKP3_Y <=> KpT + MKP3, h3 = 0.31, h-3 = 0.01
        Model = Model.addReaction(struct('propensity',{'h3*KpT_MKP3_Y'},'stoichiometry',{{'KpT_MKP3_Y',-1;'KpT',1;'MKP3',1}}));
        Model = Model.addReaction(struct('propensity',{'h_3*KpT*MKP3'},'stoichiometry',{{'KpT',-1;'MKP3',-1;'KpT_MKP3_Y',1}}));

        % R16: KpT + MKP3 <=> KpT_MKP3_T, h4 = 0.01, h-4 = 1.0
        Model = Model.addReaction(struct('propensity',{'h4*KpT*MKP3'},'stoichiometry',{{'KpT',-1;'MKP3',-1;'KpT_MKP3_T',1}}));
        Model = Model.addReaction(struct('propensity',{'h_4*KpT_MKP3_T'},'stoichiometry',{{'KpT_MKP3_T',-1;'KpT',1;'MKP3',1}}));

        % R17: KpT_MKP3_T -> K_MKP3_T, h5 = 0.5
        Model = Model.addReaction(struct('propensity',{'h5*KpT_MKP3_T'},'stoichiometry',{{'KpT_MKP3_T',-1;'K_MKP3_T',1}}));

        % R18: K_MKP3_T <=> K + MKP3, h6 = 0.086, h-6 = 0.0011
        Model = Model.addReaction(struct('propensity',{'h6*K_MKP3_T'},'stoichiometry',{{'K_MKP3_T',-1;'K',1;'MKP3',1}}));
        Model = Model.addReaction(struct('propensity',{'h_6*K*MKP3'},'stoichiometry',{{'K',-1;'MKP3',-1;'K_MKP3_T',1}}));

        % R19: KpY + MKP3 <=> KpY_MKP3, h7 = 0.01, h-7 = 1.0
        Model = Model.addReaction(struct('propensity',{'h7*KpY*MKP3'},'stoichiometry',{{'KpY',-1;'MKP3',-1;'KpY_MKP3',1}}));
        Model = Model.addReaction(struct('propensity',{'h_7*KpY_MKP3'},'stoichiometry',{{'KpY_MKP3',-1;'KpY',1;'MKP3',1}}));

        % R20: KpY_MKP3 -> K_MKP3_Y, h8 = 0.47
        Model = Model.addReaction(struct('propensity',{'h8*KpY_MKP3'},'stoichiometry',{{'KpY_MKP3',-1;'K_MKP3_Y',1}}));

        % R21: K_MKP3_Y <=> K + MKP3, h9 = 0.14, h-9 = 0.0018
        Model = Model.addReaction(struct('propensity',{'h9*K_MKP3_Y'},'stoichiometry',{{'K_MKP3_Y',-1;'K',1;'MKP3',1}}));
        Model = Model.addReaction(struct('propensity',{'h_9*K*MKP3'},'stoichiometry',{{'K',-1;'MKP3',-1;'K_MKP3_Y',1}}));

        % Add custom constraints based on sum of each monomer.
        Model.customConstraintFuns = ...
            {'MKP3+Kpp_MKP3+KpT_MKP3_Y+KpT_MKP3_T+K_MKP3_T+KpY_MKP3+K_MKP3_Y';  % Sum of MKP3
            'K+K_MEK_Y+K_MEK_T+K_MKP3_T+K_MKP3_Y';  % Sum of K
            'Kpp+Kpp_MKP3'; % Sum of Kpp
            'KpY+KpY_MEK+KpY_MKP3'; % Sum of KpY
            'KpT+KpT_MEK+KpT_MKP3_Y+KpT_MKP3_T'; % Sum of KpT
            'K+K_MEK_Y+K_MEK_T+K_MKP3_T+K_MKP3_Y+KpY+KpY_MEK+KpY_MKP3+KpT+KpT_MEK+KpT_MKP3_Y+KpT_MKP3_T'}; % Sum of all

        Model.tSpan = linspace(0,10,101);
        Model.fspOptions.fspTol = 1e-3;

        %% Models without Complete Parameters/Reactions
        % case 'MAPK2'
        % 	Model = SSIT('Empty');
        %     Nmapk = 50;
        %     Model.parameters = {'a1',1/Nmapk;'d1',150;'k1',150;'a2',1/Nmapk;...
        %         'd2',150;'k2',150;'a3',1/Nmapk;'d3',150;'k3',150;'a4',1/Nmapk;...
        %         'd4',150;'k4',150;'a5',1/Nmapk;'d5',150;'k5',150;'a6',1/Nmapk;...
        %         'd6',150;'k6',150;'a7',1/Nmapk;'d7',150;'k7',150;'a8',1/Nmapk;...
        %         'd8',150;'k8',150;'a9',1/Nmapk;'d9',150;'k9',150;'a10',1/Nmapk;...
        %         'd10',150;'k10', 150};
        %     Model.species = {'E1';'E2';'KKPase';'KPase';'KKK';'KK';'K';...
        %         'KKKp';'KKp';'KKpp';'Kp';'Kpp';'KKK_E1';'KKKp_E2';'KK_KKKp';...
        %         'KKp_KKPase';'KKp_KKKp';'KKpp_KKPase';'KKpp_K';'Kp_KPase';...
        %         'Kp_KKpp';'Kpp_KPase'};
        %
        %     x0 = zeros(numel(Model.species),1);
        %     x0(strcmp(Model.species,'E1'))     = 50;
        %     x0(strcmp(Model.species,'E2'))     = 50;
        %     x0(strcmp(Model.species,'KKPase')) = 50;
        %     x0(strcmp(Model.species,'KPase'))  = 50;
        %     x0(strcmp(Model.species,'KKK'))    = 50;
        %     x0(strcmp(Model.species,'KK'))     = 50;
        %     x0(strcmp(Model.species,'K'))      = 50;
        %     Model.initialCondition = x0;
        %
        %     % MAPK reactions:
        %     propensities = {
        %         'a1*KKK*E1'
        %         'd1*KKK_E1'
        %         'k1*KKK_E1'
        %
        %         'a2*KKKp*E2'
        %         'd2*KKKp_E2'
        %         'k2*KKKp_E2'
        %
        %         'a3*KK*KKKp'
        %         'd3*KK_KKKp'
        %         'k3*KK_KKKp'
        %
        %         'a4*KKp*KKPase'
        %         'd4*KKp_KKPase'
        %         'k4*KKp_KKPase'
        %
        %         'a5*KKp*KKKp'
        %         'd5*KKp_KKKp'
        %         'k5*KKp_KKKp'
        %
        %         'a6*KKpp*KKPase'
        %         'd6*KKpp_KKPase'
        %         'k6*KKpp_KKPase'
        %
        %         'a7*KKpp*K'
        %         'd7*KKpp_K'
        %         'k7*KKpp_K'
        %
        %         'a8*Kp*KPase'
        %         'd8*Kp_KPase'
        %         'k8*Kp_KPase'
        %
        %         'a9*Kp*KKpp'
        %         'd9*Kp_KKpp'
        %         'k9*Kp_KKpp'
        %
        %         'a10*Kpp*KPase'
        %         'd10*Kpp_KPase'
        %         'k10*Kpp_KPase'
        %     };
        %
        %     stoichiometries = {
        %         {'KKK',-1; 'E1',-1; 'KKK_E1',1}
        %         {'KKK_E1',-1; 'KKK',1; 'E1',1}
        %         {'KKK_E1',-1; 'KKKp',1; 'E1',1}
        %
        %         {'KKKp',-1; 'E2',-1; 'KKKp_E2',1}
        %         {'KKKp_E2',-1; 'KKKp',1; 'E2',1}
        %         {'KKKp_E2',-1; 'KKK',1; 'E2',1}
        %
        %         {'KK',-1; 'KKKp',-1; 'KK_KKKp',1}
        %         {'KK_KKKp',-1; 'KK',1; 'KKKp',1}
        %         {'KK_KKKp',-1; 'KKp',1; 'KKKp',1}
        %
        %         {'KKp',-1; 'KKPase',-1; 'KKp_KKPase',1}
        %         {'KKp_KKPase',-1; 'KKp',1; 'KKPase',1}
        %         {'KKp_KKPase',-1; 'KK',1; 'KKPase',1}
        %
        %         {'KKp',-1; 'KKKp',-1; 'KKp_KKKp',1}
        %         {'KKp_KKKp',-1; 'KKp',1; 'KKKp',1}
        %         {'KKp_KKKp',-1; 'KKpp',1; 'KKKp',1}
        %
        %         {'KKpp',-1; 'KKPase',-1; 'KKpp_KKPase',1}
        %         {'KKpp_KKPase',-1; 'KKpp',1; 'KKPase',1}
        %         {'KKpp_KKPase',-1; 'KKp',1; 'KKPase',1}
        %
        %         {'KKpp',-1; 'K',-1; 'KKpp_K',1}
        %         {'KKpp_K',-1; 'KKpp',1; 'K',1}
        %         {'KKpp_K',-1; 'KKpp',1; 'Kp',1}
        %
        %         {'Kp',-1; 'KPase',-1; 'Kp_KPase',1}
        %         {'Kp_KPase',-1; 'Kp',1; 'KPase',1}
        %         {'Kp_KPase',-1; 'K',1; 'KPase',1}
        %
        %         {'Kp',-1; 'KKpp',-1; 'Kp_KKpp',1}
        %         {'Kp_KKpp',-1; 'Kp',1; 'KKpp',1}
        %         {'Kp_KKpp',-1; 'Kpp',1; 'KKpp',1}
        %
        %         {'Kpp',-1; 'KPase',-1; 'Kpp_KPase',1}
        %         {'Kpp_KPase',-1; 'Kpp',1; 'KPase',1}
        %         {'Kpp_KPase',-1; 'Kp',1; 'KPase',1}
        %     };
        %
        %     assert(numel(propensities) == numel(stoichiometries), ...
        %     'Number of propensities and stoichiometries must match.');
        %     assert(numel(propensities) == 30, ...
        %     'Expected 30 MAPK reactions.');
        %
        %     % Add reactions to model
        %     for i = 1:numel(propensities)
        %
        %         newReaction = struct();
        %         newReaction.propensity = propensities{i};
        %         newReaction.stoichiometry = stoichiometries{i};
        %
        %         Model = Model.addReaction(newReaction);
        %
        %     end
        %     Model.tSpan = linspace(0,10,11);
        %     %Model.fspTol = 1e-1;
        % case 'EnzymaticFutile'
        %     Model = SSIT('Empty');
        %     Model.parameters = {'kplus1',40;'kplus2',1e4;'kminus1',200;...
        %         'kminus2',100;'kplus3',1e4;'kminus3',5000};
        %     Model.species = {
        %         'Xstar'     % phosphorylated substrate
        %         'X'         % unphosphorylated substrate
        %         'Efplus'    % free forward enzyme
        %         'Efminus'   % free reverse enzyme
        %         'Ebplus'    % bound forward enzyme complex
        %         'Ebminus'   % bound reverse enzyme complex
        %     };
        %     Model.initialCondition = [90;30;2;2;0;0];
        %     Model.stoichiometry = [
        %          0,  0,  1, -1,  1, -1
        %         -1,  1,  0,  0,  0,  1
        %         -1,  1,  1,  0,  0,  0
        %          0,  0,  0, -1,  1,  1
        %          1, -1, -1,  0,  0,  0
        %          0,  0,  0,  1, -1, -1
        %     ];
        %
        %     Model.propensityFunctions = {
        %         'kplus1*X*Efplus'       % X + Efplus -> Ebplus
        %         'kminus1*Ebplus'        % Ebplus -> X + Efplus
        %         'kplus2*Ebplus'         % Ebplus -> Xstar + Efplus
        %         'kplus3*Xstar*Efminus'  % Xstar + Efminus -> Ebminus
        %         'kminus3*Ebminus'       % Ebminus -> Xstar + Efminus
        %         'kminus2*Ebminus'       % Ebminus -> X + Efminus
        %     };
        %
        %     Model.tSpan = linspace(0,1);
        %     %Model.fspTol = 1e-6;
        % case 'p53'
        % 	Model = SSIT('Empty');
        %     Model.parameters = {'kp',0.5;'k1',9.963e-6;'dp',1.925e-5;'km',1.5e-3;...
        %     	'k2',1.5e-2;'kD',740;'k0',8e-4;'drc',1.444e-4;'kT',1.66e-2;'ki',9e-4;...
        %     	'dmn',1.66e-7;'k3',9.963e-6;'ka',0.5;'da',3.209e-5};
        %     Model.species = {'p53';'MDM2nuc';'ARF';'RNAnuc';'RNAcyt'; 'MDM2cyt';'MDM2nucARF'};
        %     Model.stoichiometry = [1,-1,0,0,0,0,0,0,0,0,0;...
        %                            0,0,0,0,0,0,1,-1,-1,0,0;...
        %        					   0,0,0,0,0,0,0,0,-1,1,-1;...
        %        					   0,0,1,-1,0,0,0,0,0,0,0;...
        %        					   0,0,0,1,-1,0,0,0,0,0,0;...
        %        					   0,0,0,0,0,1,-1,0,0,0,0;...
        %        					   0,0,0,0,0,0,0,0,1,0,0];
        %     Model.propensityFunctions = {'kp';'dp*p53+k1*p53*MDM2nuc';...
        %     	'km+k2*(p53^(1.8)/(kD^(1.8)+p53^(1.8)))';'k0*RNAnuc';'drc*RNAcyt';'kT*RNAcyt';...
        %     	'ki*MDM2cyt';'dmn*MDM2nuc*MDM2nuc';'k3*MDM2nuc*ARF';'ka';'da*ARF'};
        %     Model.initialCondition = [100;100;100;0;0;0;0];
        %     Model.tSpan = linspace(0,1000,101);
        case 'Wang2StateNC'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Wang et al. PRL 2025 Fig. 3a speed-comparison model
        %% 2-state extrinsic-noise-corrected nuclear/cytoplasmic mRNA model
        %%
        %% Original paper model:
        %%   offGene -> onGene                      rate sigma_on
        %%   onGene  -> offGene                     rate sigma_off
        %%   onGene  -> onGene + nucRNA             rate rho * beta
        %%   nucRNA  => cytRNA after fixed delay    delay tau
        %%   cytRNA  -> null                        rate d
        %%
        %% SSIT-compatible approximation below:
        %%   nucRNA -> cytRNA                       rate kexp*nucRNA, kexp = 1/tau
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Model = SSIT;
        
        % Species:
        % offGene + onGene = number of gene copies.
        % Wang et al. use two independent gene copies, so initial total gene copy
        % number is 2.
        Model.species = {'offGene'; 'onGene'; 'nucRNA'; 'cytRNA'};
        
        % Initial condition:
        % For steady-state FSP/MLE this usually should not matter much after
        % sufficient integration. Both copies are initialized inactive here.
        Model.initialCondition = [2; 0; 0; 0];
        
        % Reactions:
        %   1. offGene -> onGene
        %   2. onGene  -> offGene
        %   3. onGene  -> onGene + nucRNA
        %   4. nucRNA  -> cytRNA          % Markovian approximation to fixed delay
        %   5. cytRNA  -> null
        Model.stoichiometry = [...
            -1,  1,  0,  0,  0;  % offGene
             1, -1,  0,  0,  0;  % onGene
             0,  0,  1, -1,  0;  % nucRNA
             0,  0,  0,  1, -1]; % cytRNA
        
        % Propensity functions:
        % beta is the extrinsic-noise multiplier for cell volume.
        % For a single fixed-beta run, set beta = 1.
        % To mimic Wang et al.'s ENC mixture, sample beta externally for each cell
        % or quadrature component and rerun/mixture-average the model.
        Model.propensityFunctions = {...
            'sigma_on * offGene'; ...
            'sigma_off * onGene'; ...
            'rho * beta * onGene'; ...
            'kexp * nucRNA'; ...
            'd * cytRNA'};
        
        % Example parameter values / placeholders.
        % Wang et al. used multiple parameter sets in their supplement; the main
        % paper states d = 1 for the related total-mRNA telegraph comparison.
        % Replace these with the relevant parameter set if you obtain SM Table S2.
        tau = 0.8;
        
        Model.parameters = ({...
            'sigma_on',  0.2; ...
            'sigma_off', 2.0; ...
            'rho',       10; ...
            'kexp',      1/tau; ...
            'd',         1; ...
            'beta',      1});
        
        Model.summarizeModel
end
if ~exist("timeSets","var")
    timeSets = Model.tSpan(end);
end
end

function benchmarks = run_benchmarks(Model,opts)
arguments
    Model
    opts.nSims = 1000;
    opts.runReductions = false;
    opts.verbose = false;
    opts.runParallel = false;
    opts.verificationCode = '';
    opts.ssaInitialize = false;
    opts.addCustomConstraints = false;
    opts.followUp = [];
end

% FSP solutions
Model.solutionScheme = 'fsp';
Model.fspOptions.verbose = opts.verbose;

% Add list of common FSP constraints.
if opts.addCustomConstraints
    Model = Model.addFSPConstraints(anticorrelatedPairs='all',correlatedPairs='all');
end

% Write FSP Propensity function codes
tic
Model = Model.formPropensitiesGeneral(Model.propensityFilePrefix);
disp('FSP propensity function formed:')
benchmarks.writeFSPcodes = toc;

% Run SSA to Initialize FSP projections
if opts.ssaInitialize
    tic
    Model = Model.ssaInitializeConstraints(100);
    disp('SSA initialization solve:')
    benchmarks.SSAinitialization = toc;
end

tic
% Run first FSP Solve
[~,~,Model] = Model.solve;
disp('FSP initial solve:')
if ~isempty(opts.followUp)
    eval(opts.followUp);
end
benchmarks.initialFSPSolve = toc;

% Run another solve using the identified state space.
tic
[~,~,Model] = Model.solve;
disp('FSP subsequent solve:')
if ~isempty(opts.followUp)
    eval(opts.followUp);
end
benchmarks.subsequentFSPSolve = toc;

% Run verification code (to make plots to check results)
if ~isempty(opts.verificationCode)
    eval(opts.verificationCode);
end

% Record size of FSP projection
benchmarks.fspSize = size(Model.Solutions.stateSpace.states,2);
benchmarks.bounds = Model.fspOptions.bounds;

%% SSA Solutions
Model.solutionScheme = 'ssa';

% Run first SSA solution to write analysis codes.
Model.ssaOptions.Nsims = 1;
tic
[~,~,Model] = Model.solve;
benchmarks.initialSSASolve_1run = toc;

% Run list of SSA trajectories.
Model.ssaOptions.Nsims = opts.nSims;
Model.ssaOptions.useParallel = false;
tic
[~,~,Model] = Model.solve;
benchmarks.(['subsequentSSASolve_',num2str(opts.nSims),'runs_serial']) = toc;

% Run again in parallel.
if opts.runParallel
    Model.ssaOptions.Nsims = opts.nSims;
    Model.ssaOptions.useParallel = true;
    tic
    [~,~,Model] = Model.solve;
    benchmarks.(['subsequentSSASolve_',num2str(opts.nSims),'runs_parallel']) = toc;
end

%% ODE Solver
if length(Model.tSpan)>1
    Model.solutionScheme = 'ode';
    tic
    [~,~,Model] = Model.solve;
    benchmarks.initialODEsolve = toc;

    tic
    [~,~,Model] = Model.solve;
    benchmarks.subsequentODEsolve = toc;
end

%% Model Reduction FSP
if opts.runReductions
    Model.solutionScheme = 'fsp';

    Model.tSpan = linspace(min(Model.tSpan),max(Model.tSpan),150);
    [~,~,Model] = Model.solve;

    for redOrder = [20,30,40,50]
        Model2 = Model;
        Model2.modelReductionOptions.useModReduction = true;
        Model2.fspOptions.fspTol = inf;
        Model2.modelReductionOptions.reductionType = 'POD2';

        Model2.modelReductionOptions.reductionOrder = redOrder;

        tic
        Model2 = Model2.computeModelReductionTransformMatrices();
        benchmarks.(['PODModelReductionTime_',num2str(redOrder)]) = toc;

        tic
        [fspSoln2] = Model2.solve();
        benchmarks.(['ReducedModelSolveTime_',num2str(redOrder)]) = toc

        dims = 1:fspSoln2.fsp{end}.p.dim;
        for i = 1:length(dims)
            PODfinalError(i) = sum(sum(abs((double(fspSoln2.fsp{end}.p.data - fspSoln.fsp{end}.p.data))),setdiff(dims,i)));
        end
        benchmarks.(['ReducedModelError_',num2str(redOrder)]) = sum(PODfinalError);
    end
end

% Remove model to return memory
Model = []; Model2 = [];
clear Model Model2

end
