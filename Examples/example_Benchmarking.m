% Benchmark Examples

Models = {'Toggle'};
clear benchmarks
for iM = 1:length(Models)
    Model = Generate_Model_from_Benchmark_Library(Models{iM});
    Model.propensityFilePrefix = Models{iM};
    benchmarks.(Models{iM}) = run_benchmarks(Model);
end

function Model = Generate_Model_from_Benchmark_Library(Name)
arguments
    Name
end

switch Name
    case 'Toggle'
        Model = SSIT;
        Model.parameters = {'d1',1;'a1',5000;'b',2.5;'d2',1;'a2',1600;'g',1.5};
        Model.species = {'U';'V'};
        Model.stoichiometry = [-1,1,0, 0;...
           					    0, 0,-1,1];
        Model.propensityFunctions = {'d1*U';'a1/(1+V^b)';'d2*V';'a2/(1+U^g)'};
        Model.initialCondition = [0;0];
        Model.tSpan = linspace(0,100);
        Model.fspTol = 1e-6;
	case 'Pap'
		Model = SSIT;
        Model.parameters = {'x1',0;'x2',0;'x3',0;'x4',0;'x5',0};
        Model.species = {'G1';'G2';'G3';'G4';'LRP';'PapI'};
        Model.stoichiometry = [-1,1,-1,1,0,0,0,0,0,0;...
           					    1,-1,0,0,-1,1,0,0,0,0;...
           					    0,0,1,-1,0,0,-1,1,0,0;...
           					    0,0,0,0,1,-1,1,-1,0,0;...
           					    -1,1,-1,1,-1,1,-1,1,0,0;...
           					    0,0,0,0,0,0,0,0,1,-1];
        Model.propensityFunctions = {'LRP*x1';'(0.25+2.25/(1+x5))*x2';'LRP*x1';...
        	'(1+0.2/(1+x5))*x3';'0.01*(LRP-1)*x2';'(0.25+2.25/(1+x5))*x4';'10*x2';'x4'};
        Model.initialCondition = [1;1;0;0;0;0];
        Model.tSpan = linspace(0,10);
        Model.fspTol = 1e-6;

        % 1: LRP + G1 -> G2
        % 2: G2 -> LRP + G1
        % 3: LRP + G1 -> G3
        % 4: G3 -> G1 + LRP
        % 5: G2 -> G4
        % 6: G4 -> G2 + LRP
        % 7: G3 -> G4
        % 8: G4 -> G3 + LRP
        % 9: G2 -> G2 + PapI
        % 10: PapI -> 0
	case 'Goutsias'
		Model = SSIT;
        Model.parameters = {'k1',0.043;'k2',0.0007;'k3',0.0715;'k4',0.0039;'k5',0.012e9;...
        	'k6',0.4791;'k7',0.00012e9;'k8',0.8765e-11;'k9',0.05e9;'k10',0.5};
        Model.species = {'RNA','M','DNAD','DNA','D','DNA2D'};
        Model.stoichiometry = [0,0,1,-1,0,0,0,0,0,0;...
        					   1,-1,0,0,0,0,0,0,-2,2;...
        					   0,0,0,0,1,-1,-1,1,0,0;...
        					   0,0,0,0,-1,1,0,0,0,0;...
        					   0,0,0,0,-1,1,-1,1,1,-1;...
        					   0,0,0,0,0,0,1,-1,0,0];
        Model.propensityFunctions = {'k1*RNA';'k2*M';'k3*DNAD';'k4*RNA';'k5*DNA*D';...
        	'k6*DNAD';'k7*DNAD*D';'k8*DNA2D';'(k9/2)*M*(M-1)';'k10*D'};
        Model.initialCondition = [0;2;0;2;6;0];
        Model.tSpan = linspace(0,300);
        Model.fspTol = 1e-6;
	case 'BirthDeath'
		Model = SSIT;
        Model.parameters = {'bk',1000;'dk',1};
        Model.species = {'X'};
        Model.stoichiometry = [1,-1];
        Model.propensityFunctions = {'bk';'dk*X'};
        Model.initialCondition = 0;
        Model.tSpan = linspace(0,10);
        Model.fspTol = 1e-6;
	case 'GeneExpression'
		Model = SSIT;
        Model.parameters = {'t1',50;'t2',4;'t3',0.5;'t4',0.2};
        Model.species = {'M';'P'};
        Model.stoichiometry = [1,0,-1,0;...
           					   0,1,0,-1];
        Model.propensityFunctions = {'t1';'t2*M';'t3*M';'t4'};
        Model.initialCondition = [0;0];
	case 'MichaelisMenten'
	    Model = SSIT;
        Model.parameters = {'k1',1;'k2',1;'k3',0.1};
        Model.species = {'S';'E';'ES'; 'P'};
        Model.stoichiometry = [-1,1,0;...
           					   -1,1,1;...
           					    1,-1,-1;...
           					    0,0,1];
        Model.propensityFunctions = {'k1*S*E';'k2*ES';'k3*ES'};
        Model.initialCondition = [1;1;0;0];
        Model.tSpan = linspace(0,10);
        Model.fspTol = 1e-6;
	case 'MAPK'
		Model = SSIT;
        Model.parameters = {'k1',0;};
        species = {'E1';'E2';'KKPase';'KPase';'KKK';'KKKp';'KK';...
            'KKp';'KKpp';'K';'Kp';'Kpp';'KKK_E1';'KKKp_E2';'KK_KKKp';...
            'KKp_KKPase';'KKp_KKKp';'KKpp_KKPase';'KKpp_K';'Kp_KPase';...
            'Kp_KKpp';'Kpp_KPase'};

        x0 = zeros(numel(species),1);     
        x0(strcmp(species,'E1'))     = 50;
        x0(strcmp(species,'E2'))     = 50;
        x0(strcmp(species,'KKPase')) = 50;
        x0(strcmp(species,'KPase'))  = 50;
        x0(strcmp(species,'KKK'))    = 50;
        x0(strcmp(species,'KK'))     = 50;
        x0(strcmp(species,'K'))      = 50;
        Model.initialCondition = x0;

        % MAPK reactions:
        propensities = {
            'a1*KKK*E1'
            'd1*KKK_E1'
            'k1*KKK_E1'
        
            'a2*KKKp*E2'
            'd2*KKKp_E2'
            'k2*KKKp_E2'
        
            'a3*KK*KKKp'
            'd3*KK_KKKp'
            'k3*KK_KKKp'
        
            'a4*KKp*KKPase'
            'd4*KKp_KKPase'
            'k4*KKp_KKPase'
        
            'a5*KKp*KKKp'
            'd5*KKp_KKKp'
            'k5*KKp_KKKp'
        
            'a6*KKpp*KKPase'
            'd6*KKpp_KKPase'
            'k6*KKpp_KKPase'
        
            'a7*KKpp*K'
            'd7*KKpp_K'
            'k7*KKpp_K'
        
            'a8*Kp*KPase'
            'd8*Kp_KPase'
            'k8*Kp_KPase'
        
            'a9*Kp*KKpp'
            'd9*Kp_KKpp'
            'k9*Kp_KKpp'
        
            'a10*Kpp*KPase'
            'd10*Kpp_KPase'
            'k10*Kpp_KPase'
        };
        
        stoichiometries = {
            {'KKK',-1; 'E1',-1; 'KKK_E1',1}
            {'KKK_E1',-1; 'KKK',1; 'E1',1}
            {'KKK_E1',-1; 'KKKp',1; 'E1',1}
        
            {'KKKp',-1; 'E2',-1; 'KKKp_E2',1}
            {'KKKp_E2',-1; 'KKKp',1; 'E2',1}
            {'KKKp_E2',-1; 'KKK',1; 'E2',1}
        
            {'KK',-1; 'KKKp',-1; 'KK_KKKp',1}
            {'KK_KKKp',-1; 'KK',1; 'KKKp',1}
            {'KK_KKKp',-1; 'KKp',1; 'KKKp',1}
        
            {'KKp',-1; 'KKPase',-1; 'KKp_KKPase',1}
            {'KKp_KKPase',-1; 'KKp',1; 'KKPase',1}
            {'KKp_KKPase',-1; 'KK',1; 'KKPase',1}
        
            {'KKp',-1; 'KKKp',-1; 'KKp_KKKp',1}
            {'KKp_KKKp',-1; 'KKp',1; 'KKKp',1}
            {'KKp_KKKp',-1; 'KKpp',1; 'KKKp',1}
        
            {'KKpp',-1; 'KKPase',-1; 'KKpp_KKPase',1}
            {'KKpp_KKPase',-1; 'KKpp',1; 'KKPase',1}
            {'KKpp_KKPase',-1; 'KKp',1; 'KKPase',1}
        
            {'KKpp',-1; 'K',-1; 'KKpp_K',1}
            {'KKpp_K',-1; 'KKpp',1; 'K',1}
            {'KKpp_K',-1; 'KKpp',1; 'Kp',1}
        
            {'Kp',-1; 'KPase',-1; 'Kp_KPase',1}
            {'Kp_KPase',-1; 'Kp',1; 'KPase',1}
            {'Kp_KPase',-1; 'K',1; 'KPase',1}
        
            {'Kp',-1; 'KKpp',-1; 'Kp_KKpp',1}
            {'Kp_KKpp',-1; 'Kp',1; 'KKpp',1}
            {'Kp_KKpp',-1; 'Kpp',1; 'KKpp',1}
        
            {'Kpp',-1; 'KPase',-1; 'Kpp_KPase',1}
            {'Kpp_KPase',-1; 'Kpp',1; 'KPase',1}
            {'Kpp_KPase',-1; 'Kp',1; 'KPase',1}
        };

        assert(numel(propensities) == numel(stoichiometries), ...
        'Number of propensities and stoichiometries must match.');
        assert(numel(propensities) == 30, ...
        'Expected 30 MAPK reactions.');

        % Add reactions to model
        for i = 1:numel(propensities)
        
            newReaction = struct();
            newReaction.propensity = propensities{i};
            newReaction.stoichiometry = stoichiometries{i};
        
            Model = Model.addReaction(newReaction);
        
        end
        Model.tSpan = linspace(0,10);
        Model.fspTol = 1e-1;
	case 'EnzymaticFutile'
        Model = SSIT;
        Model.parameters = {'kplus1',40;'kplus2',1e4;'kminus1',200;...
            'kminus2',100;'kplus3',1e4;'kminus3',5000};  
        Model.species = {
            'X'         % unphosphorylated substrate
            'Efplus'    % free forward enzyme
            'Ebplus'    % bound forward enzyme complex
            'Xstar'     % phosphorylated substrate
            'Efminus'   % free reverse enzyme
            'Ebminus'   % bound reverse enzyme complex
        };    
        Model.stoichiometry = [
            -1,  1,  0,  0,  0,  1
            -1,  1,  1,  0,  0,  0
             1, -1, -1,  0,  0,  0
             0,  0,  1, -1,  1, -1
             0,  0,  0, -1,  1,  1
             0,  0,  0,  1, -1, -1
        ];  
        Model.propensityFunctions = {
            'kplus1*X*Efplus'       % X + Efplus -> Ebplus
            'kminus1*Ebplus'        % Ebplus -> X + Efplus
            'kplus2*Ebplus'         % Ebplus -> Xstar + Efplus
            'kplus3*Xstar*Efminus'  % Xstar + Efminus -> Ebminus
            'kminus3*Ebminus'       % Ebminus -> Xstar + Efminus
            'kminus2*Ebminus'       % Ebminus -> X + Efminus 
        };  

        Model.initialCondition = [30;2;0;90;2;0];
        Model.tSpan = linspace(0,1);
        Model.fspTol = 1e-6;
    case 'Toggle2'
        Model = SSIT;
        Model.parameters = {'k1',40;'k2',20;'k3',1;'k4',1;'k5',1e-5;...
            'k6',3.5e-5;'k7',1;'k8',1};  
        Model.species = {'GeneA';'A';'GeneB';'B';'bGeneB';'bGeneA'};
        Model.stoichiometry = [
             0,  0,  0,  0,  0, -1,  0,  1
             1,  0, -1,  0, -2,  0,  2,  0
             0,  0,  0,  0, -1,  0,  1,  0
             0,  1,  0, -1,  0, -2,  0,  2
             0,  0,  0,  0,  1,  0, -1,  0
             0,  0,  0,  0,  0,  1,  0, -1
        ];
        Model.propensityFunctions = {'k1*GeneA';'k2*GeneB';'k3*A';'k4*B';...
            '(k5/2)*A*(A-1)*GeneB';'(k6/2)*B*(B-1)*GeneA';'k7*bGeneB';'k8*bGeneA'};  
        Model.initialCondition = [1;0;1;0;0;0];
        Model.tSpan = linspace(0,100);
        Model.fspTol = 1e-6;
	case 'Phage'
		Model = SSIT;
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
        Model.tSpan = linspace(0,100);
        Model.fspTol = 1e-6;
	case 'p53'
		Model = SSIT;
        Model.parameters = {'kp',0.5;'k1',9.963e-6;'dp',1.925e-5;'km',1.5e-3;...
        	'k2',1.5e-2;'kD',740;'k0',8e-4;'drc',1.444e-4;'kT',1.66e-2;'ki',9e-4;...
        	'dmn',1.66e-7;'k3',9.963e-6;'ka',0.5;'da',3.209e-5};
        Model.species = {'p53';'RNAnuc';'RNAcyt'; 'MDM2cyt';'MDM2nuc';'ARF';'MDM2nucARF'};
        Model.stoichiometry = [1,-1,0,0,0,0,0,0,0,0,0;...
           					   0,0,1,-1,0,0,0,0,0,0,0;...
           					   0,0,0,1,-1,0,0,0,0,0,0;...
           					   0,0,0,0,0,1,-1,0,0,0,0;...
           					   0,0,0,0,0,0,1,-1,-1,0,0;...
           					   0,0,0,0,0,0,0,0,-1,1,-1;...
           					   0,0,0,0,0,0,0,0,1,0,0];
        Model.propensityFunctions = {'kp';'dp*p53+k1*p53*MDM2nuc';...
        	'km+k2*(p53^(1.8)/(kD^(1.8)+p53^(1.8)))';'k0*RNAnuc';'drc*RNAcyt';'kT*RNAcyt';...
        	'ki*MDM2cyt';'dmn*MDM2nuc*MDM2nuc';'k3*MDM2nuc*ARF';'ka';'da*ARF'};
        Model.initialCondition = [100;0;0;0;100;100;0];
        Model.tSpan = linspace(0,1000);
        Model.fspTol = 1e-4;
	case 'TripleRepressor'
		Model = SSIT;
        Model.parameters = {};
        Model.species = {'G0A';'G1A';'G0B'; 'G1B';'G0C';'G1C';'MA';'MB';'MC';'PA';'PB';'PC'};
        Model.stoichiometry = [-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
           					    1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
           					    0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
           					    0,0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
           					    0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0;...
           					    0,0,0,0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0;...
           					    0,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,0,0,0;...
           					    0,0,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,0,0;...
           					    0,0,0,0,0,0,0,0,1,0,0,-1,0,0,0,0,0,0;...
           					    0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1,0,0;...
           					    0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1,0;...
           					    0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1];
        Model.propensityFunctions = {'(10+1.5*x7)*(1-x1)';'(7+2*x9)*x1';'(9+4*x8)*(1-x2)';...
        	'(10+4*x7)*x2';'(11+1.5*x9)*(1-x3)';'(9+2*x8)*x3';'1.5*x1';'x2';'1.1*x3';...
        	'0.5*x4';'0.3*x5';'0.425*x6';'9.5*x4';'11*x5';'10*x6';'14.5*x7';'15x8';'11*x9'};
        Model.initialCondition = [1;0;0;0;0;0;0;0;0;0;0;0];
        Model.tSpan;
	end
end

function benchmarks = run_benchmarks(Model,opts)
arguments
    Model
    opts.nSims = 100000;
end
% FSP solutions
Model.solutionScheme = 'fsp';
tic
Model = Model.formPropensitiesGeneral(Model.propensityFilePrefix);
benchmarks.writeFSPcodes = toc;

tic
[~,~,Model] = Model.solve;
benchmarks.initialFSPSolve = toc;

tic
[fspSoln,~,Model] = Model.solve;
benchmarks.subsequentFSPSolve = toc;

% SSA Solutions
Model.solutionScheme = 'ssa';
Model.ssaOptions.Nsims = 1;
tic
[~,~,Model] = Model.solve;
benchmarks.initialSSASolve_1run = toc;

Model.ssaOptions.Nsims = opts.nSims;
Model.ssaOptions.useParallel = false;
tic
[~,~,Model] = Model.solve;
benchmarks.(['subsequentSSASolve_',num2str(opts.nSims),'runs_serial']) = toc;

Model.ssaOptions.Nsims = opts.nSims;
Model.ssaOptions.useParallel = true;
tic
[~,~,Model] = Model.solve;
benchmarks.(['subsequentSSASolve_',num2str(opts.nSims),'runs_parallel']) = toc;

%% ODE Solver
Model.solutionScheme = 'ode';
tic
[~,~,Model] = Model.solve;
benchmarks.initialODEsolve = toc;

tic
[~,~,Model] = Model.solve;
benchmarks.subsequentODEsolve = toc;

%% Model Reduction FSP
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
    benchmarks.(['ReducedModelSolveTime_',num2str(redOrder)]) = toc;

    dims = 1:fspSoln2.fsp{end}.p.dim;
    for i = 1:length(dims)
        PODfinalError(i) = sum(sum(abs((double(fspSoln2.fsp{end}.p.data - fspSoln.fsp{end}.p.data))),setdiff(dims,i)));
    end
    benchmarks.(['ReducedModelError_',num2str(redOrder)]) = sum(PODfinalError);
end

end
