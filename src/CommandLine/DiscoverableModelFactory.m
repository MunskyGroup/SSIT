classdef DiscoverableModelFactory
    methods
        function m = createModel(arch, dataFilename)
            arguments
                arch (1,1) DiscoverableModelArchitecture
                dataFilename (1,1) string {mustBeNonempty}
            end
            m = DiscoverableModel();
            switch arch
                case {"Poisson","PoissonFewer"}
                    m.ArchitectureName = "Poisson";
                    m.species = {'rna'};
                    m.initialCondition = 0;
                    m.propensityFunctions = {'kr*IDex/(kD+IDex)';'gr*rna'};
                    m.stoichiometry = [1,-1];
                    m.TrueParameters = {'kr',10;'gr',0.3;'kD',5};
                    m.parameters = {'kr',10;'gr',0.3;'kD',5};
                    m.DataToFit = {'rna','exp1_s1'};
                    m.FitParameters = 1:3;
                    m.NumberOfTimepoints = 21;
                    m.tSpan = linspace(0,10,nT);
                    m.FIMMetric = FIMMetric();
                    m.FIMMetric.MetricKind = "DeterminantCoveriance";
                    m.NumberOfMHSamples = 1000;

                    m.MuLog10Prior = [0,0,0];
                    m.SigmaLog10Prior = [1 1 1];
                    m.fittingOptions.modelVarsToFit = m.FitParameters;
                    m.fittingOptions.logPrior = ...
                        @(p)-(log10(p(:))-m.MuLog10Prior').^2./(2*m.SigmaLog10Prior'.^2);
                    m.fittingOptions.logPriorCovariance = ...
                        diag(m.SigmaLog10Prior.^2*log(10^2));
                    m.inputExpressions = {'IDex','5'};
                    m = m.formPropensitiesGeneral([dataFilename,'_S',num2str(ind)],true);
                    m = {m};
                % case {"Poisson","PoissonFewer"} ...
                case {"Burst", "BurstFewer"}
                    m.ArchitectureName = "Burst";
                    m.species = {'on';'off';'rna'};
                    m.initialCondition = [0;1;0];
                    m.propensityFunctions = {...
                        'kon*((1-2*atan(alph)/pi) + 2*atan(alph)/pi*((1e-6+IDex)/(M+((1e-6+IDex)))))*off';...
                        'koff*(2*atan(alph)/pi + (1-2*atan(alph)/pi) / (((1e-6+IDex)/(M+((1e-6+IDex))))))*on';...
                        'kr*on';'gr*rna'};
                    m.stoichiometry = [1,-1,0,0;...
                        -1,1,0,0;...
                        0,0,1,-1];
            
                    m.fspOptions.bounds = [0;0;0;1;1;75];
                    
                    m.parameters = {'kon',0.1;'koff',0.2;'kr',10;'gr',0.3;'M',4;'alph',1e-4};
                    m.DataToFit = {'rna','exp1_s3'};
                    m.FitParameters = 1:6;
            
                    m.pdoOptions.unobservedSpecies = {'on','off'};

                    m.NumberOfTimepoints = 31;
                    m.tSpan = linspace(0,30,m.NumberOfTimepoints);
                    m.FIMMetric = FIMMetric();
                    m.FIMMetric.MetricKind = "Determinant";
                    m.FIMMetric.Indices = 1:5;
                    m.NumberOfMHSamples = 5000;
              
                    m.MuLog10Prior = [0,0,0,0,0,0];
                    m.SigmaLog10Prior = [2 2 2 2 2 4];
                    m.fittingOptions.modelVarsToFit = m.FitParameters;
                    m.fittingOptions.logPrior = ...
                        @(p)-(log10(p(:))-m.MuLog10Prior').^2./(2*m.SigmaLog10Prior'.^2);
                    m.fittingOptions.logPriorCovariance = ...
                        diag(m.SigmaLog10Prior.^2*log(10^2));

                    m.inputExpressions = {'IDex','1'};
                    m = m.formPropensitiesGeneral([dataFilename,'_S',num2str(ind)],true);
                    m = {m};
                % case {"Burst", "BurstFewer"} ...
                case {"GR", "GRFewer", "GRFewest"}
                    m.ArchitectureName = "GR";
                    m.species = {'cytGR';'nucGR'};
                    m.initialCondition = [20;1];
                    m.fspOptions.bounds = [0,0,30,30];
                    m.fspOptions.verbose = false;
                    m.fspOptions.fspIntegratorAbsTol = 1e-10;
                    m.sensOptions.solutionMethod = 'finiteDifference';
                    m.propensityFunctions = ...
                        {'kcn0*(1 + (t>0)*kcn1*IDex/(MDex+IDex)) * cytGR';...
                        'knc*nucGR'; 'kg1';'gnuc*nucGR'};
                    m.stoichiometry = [-1,1,1,0;...
                        1,-1,0,-1];
                    m.customConstraintFuns = {'x1+x2'};
            
                    m.parameters = {'kcn0',0.005;'kcn1',0.08;'knc',0.014;...
                        'kg1',0.012;'gnuc',0.005;'MDex',10.44};
                                    
                    m.fspOptions.initApproxSS = true;
                    m.fspOptions.usePiecewiseFSP = true;
                    m.fspOptions.constantJacobian = true;
                    m.fspOptions.constantJacobianTime = 0.1;
            
                    m.DataToFit = {'cytGR','exp1_s1';'nucGR','exp1_s2'};
                    m.FitParameters = 1:6;
                    m.NumberOfTimepoints = 6;
                    m.tSpan = [0,10,30,50,75,180];
                    m.FIMMetric = FIMMetric();
                    m.FIMMetric.MetricKind = "Determinant";
                    m.FIMMetric.Indices = [2, 5, 6];

                    m.NumberOfMHSamples = 5000;
              
                    m.MuLog10Prior = [-2 1 -2 -2 -2 1];
                    m.SigmaLog10Prior = 2*ones(1,6);
                    m.fittingOptions.modelVarsToFit = m.FitParameters;
                    m.fittingOptions.logPrior = ...
                        @(p)-(log10(p(:))-m.MuLog10Prior').^2./(2*m.SigmaLog10Prior'.^2);
                    m.fittingOptions.logPriorCovariance = ...
                        diag(m.SigmaLog10Prior.^2*log(10^2));

                    m.inputExpressions = {'IDex','100'};
                    m = m.formPropensitiesGeneral([dataFilename,'_S',num2str(ind)],true);
                    m = {m};
                % case {"GR", "GRFewer", "GRFewest"} ...
            end % switch arch
        end % createModel
    end % methods
end