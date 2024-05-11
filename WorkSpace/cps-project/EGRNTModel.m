classdef EGRNTModel < GeneExprSysModel
    methods
        % Need constructor to set model specifics
        function model = EGRNTModel(rng)
            % Call superclass constructor with no arguments
            model@GeneExprSysModel({});
            model.species = {'GRc';'GRn';'on';'off';'mRNA'};
            model.initialCondition = [20;10;0;1;0];
            model.propensityFunctions = {...
                '(kcn0+kcn1*IDex';... % GR translocation to nucleus
                'knc*GRn';... % GR translocation to cytoplasm
                'kon*GRn*on';... % Gene activation (off to on)
                'koff*off';... % Gene deactivation (on to off)
                'kr*on';... % Transcription (production of mRNA)
                'gr*mRNA'... % mRNA degradation
                };
            model.parameters = ({'koff',0.14;'kon',0.01;'kr',1;'gr',0.01;...
                'kcn0',0.01;'kcn1',0.1;'knc',1;'r1',0.01});
            model.inputExpressions = {'IDex','(t>0)*exp(-r1*t)'};
            model.stoichiometry = [...
                -1, 1, 0, 0, 0;... % GR translocation to nucleus
                1, -1, 0, 0, 0;... % GR translocation to cytoplasm
                0, 0,  1,-1, 0;... % Gene activation (off to on)
                0, 0, -1, 1, 0;... % Gene deactivation (on to off)
                0, 0, 0,  0, 1;... % Transcription (production of mRNA)
                0, 0, 0,  0,-1 ... % mRNA degradation
                ];
            model.pdoOptions.unobservedSpecies = {'on','off'};

            %% Prior
            muLog10Prior = [0,0,0,0,0];
            sigLog10Prior = [2 2 2 2 2];
            model.fittingOptions.modelVarsToFit = ...
                [1:length(model.parameters)];
            model.fittingOptions.logPrior = ...
                @(p)-(log10(p(:))-muLog10Prior').^2./(2*sigLog10Prior'.^2);
            model.fittingOptions.logPriorCovariance = ...
                diag(sigLog10Prior.^2*log(10^2));
        end
    end
end