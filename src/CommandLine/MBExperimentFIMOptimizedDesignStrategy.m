classdef MBExperimentFIMOptimizedDesignStrategy < ...
        MBAbstractExperimentDesignStrategy

    methods (Access = private)
        function k = findBestMove(obj,Ncp,NcMax)
            arguments
                obj (1, 1) MBExperimentFIMOptimizedDesignStrategy
                Ncp
                NcMax = []              
            end

            info = obj.FIMInfo();            
            FIMs = info.FIMs;
            met = info.Metric;
            statistic = info.Statistic;

            Nt = size(FIMs, 1);
            if isempty(NcMax)
                NcMax = inf * ones(1, Nt);
            end

            Ns = size(FIMs, 2);
            objFun = zeros(Nt, Ns);
            FIM0 = obj.totalFIM(Ncp);

            quantum = obj.Info.ObservationQuantum;
            for is = 1:Ns
                for it = 1:Nt
                    if Ncp(it) + quantum <= NcMax(it)
                        % If one can do that experiment.
                        FIM = FIM0{is} + quantum * FIMs{it, is};
                    else
                        % If there are no more cells available for this
                        % time point.
                        FIM = FIM0{is};
                    end
                    objFun(it, is) = met.calculateForFIM(FIM);
                end
            end

            [~, k] = statistic.calculateStatisticOfObjective(objFun);
        end % findBestMove

        function FIM = totalFIM(obj,Nc)
            arguments 
                obj (1, 1) MBExperimentFIMOptimizedDesignStrategy
                Nc
            end

            info = obj.FIMInfo();            
            FIMs = info.FIMs;
            covPrior = info.CovariancePrior;
            
            if isempty(covPrior)
                fimPrior = zeros(size(FIMs{1}));
            else
                fimPrior = inv(covPrior);
            end

            Nt = size(FIMs, 1);
            Ns = size(FIMs, 2);
            FIM = cell(1, Ns);
            for is = 1:Ns
                FIM{is} = fimPrior;
                for it = 1:Nt
                    FIM{is} = FIM{is} + Nc(it) * FIMs{it, is};
                end
            end
        end % totalFIM
    end % Private methods

    methods (Static, Access = protected)
        function info = createInfoInternal()
            % Descendants can override to use a descendant strategy info.

            info = MBExperimentFIMOptimizedDesignStrategyInfo();
        end
    end % Static protected methods

    methods (Access = protected)
        function [design, cellVecForOptimalFIMCalculation] = ...
                apportionObservationsInternal(obj, design)
            arguments
                obj (1, 1) MBExperimentFIMOptimizedDesignStrategy
                design (1, 1) MBExperimentDesign
            end

            info = obj.FIMInfo();            
            FIMs = info.FIMs;
            
            NcFixed = info.NcFixed;
            NcGuess = info.NcGuess;
            NcMax = design.getTheoreticalObservationMatrix();

            NT = size(FIMs,1);
            NS = size(FIMs,2);

            if isempty(NcFixed)
                NcFixed = zeros(1,NT);
            end

            if isempty(NcMax)
                NcMax = inf*ones(1,NT);
            end

            allCellsHaveBeenDistributed = false;

            if isempty(NcGuess)
                % Distributed available cells among experiments.
                NcGuess = NcFixed;
                iExpt = 1;

                nCellsTotalNew = obj.Info.NumberOfObservationsPerExperiment;

                while nCellsTotalNew>0&&iExpt<=length(NcGuess)
                    avblSlots = NcMax(iExpt) - NcFixed(iExpt);
                    if avblSlots>=nCellsTotalNew
                        NcGuess(iExpt) = NcGuess(iExpt) + nCellsTotalNew;
                        iExpt = inf;
                    else
                        quantum = obj.Info.ObservationQuantum;
                        while avblSlots >= quantum
                            NcGuess(iExpt) = NcGuess(iExpt) + quantum;
                            nCellsTotalNew = nCellsTotalNew - quantum;
                            avblSlots = avblSlots - quantum;
                        end
                        iExpt = iExpt + 1;
                        if iExpt>length(NcGuess)&&nCellsTotalNew>=0
                            NcDNewDesign = NcGuess - NcFixed;
                            warning('All cells have been distributed.')
                            allCellsHaveBeenDistributed = true;
                        end
                    end
                end
            else
                NcGuess = NcFixed+NcGuess;
            end

            if ~allCellsHaveBeenDistributed
                % Process to search for optimal experiment

                converged = false;
                quantum = obj.Info.ObservationQuantum;
                while ~converged
                    converged = true;
                    % Iterate through the time points
                    for i = 1:NT
                        % If
                        while NcGuess(i) - quantum >= NcFixed(i)
                            Ncp = NcGuess;
                            Ncp(i) = Ncp(i) - quantum;

                            k = obj.findBestMove(Ncp, NcMax);
                            if k==i
                                break
                            end
                            NcGuess = Ncp;
                            NcGuess(k) = NcGuess(k) + quantum;
                            converged = false;
                        end
                    end
                end % while ~converged
                NcDNewDesign = NcGuess - NcFixed;
            end

            % Reshape the new design to match the shape expected by the
            % Design class and then update the design accordingly.

            cellVecForOptimalFIMCalculation = NcDNewDesign;
            NcDNewDesign = reshape(NcDNewDesign', [NT, NS])';
            [~, inputConfigurations, times] = ...
                design.getAsObservationMatrix();
            design = design.setFromObservationMatrix(...
                NcDNewDesign, inputConfigurations, times);        
        end % apportionObservationsInternal
        
        function obj = incorporateRoundDetailsInternal(obj, round)
            arguments
                obj (1, 1) MBAbstractExperimentDesignStrategy
                round (1, 1) MBExperimentRound
            end

            info = MBAbstractExperimentFIMOptimizedDesignStrategyInfo();
            info.FIMs = round.FIMResults;
            info.NcGuess = [];
            info.NcFixed = round.CumulativeNumbersOfObservations;
        end % incorporateRoundDetailsInternal
    end % Protected methods

    methods
        function info = FIMInfo(obj)
            arguments (Output)
                info (1, 1) MBExperimentFIMOptimizedDesignStrategyInfo
            end
            % FIMInfo provides a cast of the strategy's info so that all of
            % the properties can be readily obtained.

            %info = MBExperimentFIMOptimizedDesignStrategyInfo(obj.Info); 
            info = obj.Info;
        end
    end % Public methods
end