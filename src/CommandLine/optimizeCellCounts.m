function [NcDNewDesign] = optimizeCellCounts(obj,fims,nCellsTotalNew,FIMMetric,...
                NcGuess,NcFixed,NcMax,statistic,covPrior,incrementAdd)
            %% SSIT.optimizeCellCounts - This function optimizes the number
            %% of cells per time point according to the user-provide metric.
            %
            % Inputs:
            %   * 'fims' - either an [Nt x 1] cell array containing the FIM
            %      matrices for each of the Nt time points, or an [Nt x Ns]
            %      cellarray containing the FIM for each combination of Nt
            %      time points and Ns different parameter sets
            %   * 'nCellsTotalNew' - the total number of cells to be
            %       measured, spread out among the Nt time points
            %   * 'FIMmetric' - type of optimization, allowable metrics are:
            %       'D-opt' - maximize the expected determinant of the FIM
            %       'D-cov' - minimize the expected determinant of the MLE
            %                 covariance
            %       'E-opt' - maximize the smallest eigenvalue of the FIM
            %       'Trace' - maximize the trace of the FIM
            %       'D-opt-sub-inv[<i1>,<i2>,...]' 
            %               - minimize the determinant of the inverse FIM 
            %                 for the specified indices, (all other 
            %                 parameters are assumed to be free)
            %       'D-opt-sub[<i1>,<i2>,...]' 
            %               - maximize the determinant of the FIM for the
            %                 specified indices, (only the parameters in
            %                 obj.fittingOptions.modelVarsToFit are assumed
            %                 to be free)
            %   * 'Nc' - an optimal guess for the optimal experiment
            %            design
            %   * 'NcFixed' - a minimal number of cells to measure at each
            %      time point; this is useful for subsequent experiment
            %      design, having already obtained measured cells from a
            %      previous experiment
            %   * 'NcMax' - maximum total number of cells allowed for each
            %      time point; this is useful in simulated experiment design,
            %      where there are only so many cells available in the real
            %      data
            %
            % Outputs:
            %   * 'Nc' is the optimized experiment design (number of cells
            %      to measure at each point in time)
            %
            % Example: Model.optimizeCellCounts(fimResults,nCellsTotal,...
            %           'D-opt',[],[],[],[],diag(log10.^2));
            arguments
                obj
                fims
                nCellsTotalNew
                FIMMetric = 'E-opt'
                NcGuess = []
                NcFixed = []
                NcMax = []
                statistic = 'mean'
                covPrior = []
                incrementAdd = 1
            end
            if mod(nCellsTotalNew,incrementAdd)~=0
                error('Number of cells must be evenly divisible by incrementAdd.')
            end
            switch FIMMetric
                case 'D-opt'
                    met = @(A)-max(0,det(A));
                case 'D-cov'
                    met = @(A)max(0,det(inv(A)));
                case 'E-opt'
                    met = @(A)-min(eig(A));
                case 'Trace'
                    met = @(A)-trace(A);
                otherwise
                    if startsWith(FIMMetric,'D-opt-sub-inv')
                        k = eval(FIMMetric(14:end));
                        met = @(A)max(0,det(inv(A(k,k))));
                    elseif startsWith(FIMMetric,'D-opt-sub')
                        k = eval(FIMMetric(10:end));
                        met = @(A)-max(0,det((A(k,k))));
                    elseif startsWith(FIMMetric,'D-cov-sub')
                        k = eval(FIMMetric(10:end));
                        ek = zeros(length(k),length(fims{1}));
                        ek(1:length(k),k) = eye(length(k));
                        met = @(A)max(0,det(ek*inv(A)*ek'));
                    else  % all parameters are free.
                        k = eval(FIMMetric);
                        ek = zeros(length(k),length(fims{1}));
                        ek(1:length(k),k) = eye(length(k));
                        met = @(A)max(0,det(ek*inv(A)*ek'));
                    end
            end
            NT = size(fims,1);
            NS = size(fims,2);

            if isempty(NcFixed)
                NcFixed = zeros(1,NT);
            end

            if isempty(NcMax)
                NcMax = inf*ones(1,NT);
            end

            if isempty(NcGuess)
                % Distributed available cells among experiments.
                NcGuess = NcFixed;
                iExpt = 1;
                while nCellsTotalNew>0&&iExpt<=length(NcGuess)
                    avblSlots = NcMax(iExpt) - NcFixed(iExpt);
                    if avblSlots>=nCellsTotalNew
                        NcGuess(iExpt) = NcGuess(iExpt) + nCellsTotalNew;
                        iExpt = inf;
                    else
                        while avblSlots >= incrementAdd
                            NcGuess(iExpt) = NcGuess(iExpt) + incrementAdd;
                            nCellsTotalNew = nCellsTotalNew - incrementAdd;
                            avblSlots = avblSlots - incrementAdd;
                        end
                        iExpt = iExpt + 1;
                        if iExpt>length(NcGuess)&&nCellsTotalNew>=0
                            NcDNewDesign = NcGuess - NcFixed;
                            warning('All cells have been distributed.')
                            return
                        end
                    end
                end
            else
                NcGuess = NcFixed+NcGuess;
            end

            % Process to search for optimal experiment
            Converged = 0;
            while Converged==0
                Converged = 1;
                % Iterate through the time points
                for i = 1:NT
                    % If
                    while NcGuess(i)-incrementAdd>=NcFixed(i)
                        Ncp = NcGuess;
                        Ncp(i) = Ncp(i)-incrementAdd;
                        k = SSIT.findBestMove(fims,Ncp,met,NcMax,statistic,covPrior,incrementAdd);
                        if k==i
                            break
                        end
                        NcGuess = Ncp;
                        NcGuess(k)=NcGuess(k)+incrementAdd;
                        Converged = 0;
                    end
                end
            end
            NcDNewDesign = NcGuess - NcFixed;
        end