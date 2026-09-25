classdef Forman2027

    properties
        kr = 100;
        gr = 1;
        N = 100;
        kD = 10;
        koff_inf = 0.01;
        star_koff = 10;
        mu1 = 2;
        mu2 = 75;
        nCellsFig2 = 10;
        kon_domain = logspace(-2,1,3000);
        Model;
        ModelKoff;
        ModelKoffSig;
        Model_BinomialPDO;
        nCellsInExperiment;
        FIM = [];
        exp1NCells =[];
        exp2NCells =[];
        exp3NCells =[];
        Ncells = 600;
        OptExperiment = [];
        indsFims = [];
        allFims = {};
        tt = [];
        FIM_Exp1
        FIM_Exp2
        FIM_Exp3
        FIM_Opt
        freeParsFig4 = [1:5]

    end

    properties (Dependent)
        star_kon
        final_kon
        final_koff
        Sarray
        star_S
        final_S
    end

    methods
        %% Preliminary Model Definitions.
        function obj = Forman2027()
            %%
            Model = SSIT('Empty');
            Model.species = {'gON','mRNA'};
            Model.initialCondition = [0;0];
            Model.parameters = {'kon0', obj.star_kon;...
                'koff0', obj.star_koff;...
                'kr',obj.kr;...
                'g',obj.gr;...
                'kon1', obj.star_kon; ...
                'koff1', obj.star_koff; ...
                };

            Model.inputExpressions = {'I', ...
                't>=1'};

            Model = Model.addReaction(struct(...
                'propensity',{'(kon0 + (kon1-kon0)*I)*(1-gON)'},...
                'stoichiometry',{{'gON',1}}));

            Model = Model.addReaction(struct(...
                'propensity',{'(koff0 + (koff1-koff0)*I)*gON'},...
                'stoichiometry',{{'gON',-1}}));

            Model = Model.addReaction(struct(...
                'propensity',{'kr*gON'},...
                'stoichiometry',{{'mRNA',1}}));

            Model = Model.addReaction(struct(...
                'propensity',{'g*mRNA'},...
                'stoichiometry',{{'mRNA',-1}}));

            Model.fspOptions.initApproxSS = true;
            Model.tSpan = linspace(0,15,31);

            Model = Model.formPropensitiesGeneral('Forman2027');
            Model = Model.solve;

            obj.Model = Model;

            %%
            % ModelKoff = SSIT('Empty');
            % ModelKoff.species = {'gON','mRNA'};
            % ModelKoff.initialCondition = [0;0];
            % ModelKoff.parameters = {'kon0', obj.star_kon;...
            %     'koff0', obj.star_koff;...
            %     'kr',obj.kr;...
            %     'g',obj.gr;...
            %     'koff1', obj.final_koff; ...
            %     };
            % 
            % ModelKoff.inputExpressions = {'I', ...
            %     't>=1'};
            % 
            % ModelKoff = ModelKoff.addReaction(struct(...
            %     'propensity',{'kon0*(1-gON)'},...
            %     'stoichiometry',{{'gON',1}}));
            % 
            % ModelKoff = ModelKoff.addReaction(struct(...
            %     'propensity',{'(koff0 + (koff1-koff0)*I)*gON'},...
            %     'stoichiometry',{{'gON',-1}}));
            % 
            % ModelKoff = ModelKoff.addReaction(struct(...
            %     'propensity',{'kr*gON'},...
            %     'stoichiometry',{{'mRNA',1}}));
            % 
            % ModelKoff = ModelKoff.addReaction(struct(...
            %     'propensity',{'g*mRNA'},...
            %     'stoichiometry',{{'mRNA',-1}}));
            % 
            % ModelKoff.fspOptions.initApproxSS = true;
            % ModelKoff.tSpan = linspace(0,25,31);
            % 
            % ModelKoff = ModelKoff.formPropensitiesGeneral('Forman2027Koff');
            % ModelKoff = ModelKoff.solve;
            % 
            % obj.ModelKoff = ModelKoff;
            % 
            % obj.nCellsInExperiment = 0*ModelKoff.tSpan;
            % obj.nCellsInExperiment([1,6,31]) = 200;

            %%
            ModelKoffSig = SSIT('Empty');
            ModelKoffSig.species = {'gON','mRNA'};
            ModelKoffSig.initialCondition = [0;0];
            ModelKoffSig.parameters = {'kon', obj.star_kon;...
                'koff_inf', obj.koff_inf;...
                'kr',obj.kr;...
                'g',obj.gr;...
                'kD',obj.kD;...
                'S1', obj.star_S; ...
                'S2', obj.final_S; ...
                };

            ModelKoffSig.inputExpressions = {'I', ...
                'S1 + (S2-S1)*(t>=0)'};

            ModelKoffSig = ModelKoffSig.addReaction(struct(...
                'propensity',{'kon*(1-gON)'},...
                'stoichiometry',{{'gON',1}}));

            ModelKoffSig = ModelKoffSig.addReaction(struct(...
                'propensity',{'koff_inf*(kD+I)/(I)*gON'},...
                'stoichiometry',{{'gON',-1}}));

            ModelKoffSig = ModelKoffSig.addReaction(struct(...
                'propensity',{'kr*gON'},...
                'stoichiometry',{{'mRNA',1}}));

            ModelKoffSig = ModelKoffSig.addReaction(struct(...
                'propensity',{'g*mRNA'},...
                'stoichiometry',{{'mRNA',-1}}));

            ModelKoffSig.fspOptions.initApproxSS = true;
            ModelKoffSig.tSpan = linspace(0,25,31);

            ModelKoffSig = ModelKoffSig.formPropensitiesGeneral('Forman2027KoffSig');
            ModelKoffSig = ModelKoffSig.solve;

            obj.ModelKoffSig = ModelKoffSig;

            obj.nCellsInExperiment = 0*ModelKoffSig.tSpan;
            obj.nCellsInExperiment([1,6,31]) = 200;           
        end
        function Sarray = get.Sarray(obj)
            Sarray = [0.001,0.01,0.1,1,10,100];
        end
        function star_kon = get.star_kon(obj)
            g = (obj.mu1*obj.gr)/obj.kr;
            star_kon = (g*obj.star_koff)/(1-g);
        end
        function final_kon = get.final_kon(obj)
            g = (obj.mu2*obj.gr)/obj.kr;
            final_kon = (g*obj.star_koff)/(1-g);
        end
        function star_S = get.star_S(obj)
            star_S = obj.kD*obj.koff_inf/(obj.star_koff-obj.koff_inf);
        end
        function final_S = get.final_S(obj)
            final_S = obj.kD*obj.koff_inf/(obj.final_koff-obj.koff_inf);
        end
        function final_koff = get.final_koff(obj)
            final_koff = (obj.kr*obj.star_kon)/(obj.mu2*obj.gr)-obj.star_kon;
        end
        %% Figure 1
        function makeFig1B(obj,ax1B)
            arguments
                obj
                ax1B = [];
            end

            if isempty(ax1B)
                figure
                ax1B = gca;
            else
                axes(ax1B);
            end

            kArr = logspace(-2,3,obj.N);
            MEAN = zeros(obj.N,obj.N);
            FANO = zeros(obj.N,obj.N);
            for i = 1:obj.N
                for j = 1:obj.N
                    F(i,j) = kArr(i)/(kArr(i)+kArr(j));
                    MEAN(i,j) = obj.kr/obj.gr*F(i,j);
                    FANO(i,j) = 1+((1-F(i,j))*obj.kr)/(kArr(i)+kArr(j)+obj.gr);
                end
            end
            [~,cF] = contourf(ax1B,log10(kArr),log10(kArr),log10(FANO),linspace(0,log10(64),30)); hold on;
            mu = 2;
            g = mu*obj.gr/obj.kr;
            kons = g*kArr/(1-g);
            plot(log10(kArr), log10(kons), 'k--', 'LineWidth', 2)

            mu = 25;
            g = mu*obj.gr/obj.kr;
            kons = g*kArr/(1-g);
            plot(log10(kArr), log10(kons), 'k--', 'LineWidth', 2)

            mu = 75;
            g = mu*obj.gr/obj.kr;
            kons = g*kArr/(1-g);
            plot(log10(kArr), log10(kons), 'k--', 'LineWidth', 2)

            colormap(ax1B,jet)
            cb = colorbar(ax1B);
            cb.Ticks = log10([1.0001,4,16,64]);
            cb.TickLabels = {'1','4','16','64'};

            set(gca,'xtick',[-2:3],'ytick',[-2:3],...
                'XTickLabel',{'10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}'},...
                'YTickLabel',{'10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}'},...
                'FontSize',16);
            xlim([-2 3]);
            ylim([-2 3]);

            % find star
            % star_koff = 10;
            mu = 2;
            g = (mu*obj.gr)/obj.kr;
            star_kon = (g*obj.star_koff)/(1-g);

            % find collision point with mu == 25 line
            % vary kon
            mu = 75;
            g = (mu*obj.gr)/obj.kr;
            final_kon = (g*obj.star_koff)/(1-g);

            % vary koff
            final_koff = (obj.kr*star_kon)/(mu*obj.gr)-star_kon;

            % plot on figure 1
            % Coordinates
            x0 = log10(obj.star_koff);
            y0 = log10(star_kon);

            x1 = log10(final_koff);
            y1 = log10(star_kon);

            x2 = log10(obj.star_koff);
            y2 = log10(final_kon);


            % Squares
            plot(ax1B, x1, y1, 's', ...
                'MarkerSize', 12, ...
                'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', 'b', ...
                'LineWidth', 3);

            plot(ax1B, x2, y2, 's', ...
                'MarkerSize', 12, ...
                'MarkerEdgeColor', 'k', ...
                'MarkerFaceColor', 'r', ...
                'LineWidth', 3);

            % Arrows
            quiver(ax1B, x0, y0, x1-x0, y1-y0, 0, ...
                Color='b', LineWidth=2, MaxHeadSize=0.25, ShowArrowHead='on');

            quiver(ax1B, x0, y0, x2-x0, y2-y0, 0, ...
                Color='r', LineWidth=2, MaxHeadSize=0.25, ShowArrowHead='on');

            % Star
            plot(ax1B, x0, y0, 'k*', ...
                'MarkerSize', 15, 'LineWidth', 3);
            % plot(ax1, x0, y0, 'k*', ...
            %     'MarkerSize', 14, 'LineWidth', 1.5);
        end
        function makeFigs1CD(obj,ax1C1,ax1C2,ax1C3,ax1D)
            arguments
                obj
                ax1C1
                ax1C2
                ax1C3
                ax1D
            end

            Model_chg = obj.Model;

            Model_chg.plotFSP(plotType='marginals', SpeciesIdx=[2], indTimes=31, figureNums=ax1C1, Title='')

            Model_chg.parameters{5,2} = obj.star_kon;
            Model_chg.parameters{6,2} = obj.final_koff;
            Model_chg = Model_chg.solve;
            Model_chg.plotFSP(plotType='marginals', SpeciesIdx=[2], indTimes=31, figureNums=ax1C2, Title='')
            Model_chg.plotFSP(plotType='meansAndDevs', SpeciesIdx=[2], Colors=[0 0 1], figureNums=ax1D, Title='')

            Model_chg.parameters{6,2} = obj.star_koff;
            Model_chg.parameters{5,2} = obj.final_kon;
            Model_chg = Model_chg.solve;
            Model_chg.plotFSP(plotType='marginals', SpeciesIdx=[2], indTimes=31, figureNums=ax1C3, Title='')
            Model_chg.plotFSP(plotType='meansAndDevs', SpeciesIdx=[2], Colors=[1 0 0], figureNums=ax1D, Title='')

            for figNum = [ax1C1,ax1C2,ax1C3]
                figure(figNum);
                xlim([0 120]);

                ax = gca;
                ax.Children.LineWidth = 4;
                ax.LineWidth = 2;
            end

            figure(ax1D);
            ax = gca;
            set(findobj(ax, '-property', 'LineWidth'), 'LineWidth', 3);
            ax.LineWidth = 2;
        end
        function exportFigs1(~,f1B,f1C1,f1C2,f1C3,f1D)
            %% Export figures for Paper Figure 1
            % Annual Review of Biochemistry
            %
            % Canvas: 6.33 x 7.9 inches
            % Each figure row: 1.975 inches high

            outputFolder = 'AnnualReview_Figures';

            if ~exist(outputFolder, 'dir')
                mkdir(outputFolder);
            end

            fullWidth = 6.33;
            thirdWidth = 6.33 / 3;
            quarterHeight = 7.9 / 4;
            sixtenthHeight = 7.9/16;

            % Figure 1B
            fig = figure(f1B);
            ax = gca;

            % Remove title and axis labels
            ax.Title.String = '';
            ax.XLabel.String = '';
            ax.YLabel.String = '';
            ax.XTickLabel = [];
            ax.YTickLabel = [];
            fig.Units = 'inches';
            fig.Position(3:4) = [fullWidth 6*sixtenthHeight];

            exportgraphics(fig, ...
                fullfile(outputFolder, 'figure1b.svg'), ...
                'ContentType', 'vector');

            % Figure 1C-I, II, III
            figNums = [f1C1,f1C2,f1C3];
            fileNames = {'figure1cI.svg', 'figure1cII.svg', 'figure1cIII.svg'};

            for i = 1:3

                fig = figure(figNums(i));
                ax = gca;

                % Explicitly remove title
                sgtitle(fig, '');

                ax.Title.String = '';
                ax.Title.Visible = 'off';

                ax.Subtitle.String = '';
                ax.Subtitle.Visible = 'off';

                % Explicitly remove axis labels
                ax.XLabel.String = '';
                ax.XLabel.Visible = 'off';

                ax.YLabel.String = '';
                ax.YLabel.Visible = 'off';

                ax.YGrid = 'off';

                ax.XTickLabel = [];
                ax.YTickLabel = [];

                % Keep requested x-axis limits
                xlim(ax, [0 120]);

                % Set physical dimensions
                fig.Units = 'inches';
                fig.Position(3:4) = [thirdWidth 2*sixtenthHeight];

                % Export
                exportgraphics(fig, ...
                    fullfile(outputFolder, fileNames{i}), ...
                    'ContentType', 'vector');

            end

            % Figure 1D

            fig = figure(f1D);
            ax = gca;

            % Remove title and axis labels
            ax.Title.String = '';
            ax.Title.Visible = 'off';

            ax.XLabel.String = '';
            ax.XLabel.Visible = 'off';

            ax.YLabel.String = '';
            ax.YLabel.Visible = 'off';

            ax.XTickLabel = [];
            ax.YTickLabel = [];

            cb = colorbar(ax);
            cb.TickLabels = [];

            % Remove legend
            leg = findobj(fig, 'Type', 'Legend');

            if ~isempty(leg)
                delete(leg);
            end

            % Set physical dimensions
            fig.Units = 'inches';
            fig.Position(3:4) = [fullWidth quarterHeight];

            % Export
            exportgraphics(fig, ...
                fullfile(outputFolder, 'figure1d.svg'), ...
                'ContentType', 'vector');

            disp('All SVG figures exported successfully.');

        end
        
        %% Figure 2
        function obj = makeFigs2A(obj,ax2A, opts)
            arguments
                obj
                ax2A
                opts.xlims = [10^1, 10^2.5]
                opts.ylims = [-10, 0]
            end
            Model_chg = obj.ModelKoffSig;
            Model_chg.fittingOptions.modelVarsToFit = [1];
            nCellsInExperiment = zeros(size(Model_chg.tSpan));
            nCellsInExperiment([1]) = obj.nCellsFig2;
            Model_chg = Model_chg.solve;
            Model_chg.ssaOptions.Nexp = 5000;

            % Model_chg.plotFSP(plotType='meansAndDevs', SpeciesIdx=[2], Title='testing steady state') % Test successful
            Model_chg.sampleDataFromFSP(saveFile='dataForFIMIntro.csv',nCells=nCellsInExperiment,species2save={'mRNA'});
            Model_chg = Model_chg.loadData('dataForFIMIntro.csv', {'mRNA', 'exp1_mRNA'});

            pars = [Model_chg.parameters{:,2}];

            likelihoods = zeros(size(obj.kon_domain));
            for i = 1:length(obj.kon_domain)
                pars(1) = obj.kon_domain(i);
                l = Model_chg.computeLikelihood(pars(1));
                likelihoods(i) = l;
            end

            axes(ax2A)
            xlim(opts.xlims)
            ylim(opts.ylims)
            plot(obj.kon_domain, likelihoods, 'LineWidth', 1.5)
            hold on

            ax2A.Box = 'on';
            ax2A.LineWidth = 1.5;
            ax2A.FontSize = 11;
            ax2A.FontWeight = 'bold';
            ax2A.XColor = 'k';
            ax2A.YColor = 'k';
            ax2A.TickLength = [0.015 0.015];

            % Find MLE
            [~, max_idx] = max(likelihoods);
            mle = obj.kon_domain(max_idx);

            % True parameter
            xline(obj.star_kon, 'k--', 'LineWidth', 2)

            % MLE
            xline(mle, 'r-', 'LineWidth', 2)

            set(ax2A, 'XScale', 'log')

            xlabel('k_{on}')
            ylabel('Log-Likelihood')
            legend('Likelihood', 'True k_{on}', 'MLE', 'Location', 'best')
            % grid on
        end
        function prepareFigs2B(obj)
            arguments
                obj
            end
            %% MLE FIM relationship - Multiple cell - Bursting Model - Compute
            Model_chg = obj.ModelKoffSig;
            Model_chg.fittingOptions.modelVarsToFit = [1];
            log_probs_v_pars = zeros(ceil(Model_chg.fspOptions.bounds(4)), length(obj.kon_domain));
            for i = 1:length(obj.kon_domain)
                Model_chg.parameters{1,2} = obj.kon_domain(i);
                Model_chg = Model_chg.solve(solver='fsp');
                log_probs_v_pars(1:Model_chg.Solutions.fsp{1}.p.data.size(2), i) =...
                    log(double(Model_chg.Solutions.fsp{1}.p.sumOver(1).data));
            end
            save('probs_v_pars.mat', 'log_probs_v_pars')

        end
        function obj = makeFigs2BtoF(obj,ax2B,ax2C,ax2D,f2E,ax2E2,ax2F1,ax2F2,opts)
            arguments
                obj
                ax2B
                ax2C
                ax2D
                f2E
                ax2E2
                ax2F1
                ax2F2
                opts.Nsamps = 20
                opts.xlims = [10^1, 10^2.5]
                opts.ylims = [-40, -1]
            end

            %% MLE FIM relationship - Multiple cell - Bursting Model - Plotting
            load("probs_v_pars.mat",'log_probs_v_pars');

            T = readtable('dataForFIMIntro.csv');
            samples = T{:, 2:end};

            Nsamps = opts.Nsamps;
            sumLogL = NaN*ones(Nsamps,length(obj.kon_domain));
            for i = 1:Nsamps
                sumLogL(i,:) =  sum(log_probs_v_pars(samples(:,i)+1,:));
            end

            [~, max_idx] = max(sumLogL, [], 2);
            mle = obj.kon_domain(max_idx);

            axes(ax2B);
            hold(ax2B,'on')

            plot(ax2B,obj.kon_domain, sumLogL', 'lineWidth', 1.5)
            set(gca, 'XScale', 'log')

            ylim(opts.ylims)
            ylims = ylim;
            ymin = ylims(1);

            plot(mle, ymin * ones(size(mle)), 'rx', ...
                'MarkerSize', 12, 'LineWidth', 2)

            xline(obj.star_kon, 'k--', 'lineWidth', 2)

            ax2B.Box = 'on';
            ax2B.LineWidth = 1.5;
            ax2B.FontSize = 11;
            ax2B.FontWeight = 'bold';
            ax2B.XColor = 'k';
            ax2B.YColor = 'k';
            ax2B.TickLength = [0.015 0.015];
            % xlim([10^1, 10^2.5])

            %% Fig 2C

            Nsamps = size(samples,2);
            sumLogL = NaN*ones(Nsamps,length(obj.kon_domain));
            for i = 1:Nsamps
                sumLogL(i,:) =  sum(log_probs_v_pars(samples(:,i)+1,:));
            end

            [~, max_idx] = max(sumLogL, [], 2);
            mle = obj.kon_domain(max_idx);
            mle_log = log(mle);
            mleVar_log = var(mle_log);

            axes(ax2C);
            % Define bins uniformly in log10 space
            nBins = 25;
            xlims = [min(mle), max(mle)];

            log_edges = linspace(log10(xlims(1)), ...
                log10(xlims(2)), nBins+1);

            bin_edges = 10.^log_edges;

            histogram(mle, bin_edges, ...
                'Normalization', 'pdf', ...
                'FaceColor', [0.2 0.5 0.8], ...
                'EdgeColor', 'none')

            hold on

            xline(obj.star_kon, 'k--', 'LineWidth', 2)
            xline(mean(mle), 'r', 'LineWidth', 2)

            set(ax2C, 'XScale', 'log')
            xlim(xlims)

            xlabel('MLE k_{on}')
            ylabel('Probability Density')
            % grid on

            ax2C.Box = 'on';
            ax2C.LineWidth = 1.5;
            ax2C.FontSize = 11;
            ax2C.FontWeight = 'bold';
            ax2C.XColor = 'k';
            ax2C.YColor = 'k';
            ax2C.TickLength = [0.015 0.015];

            %% 2D MLE FIM relationship - sensitivity and FIM prediction
            dLogL_dkon = 0*sumLogL;

            for j = 1:size(sumLogL,1)
                dLogL_dkon(j,:) = gradient(sumLogL(j,:), obj.kon_domain);
            end

            % Evaluate sensitivity at the true parameter
            sensitivity = zeros(1,size(samples,2));

            for j = 1:size(sumLogL,1)
                sensitivity(j) = interp1(obj.kon_domain, dLogL_dkon(j,:), ...
                    obj.star_kon, 'linear');
            end

            % Sensitivity squared
            sensitivity_squared = sensitivity.^2;

            % Plot histogram
            axes(ax2D);
            hold on;
            histogram(sensitivity_squared, 50, 'Normalization', 'pdf')
            xlabel('Sensitivity^2')
            ylabel('Probability Density')
            title('Sensitivity Squared at True Parameter')

            % xlim([0,0.02])

            % Empirical MLE variance -> information in log space
            % xline(1/mleVar_log, 'b--', 'LineWidth', 2)

            % Model FIM
            Model_chg = obj.ModelKoffSig;
            Model_chg.fittingOptions.modelVarsToFit = [1];
            FIM = Model_chg.computeFIM();
            FIMEstimate = FIM{1};
            xline(FIMEstimate*obj.nCellsFig2, 'r', 'LineWidth', 2)

            xline(mean(sensitivity_squared), 'k--', 'LineWidth', 2)

            ax2D.Box = 'on';
            ax2D.LineWidth = 1.5;
            ax2D.FontSize = 11;
            ax2D.FontWeight = 'bold';
            ax2D.XColor = 'k';
            ax2D.YColor = 'k';
            ax2D.TickLength = [0.015 0.015];
            % legend('Sensitivity^2', ...
            %        'Mean sensitivity^2', ...
            %        '1 / MLE variance', ...
            %        'FIM')

            %% 2 MLE and sensitivity for multiple cell numbers
            % Can contain as many cell numbers as you want
            cell_numbers = [2 4 10 50 100 200 500];
            nSets = 10000;

            % All available single-cell likelihood curves
            L_all = log_probs_v_pars(samples(:)+1,:);

            MLEs = cell(length(cell_numbers),1);
            SensitivitySquared = cell(length(cell_numbers),1);

            % Calculate MLEs and sensitivity for all cell numbers

            for n = 1:length(cell_numbers)

                nCells = cell_numbers(n);

                MLEs{n} = zeros(nSets,1);
                SensitivitySquared{n} = zeros(nSets,1);

                for k = 1:nSets

                    % Randomly select cells
                    idx = randperm(size(L_all,1), nCells);

                    % Sum log-likelihoods across cells
                    L_sum = sum(L_all(idx,:), 1);

                    % Find MLE
                    [~, max_idx] = max(L_sum);
                    MLEs{n}(k) = obj.kon_domain(max_idx);

                    % Sensitivity of summed log-likelihood
                    dL_dkon = gradient(L_sum, obj.kon_domain);

                    % Evaluate sensitivity at true kon
                    sensitivity = interp1(obj.kon_domain, dL_dkon, ...
                        obj.star_kon, 'linear');

                    % Sensitivity squared
                    SensitivitySquared{n}(k) = sensitivity^2;
                end

                fprintf('Finished %d cells\n', nCells);
            end


            % Plot ONLY the first 4 cell numbers

            % figure(105);
            % clf
            figure(f2E)

            nPlot = min(4, length(cell_numbers));

            reference_MLE = MLEs{1};

            xlims = [min(reference_MLE), max(reference_MLE)];

            nBins = 32;

            log_edges = linspace(log10(xlims(1)), ...
                log10(xlims(2)), nBins+1);

            bin_edges = 10.^log_edges;

            for n = 1:nPlot

                subplot(2,2,n)

                histogram(MLEs{n}, bin_edges, ...
                    'Normalization', 'pdf', ...
                    'FaceColor', [0.2 0.5 0.8], ...
                    'FaceAlpha', 0.45, ...
                    'EdgeColor', 'none')

                hold on

                % True parameter
                xline(obj.star_kon, 'k--', 'LineWidth', 2)

                % Mean MLE
                xline(mean(MLEs{n}), 'r-', 'LineWidth', 2)

                set(gca, 'XScale', 'log')
                xlim(xlims)

                xlabel('MLE k_{on}')
                ylabel('Probability Density')
                title(sprintf('%d Cells', cell_numbers(n)))

                % grid on
            end

            ax = gca;
            ax.Box = 'on';
            ax.LineWidth = 1.5;
            ax.FontSize = 11;
            ax.FontWeight = 'bold';
            ax.XColor = 'k';
            ax.YColor = 'k';
            ax.TickLength = [0.015 0.015];


            %% 2E MLE variance and FIM Convergence

            mleVar = zeros(size(cell_numbers));
            fimVar = zeros(size(cell_numbers));

            for n = 1:length(cell_numbers)

                % Empirical variance of MLE
                mleVar(n) = var(MLEs{n});

                % FIM prediction
                fimVar(n) = 1 / (cell_numbers(n) * FIMEstimate);

            end

            axes(ax2E2);

            plot(cell_numbers, mleVar, 'ko-', ...
                'LineWidth', 2, ...
                'MarkerFaceColor', 'k')
            hold on

            plot(cell_numbers, fimVar, 'r^-', ...
                'LineWidth', 2, ...
                'MarkerFaceColor', 'r')

            xlabel('Number of Cells')
            ylabel('Variance of k_{on}')
            legend('MLE variance', 'FIM prediction', ...
                'Location', 'best')

            % grid on

            set(ax2E2, 'XScale', 'log')
            set(ax2E2, 'YScale', 'log')

            ax2E2.Box = 'on';
            ax2E2.LineWidth = 1;
            ax2E2.FontSize = 11;
            ax2E2.FontWeight = 'bold';
            ax2E2.XColor = 'k';
            ax2E2.YColor = 'k';
            ax2E2.TickLength = [0.008 0.008];

            %% 2F - MSE of MLE vs number of cells
            mleMSE = zeros(size(cell_numbers));

            for n = 1:length(cell_numbers)

                % MLE estimates for this number of cells
                estimates = MLEs{n};

                % Mean squared error relative to true parameter
                mleMSE(n) = mean((estimates - obj.star_kon).^2);

            end

            % Plot
            axes(ax2F1)

            plot(cell_numbers, mleMSE, 'ko-', ...
                'LineWidth', 2, ...
                'MarkerFaceColor', 'k')

            set(gca, 'XScale', 'log')
            set(gca, 'YScale', 'log')

            xlabel('Number of Cells')
            ylabel('MLE MSE')
            title('MLE Mean Squared Error vs Number of Cells')

            % grid on

            ax2F1.Box = 'on';
            ax2F1.LineWidth = 1.5;
            ax2F1.FontSize = 11;
            ax2F1.FontWeight = 'bold';
            ax2F1.XColor = 'k';
            ax2F1.YColor = 'k';
            ax2F1.TickLength = [0.015 0.015];

            %% MSE between MLE variance and FIM variance
            varMSE = (mleVar - fimVar).^2;

            % figure(108);
            % clf
            axes(ax2F2)

            plot(cell_numbers, varMSE, 'ko-', ...
                'LineWidth', 2, ...
                'MarkerFaceColor', 'k')

            set(gca, 'XScale', 'log')
            set(gca, 'YScale', 'log')

            xlabel('Number of Cells')
            ylabel('MSE: MLE Variance vs FIM Variance')
            title('MSE Between Empirical Variance and FIM Estimate')

            % grid on

            ax2F2.Box = 'on';
            ax2F2.LineWidth = 1;
            ax2F2.FontSize = 11;
            ax2F2.FontWeight = 'bold';
            ax2F2.XColor = 'k';
            ax2F2.YColor = 'k';
            ax2F2.TickLength = [0.008 0.008];

        end
        function prepareFig2G(obj,opts)
            arguments
                obj
                opts.nMLE = 200;
            end
            Model_chg = obj.ModelKoffSig;

            %% Verification of FIM using CRLB (spread of MLE)
            Model_chg.fittingOptions.modelVarsToFit = [1:2];
            MLE = Model_chg.estimateMLEspread(nCells=obj.nCellsInExperiment,observableSpecies={'mRNA'},nMLE=opts.nMLE,simsSaveFile='BurstFIMSims.csv',freePars=[1:2],restart=true);
            MLE = Model_chg.estimateMLEspread(nCells=obj.nCellsInExperiment,observableSpecies={'mRNA'},nMLE=opts.nMLE,simsSaveFile='BurstFIMSims.csv',freePars=[1:2],startPars=exp(MLE.mhSamples),restart=false);
            save('MLEForCRLBVerifications.mat', 'MLE')
        end
        function makeFig2G(obj)
            Model_chg = obj.ModelKoffSig;
            Model_chg.fittingOptions.modelVarsToFit = [1:2];
            load('MLEForCRLBVerifications.mat', 'MLE')

            FIMs = Model_chg.computeFIM(scale='log',freePars=[1:2],...
                observed={'mRNA'});
            FIMTotal = Model_chg.totalFim(FIMs,obj.nCellsInExperiment);
            Model_chg.plotMHResults(MLE,FIM=FIMTotal,fimScale='log',truncateChain=false);
        end
        function makeFig2H(obj, f2g1, f2g2, f2h1, f2h2)
            f1 = figure(f2g1); % fim ellipse
            clf
            f2 = figure(150); % default fim analysis
            clf

            Model_chg = obj.ModelKoffSig;
            Model_chg.fittingOptions.modelVarsToFit = [1:2];
            load('MLEForCRLBVerifications.mat', 'MLE')

            FIMs = Model_chg.computeFIM(scale='log',freePars=[1:2],...
                observed={'mRNA'});
            FIMTotal = Model_chg.totalFim(FIMs,obj.nCellsInExperiment);

            FIM = FIMTotal{1};

            Model_chg.plotFIMResults(FIM^(-1)/log(10)^2, 'log',...
                Model_chg.parameters(1:2,1),...
                [Model_chg.parameters{1:2,2}],...
                PlotEllipses=true, ...
                EllipseFigure=f1,...
                Colors = struct('EllipseColors',[0, 0, 0],'CenterSquare',[0,0,0]), ...
                EllipsePairs=[1,2], ...
                FigureHandle=f2,...
                LogThreshold=-4,...
                HeatMapType='invfim',...
                MatrixType='invfim');

            hold on

            C = FIM^(-1)/log(10)^2;

            % Parameters corresponding to your ellipse pair
            C2 = C([1 2],[1 2]);

            % Eigenvectors/eigenvalues
            [V,D] = eig(C2);

            % Sort eigenvalues from smallest to largest
            [lambda,idx] = sort(diag(D));
            V = V(:,idx);

            % Center of ellipse
            x0 = log10(Model_chg.parameters{2,2});
            y0 = log10(Model_chg.parameters{1,2});

            % Scale factor for visualization
            scale = 2;

            % MLE estimates
            MLElog = [MLE.mhSamples(:,2)/log(10), ...
                MLE.mhSamples(:,1)/log(10)];

            % MLE mean and covariance
            muMLE = mean(MLElog,1)';
            CMLE = cov(MLElog);

            % Eigenvectors/eigenvalues of MLE covariance
            [VMLE,DMLE] = eig(CMLE);

            [lambdaMLE,idxMLE] = sort(diag(DMLE),'descend');
            VMLE = VMLE(:,idxMLE);

            % 95% confidence ellipse
            chi2val = icdf('chi2',0.9,2);

            aMLE = sqrt(chi2val*lambdaMLE(1));
            bMLE = sqrt(chi2val*lambdaMLE(2));

            t = linspace(0,2*pi,300);

            MLEellipse = VMLE * ...
                [aMLE*cos(t); bMLE*sin(t)];

            % MLE estimates
            scatter(MLElog(:,1), ...
                MLElog(:,2), ...
                10, [0.5 0.5 0.5], 'filled');

            % Small eigenvalue direction
            quiver(x0,y0,...
                V(2,1)*sqrt(lambda(1))*scale,...
                V(1,1)*sqrt(lambda(1))*scale,...
                0,...
                'LineWidth',2,...
                'Color','r',...
                'MaxHeadSize',0.5);

            % Large eigenvalue direction
            quiver(x0,y0,...
                V(2,2)*sqrt(lambda(2))*scale,...
                V(1,2)*sqrt(lambda(2))*scale,...
                0,...
                'LineWidth',2,...
                'Color','b',...
                'MaxHeadSize',0.5);

            ax = gca;
            ax.Box = 'on';
            ax.LineWidth = 1.5;
            ax.FontSize = 11;
            ax.FontWeight = 'bold';
            ax.XColor = 'k';
            ax.YColor = 'k';
            ax.TickLength = [0.015 0.015];

            % MLE covariance ellipse
            plot(muMLE(1) + MLEellipse(1,:), ...
                muMLE(2) + MLEellipse(2,:), ...
                'c-', ...
                'LineWidth',2);

            % MLE mean
            % plot(muMLE(1),muMLE(2),...
            %     'ro',...
            %     'MarkerSize',8,...
            %     'MarkerFaceColor','r',...
            %     'LineWidth',2);


            figure(f2h1); % heatmap of I^(-1)
            clf
            obj.plotHeatmap(C, {'k_{on}', 'k_{off}'}, {'k_{on}', 'k_{off}'}, 'I^{-1}')


            figure(f2h2); % heatmap of eig(I^(-1))
            clf

            [V, D] = eig(C);

            % Sort eigenvalues from largest to smallest
            [lambda, idx] = sort(diag(D), 'descend');

            % Reorder eigenvectors to match
            V = V(:, idx);

            % Rebuild diagonal eigenvalue matrix
            D = sqrt(diag(lambda));

            % Plot V*D
            obj.plotHeatmap(D, ...
                {'k_{on}', 'k_{off}'}, ...
                {'\lambda_{1}', '\lambda_{2}'}, ...
                'V(I^{-1}) \lambda(I^{-1})')

            % Make eigenvector orientation deterministic
            if V(1,1) < 0
                V(:,1) = -V(:,1);
            end

            if det(V) < 0
                V(:,2) = -V(:,2);
            end

            % Center in the SAME parameter ordering as C and V:
            % [kon, koff]
            mu = [log10(Model_chg.parameters{1,2}); ...
                log10(Model_chg.parameters{2,2})];

            % Chi-square scaling
            chi2val = icdf('chi2',0.9,2);

            % Principal-axis lengths
            a = sqrt(chi2val * lambda(1));
            b = sqrt(chi2val * lambda(2));

            % Parameterize ellipse
            t = linspace(0,2*pi,300);

            xEllipse = a*cos(t);
            yEllipse = b*sin(t);

            % Put MLE estimates into [kon, koff] ordering
            MLE_FIMorder = [MLE.mhSamples(:,1)/log(10), ...
                MLE.mhSamples(:,2)/log(10)];

            % Rotate MLE estimates into FIM eigenvector coordinates
            MLErot = V' * (MLE_FIMorder' - mu);

            % MLE mean in FIM eigenvector coordinates
            muMLE_FIMorder = mean(MLE_FIMorder,1)';

            muMLERot = V' * (muMLE_FIMorder - mu);

            % MLE covariance in FIM eigenvector coordinates
            CMLE = cov(MLE_FIMorder);

            CMLErot = V' * CMLE * V;

            % Eigenvectors/eigenvalues of rotated MLE covariance
            [VMLErot,DMLErot] = eig(CMLErot);

            [lambdaMLErot,idxMLErot] = sort(diag(DMLErot),'descend');
            VMLErot = VMLErot(:,idxMLErot);

            % MLE ellipse in rotated coordinates
            aMLERot = sqrt(chi2val*lambdaMLErot(1));
            bMLERot = sqrt(chi2val*lambdaMLErot(2));

            MLEellipseRot = VMLErot * ...
                [aMLERot*cos(t); bMLERot*sin(t)];

            % Plot in eigenvector coordinates
            figure(f2g2);
            clf;
            hold on;

            % FIM ellipse
            plot(xEllipse,...
                yEllipse,...
                'k-',...
                'LineWidth',2);

            % MLE estimates
            scatter(MLErot(1,:),...
                MLErot(2,:),...
                10,...
                [0.5 0.5 0.5],...
                'filled');

            % MLE covariance ellipse
            plot(muMLERot(1) + MLEellipseRot(1,:),...
                muMLERot(2) + MLEellipseRot(2,:),...
                'c-',...
                'LineWidth',2);

            % FIM center
            plot(0,0,...
                'ks',...
                'MarkerSize',8,...
                'MarkerFaceColor','w',...
                'LineWidth',2);

            % MLE mean
            % plot(muMLERot(1),muMLERot(2),...
            %     'ro',...
            %     'MarkerSize',8,...
            %     'MarkerFaceColor','r',...
            %     'LineWidth',2);

            % Principal axes
            quiver(0,0,...
                a,0,...
                0,...
                'b',...
                'LineWidth',2,...
                'MaxHeadSize',0.5);

            quiver(0,0,...
                0,b,...
                0,...
                'r',...
                'LineWidth',2,...
                'MaxHeadSize',0.5);

            xlabel('Largest variance eigenvector');
            ylabel('Smallest variance eigenvector');

            grid on

            ax = gca;
            ax.Box = 'on';
            ax.LineWidth = 1.5;
            ax.FontSize = 11;
            ax.FontWeight = 'bold';
            ax.XColor = 'k';
            ax.YColor = 'k';
            ax.TickLength = [0.015 0.015];

            axis equal
        end
        function makeFigs3A(obj,f3A1,f3A2,f3A3,f3A4,f3A5)

            Model_chg = obj.ModelKoffSig;

            Model_chg.fittingOptions.modelVarsToFit = [1:5];
            FIMs = Model_chg.computeFIM(scale='log',freePars=[1:5],...
                observed={'mRNA'});
            FIMTotal = Model_chg.totalFim(FIMs,obj.nCellsInExperiment);

            FIM = FIMTotal{1};

            Model_chg.plotFIMResults(FIM^(-1)/log(10)^2, 'log',...
                Model_chg.parameters(1:5,1),...
                [Model_chg.parameters{1:5,2}],...
                PlotEllipses=true, ...
                EllipseFigure=f3A1,...
                Colors = struct('EllipseColors',[0, 0, 0],'CenterSquare',[0,0,0]), ...
                EllipsePairs=[1,2], ...
                FigureHandle=f3A2,...
                LogThreshold=-4,...
                HeatMapType='invfim',...
                MatrixType='invfim');

            hold on

            C = FIM^(-1)/log(10)^2;

            % Parameters corresponding to your ellipse pair
            C2 = C([1 2],[1 2]);

            % Eigenvectors/eigenvalues
            [V,D] = eig(C2);

            % Sort eigenvalues from smallest to largest
            [lambda,idx] = sort(diag(D));
            V = V(:,idx);

            % Center of ellipse
            x0 = log10(Model_chg.parameters{2,2});
            y0 = log10(Model_chg.parameters{1,2});

            % Scale factor for visualization
            scale = 2;

            % Small eigenvalue direction
            quiver(x0,y0,...
                V(2,1)*sqrt(lambda(1))*scale,...
                V(1,1)*sqrt(lambda(1))*scale,...
                0,...
                'LineWidth',2,...
                'Color','r',...
                'MaxHeadSize',0.5);

            % Large eigenvalue direction
            quiver(x0,y0,...
                V(2,2)*sqrt(lambda(2))*scale,...
                V(1,2)*sqrt(lambda(2))*scale,...
                0,...
                'LineWidth',2,...
                'Color','b',...
                'MaxHeadSize',0.5);
            ax = gca;
            ax.Box = 'on';
            ax.LineWidth = 1.5;
            ax.FontSize = 11;
            ax.FontWeight = 'bold';
            ax.XColor = 'k';
            ax.YColor = 'k';
            ax.TickLength = [0.015 0.015];

            figure(f3A3); % heatmap of I^(-1)
            clf
            obj.plotHeatmap(C, {'k_{on,init}', 'k_{off}', 'k_r', '\gamma', 'k_{on,final}'}, {'k_{on,init}', 'k_{off}', 'k_r', '\gamma', 'k_{on,final}'}, 'I^{-1}')


            figure(f3A4); % heatmap of eig(I^(-1))
            clf

            [V, D] = eig(C);

            % Sort eigenvalues from largest to smallest
            [lambda, idx] = sort(diag(D), 'descend');

            % Reorder eigenvectors to match
            V = V(:, idx);

            % Rebuild diagonal eigenvalue matrix
            D = sqrt(diag(lambda));

            % Plot V*D
            obj.plotHeatmap(D, ...
                {'k_{on,init}', 'k_{off}', 'k_r', '\gamma', 'k_{on,final}'}, ...
                {'\lambda_{1}', '\lambda_{2}', '\lambda_{3}', '\lambda_{4}', '\lambda_{5}'}, ...
                'V(I^{-1}) \lambda(I^{-1})')

            % Make eigenvector orientation deterministic
            % Largest eigenvector should point generally in +x direction
            if V(1,1) < 0
                V(:,1) = -V(:,1);
            end

            % Make second eigenvector form a right-handed coordinate system
            if det(V) < 0
                V(:,2) = -V(:,2);
            end

            % Center
            x0 = log10(Model_chg.parameters{2,2});
            y0 = log10(Model_chg.parameters{1,2});
            mu = [x0; y0];

            % Chi-square scaling
            chi2val = icdf('chi2',0.95,2);

            % Principal-axis lengths
            a = sqrt(chi2val * lambda(1));   % LARGE variance -> x
            b = sqrt(chi2val * lambda(2));   % SMALL variance  -> y

            % Parameterize ellipse DIRECTLY in eigenvector coordinates
            t = linspace(0,2*pi,300);

            xEllipse = a*cos(t);
            yEllipse = b*sin(t);

            % Plot in eigenvector coordinates
            figure(f3A5);
            clf;
            hold on;

            plot(xEllipse,yEllipse,...
                'k-',...
                'LineWidth',2);

            plot(0,0,...
                'ks',...
                'MarkerSize',8,...
                'MarkerFaceColor','w',...
                'LineWidth',2);

            % Principal axes
            quiver(0,0,...
                a,0,...
                0,...
                'b',...
                'LineWidth',2,...
                'MaxHeadSize',0.5);

            quiver(0,0,...
                0,b,...
                0,...
                'r',...
                'LineWidth',2,...
                'MaxHeadSize',0.5);

            xlabel('Largest variance eigenvector');
            ylabel('Smallest variance eigenvector');

            axis equal;
            grid on;

            ax = gca;
            ax.Box = 'on';
            ax.LineWidth = 1.5;
            ax.FontSize = 11;
            ax.FontWeight = 'bold';
            ax.XColor = 'k';
            ax.YColor = 'k';
            ax.TickLength = [0.015 0.015];

        end

        function obj = prepareFIMS(obj)
            Model_chg = obj.ModelKoffSig;
            Model_chg.solutionScheme = 'fspsens';
            Model_chg.fspOptions.fspTol = 1e-5;
            obj.FIM = cell(length(obj.Sarray),length(obj.Sarray),length(Model_chg.tSpan));
            for iS1 = 1:length(obj.Sarray)
                for iS2 = 1:length(obj.Sarray)
                    Model_chg = Model_chg.changeParameter({'S1',obj.Sarray(iS1);'S2',obj.Sarray(iS2)});
                    obj.FIM(iS1,iS2,:) = Model_chg.computeFIM(freePars=(1:5),scale='log');
                end
            end

        end

        function obj = makeFig3B(obj,f3B1,f3B2,f3B3,f3B4)
            Model_chg = obj.ModelKoffSig;
            obj.tt = linspace(min(Model_chg.tSpan), max(Model_chg.tSpan), 100);

            % exp 1 - steady states
            obj.exp1NCells = zeros(size(obj.FIM));
            obj.exp1NCells(5,5, 1) = obj.Ncells/2;
            obj.exp1NCells(1,1, 1) = obj.Ncells/2;

            figure(f3B1)
            hold on

            f = @(t)(t > 1)*obj.Sarray(end) + (t <= 1)*obj.Sarray(1);

            xx = [Model_chg.tSpan(1), Model_chg.tSpan(31)];
            yy = [f(Model_chg.tSpan(1)), f(Model_chg.tSpan(31))];

            plot(obj.tt, f(obj.tt), 'k-', 'LineWidth', 3)

            dx = obj.tt(end) - obj.tt(end-1);
            dy = f(obj.tt(end)) - f(obj.tt(end-1));

            quiver(obj.tt(end-1), f(obj.tt(end-1)), dx, dy, 0, ...
                'Color', 'k', ...
                'LineWidth', 3, ...
                'MaxHeadSize', 100);

            % Alpha proportional to number of cells
            alpha1 = (obj.Ncells) / obj.Ncells;
            alpha = min(alpha1 + 0.2, 1);

            scatter(xx, yy, 1000, 'r', 'x', ...
                'LineWidth', 3, ...
                'MarkerEdgeAlpha', alpha);

            ylim([-5,35])
            ax = gca;
            ax.LineWidth = 1.5;

            % exp 2 - steady state plus one during transition
            obj.exp2NCells = zeros(size(obj.FIM));
            obj.exp2NCells(1,5, [1,6,31]) = obj.Ncells/3;

            figure(f3B2)
            hold on

            f = @(t)(t > 1)*obj.Sarray(end) + (t <= 1)*obj.Sarray(1);

            xx = [Model_chg.tSpan(1), Model_chg.tSpan(6), Model_chg.tSpan(31)];
            yy = [f(Model_chg.tSpan(1)), f(Model_chg.tSpan(6)), f(Model_chg.tSpan(31))];

            plot(obj.tt, f(obj.tt), 'k-', 'LineWidth', 3)
            dx = obj.tt(end) - obj.tt(end-1);
            dy = f(obj.tt(end)) - f(obj.tt(end-1));
            quiver(obj.tt(end-1), f(obj.tt(end-1)), dx, dy, 0, ...
                'Color', 'k', ...
                'LineWidth', 3, ...
                'MaxHeadSize', 100);

            % Alpha proportional to number of cells
            alpha2 = (obj.Ncells*(2/3)) / obj.Ncells;
            alpha = min(alpha2 + 0.2, 1);

            scatter(xx, yy, 1000, 'r', 'x', ...
                'LineWidth', 3, ...
                'MarkerEdgeAlpha', alpha);

            ylim([-5,35])
            ax = gca;
            ax.LineWidth = 1.5;

            % exp 3 - multiple of exp 2
            obj.exp3NCells = zeros(size(obj.FIM));
            obj.exp3NCells(1, 5, [1,6,31]) = obj.Ncells/6;
            obj.exp3NCells(5, 3, [1,6,31]) = obj.Ncells/6;

            figure(f3B3)
            hold on

            f = @(t)(t > 1)*obj.Sarray(end) + (t <= 1)*obj.Sarray(1);

            plot(obj.tt, f(obj.tt), 'k-', 'LineWidth', 3)
            dx = obj.tt(end) - obj.tt(end-1);
            dy = f(obj.tt(end)) - f(obj.tt(end-1));
            quiver(obj.tt(end-1), f(obj.tt(end-1)), dx, dy, 0, ...
                'Color', 'k', ...
                'LineWidth', 3, ...
                'MaxHeadSize', 100);

            xx = [Model_chg.tSpan(1), Model_chg.tSpan(6), Model_chg.tSpan(31)];
            yy = [f(Model_chg.tSpan(1)), f(Model_chg.tSpan(6)), f(Model_chg.tSpan(31))];

            % Alpha for Ncells/6
            alpha3 = (obj.Ncells/3) / obj.Ncells;
            alpha = min(alpha3 + 0.2, 1);

            scatter(xx, yy, 1000, 'r', 'x', ...
                'LineWidth', 3, ...
                'MarkerEdgeAlpha', alpha3);

            f = @(t)(t > 1)*obj.Sarray(3) + (t <= 1)*obj.Sarray(end);

            plot(obj.tt, f(obj.tt), 'b-', 'LineWidth', 3)
            dx = obj.tt(end) - obj.tt(end-1);
            dy = f(obj.tt(end)) - f(obj.tt(end-1));
            quiver(obj.tt(end-1), f(obj.tt(end-1)), dx, dy, 0, ...
                'Color', 'b', ...
                'LineWidth', 3, ...
                'MaxHeadSize', 100);


            xx = [Model_chg.tSpan(1), Model_chg.tSpan(6), Model_chg.tSpan(31)];
            yy = [f(Model_chg.tSpan(1)), f(Model_chg.tSpan(6)), f(Model_chg.tSpan(31))];

            scatter(xx, yy, 1000, 'r', 'x', ...
                'LineWidth', 3, ...
                'MarkerEdgeAlpha', alpha);

            ylim([-5,35])
            ax = gca;
            ax.LineWidth = 1.5;

            figure(f3B4)
            clf
            hold on

            J = find(obj.OptExperiment);
            experiments = [];
            ssExperiments = [];

            for j = 1:length(J)
                optimizedParams = obj.OptExperiment(J(j));
                paramIndices = obj.indsFims(J(j), :);

                if paramIndices(1) == paramIndices(2)
                    % steady state
                    ssExperiments = [ssExperiments; paramIndices, optimizedParams];
                else
                    % new experiment
                    experiments = [experiments; paramIndices, optimizedParams];
                end
            end
            experiments

            [uniqueExperiments, ~, ic] = unique(experiments(:, 1:2), 'rows');

            nUniqueExperiments = size(uniqueExperiments, 1);

            colors = {'k', 'b', 'g'};

            for i = 1:size(uniqueExperiments, 1)

                % First two columns define the experiment
                startState = uniqueExperiments(i, 1);
                endState   = uniqueExperiments(i, 2);

                % Step function
                f = @(t) (t > 1)*obj.Sarray(endState) + ...
                    (t <= 1)*obj.Sarray(startState);

                % Draw step function
                plot(obj.tt, f(obj.tt), '-', ...
                    'Color', colors{mod(i, length(colors))+1}, ...
                    'LineWidth', 3);

                % Add arrow at the end of the curve
                dx = obj.tt(end) - obj.tt(end-1);
                dy = f(obj.tt(end)) - f(obj.tt(end-1));

                quiver(obj.tt(end-1), f(obj.tt(end-1)), dx, dy, 0, ...
                    'Color', colors{mod(i, length(colors))+1}, ...
                    'LineWidth', 3, ...
                    'MaxHeadSize', 100);

            end

            ylim([-5, 35])

            for i = 1:size(ssExperiments, 1)
                ssConc = ssExperiments(i,1);

                possibleExperiments = find(any(uniqueExperiments == ssConc, 2));
                matches = uniqueExperiments(possibleExperiments,:) == ssConc;
                col = zeros(size(possibleExperiments));

                col(matches(:,1)) = 1;
                col(matches(:,2)) = 2;

                if ~isempty(possibleExperiments)
                    for j = 1:length(possibleExperiments)

                        tIndex = 1 + (col(j) == 2)*30;

                        newExperiment = [uniqueExperiments(possibleExperiments(j),:), ...
                            tIndex, ...
                            ssExperiments(i,4)/length(possibleExperiments)];

                        % Look for an existing entry matching the first 3 columns
                        match = all(experiments(:,1:3) == newExperiment(1:3), 2);

                        if any(match)
                            % Add cells to existing entry
                            experiments(match,4) = experiments(match,4) + newExperiment(4);
                        else
                            % Append as a new entry
                            experiments = [experiments; newExperiment];
                        end

                    end
                end
            end
            experiments

            for i = 1:size(experiments, 1)

                initialConc = experiments(i,1);
                finalConc   = experiments(i,2);
                tIndex      = experiments(i,3);
                nCells      = experiments(i,4);

                f = @(t)(t > 1)*obj.Sarray(finalConc) + (t <= 1)*obj.Sarray(initialConc);

                x = Model_chg.tSpan(tIndex);

                alpha = min(nCells / obj.Ncells + 0.2, 1);

                scatter(x, f(x), 1000, 'r', 'x', ...
                    'LineWidth', 3, ...
                    'MarkerEdgeAlpha', alpha);

            end
            ax = gca;
            ax.LineWidth = 1.5;
            ylim([-5,35])

        end

        function obj = computeFIMDeterminants(obj)
            obj.FIM_Exp1 = zeros(size(obj.FIM{1}));
            for i = 1:size(obj.FIM,1)
                for j = 1:size(obj.FIM,2)
                    for k = 1:size(obj.FIM,3)
                        obj.FIM_Exp1 = obj.FIM_Exp1 + ...
                            obj.exp1NCells(i,j,k) .* obj.FIM{i,j,k};
                    end
                end
            end
            disp(['Determinant of FIM for experiment 1: ',num2str(det(obj.FIM_Exp1))])

            obj.FIM_Exp2 = zeros(size(obj.FIM{1}));
            for i = 1:size(obj.FIM,1)
                for j = 1:size(obj.FIM,2)
                    for k = 1:size(obj.FIM,3)
                        obj.FIM_Exp2 = obj.FIM_Exp2 + ...
                            obj.exp2NCells(i,j,k) .* obj.FIM{i,j,k};
                    end
                end
            end
            disp(['Determinant of FIM for experiment 2: ',num2str(det(obj.FIM_Exp2))])

            % Measurement at change from one SS to another at three time points.
            obj.FIM_Exp3 = zeros(size(obj.FIM{1}));
            for i = 1:size(obj.FIM,1)
                for j = 1:size(obj.FIM,2)
                    for k = 1:size(obj.FIM,3)
                        obj.FIM_Exp3 = obj.FIM_Exp3 + ...
                            obj.exp3NCells(i,j,k) .* obj.FIM{i,j,k};
                    end
                end
            end
            disp(['Determinant of FIM for experiment 3: ',num2str(det(obj.FIM_Exp3))])

        end

        function obj = optimizeExperiment(obj)
            Model_chg = obj.ModelKoffSig;

            obj.allFims = {}; %cell(numel(FIM),1);
            obj.indsFims = []; zeros(numel(obj.FIM),3);
            k = 0;
            for iS0 = 1:length(obj.Sarray)
                for iS1 = 1:length(obj.Sarray)
                    for iT = 1:length(Model_chg.tSpan)
                        k = k+1;
                        obj.allFims(k,1) = obj.FIM(iS0,iS1,iT);
                        obj.indsFims(k,:) = [iS0,iS1,iT];
                    end
                end
            end
            obj.OptExperiment = Model_chg.optimizeCellCounts(obj.allFims,obj.Ncells,'D-opt');
            J = find(obj.OptExperiment);
            disp(['Optimized Experiment Design:'])
            for j = 1:length(J)
                % Store optimized parameters and their corresponding indices
                optimizedParams = obj.OptExperiment(J(j));
                paramIndices = obj.indsFims(J(j), :);
                if paramIndices(3)==1||paramIndices(1)==paramIndices(2) % SS experiment
                    disp(['   ',num2str(optimizedParams),' cells at steady state for S0 = ',num2str(obj.Sarray(paramIndices(1)))])
                else
                    disp(['   ',num2str(optimizedParams),' cells at time ',num2str(Model_chg.tSpan(paramIndices(3))),' for S0 = ',num2str(obj.Sarray(paramIndices(1))),' and S1 = ',num2str(obj.Sarray(paramIndices(2)))])
                end
            end

            obj.FIM_Opt = 0;
            for i = 1:length(obj.OptExperiment)
                obj.FIM_Opt = obj.FIM_Opt + obj.OptExperiment(i)*obj.allFims{i};
            end
            disp(['Determinant of FIM for optimized measurements: ',num2str(det(obj.FIM_Opt))])

        end

        function obj = makeFig3C(obj,f3C1,f3C2,f3C3,f3C4)
            figure(f3C1);
            clf;
            fprintf('steady state experiment determinate %e\n', det(obj.FIM_Exp1^(-1)))
            obj.plotHeatmap(obj.FIM_Exp1^(-1), {'k_{on,init}', 'k_{off}', 'k_r', '\gamma', 'k_{on,final}'}, ...
                {'k_{on,init}', 'k_{off}', 'k_r', '\gamma', 'k_{on,final}'}, ...
                'I^{-1} - Exp 1')

            figure(f3C2);
            fprintf('transient experiment determinate %e\n', det(obj.FIM_Exp1^(-1)))
            obj.plotHeatmap(obj.FIM_Exp2^(-1), {'k_{on}', 'k_{off}', 'k_r', '\gamma', 'k_{on,final}'}, ...
                {'k_{on,init}', 'k_{off}', 'k_r', '\gamma', 'k_{on,final}'}, ...
                'I^{-1} - Exp 2')

            figure(f3C3);
            clf;
            fprintf('intuitive experiment determinate %e\n', det(obj.FIM_Exp1^(-1)))
            obj.plotHeatmap(obj.FIM_Exp3^(-1), {'k_{on}', 'k_{off}', 'k_r', '\gamma', 'k_{on,final}'}, ...
                {'k_{on,init}', 'k_{off}', 'k_r', '\gamma', 'k_{on,final}'}, ...
                'I^{-1} - Exp 3')

            figure(f3C4);
            clf;
            fprintf('optimal experiment determinate %e\n', det(obj.FIM_Exp1^(-1)))
            obj.plotHeatmap(obj.FIM_Opt^(-1), {'k_{on,init}', 'k_{off}', 'k_r', '\gamma', 'k_{on,final}'}, ...
                {'k_{on,init}', 'k_{off}', 'k_r', '\gamma', 'k_{on,final}'}, ...
                'I^{-1} - Exp Opt')

            % Use figure 209 as the reference
            refFig = figure(f3C2);
            refAx = gca;
            refCb = colorbar(refAx);

            % Get reference settings
            refCLim = refAx.CLim;
            refCMap = colormap(refAx);

            % Apply to the other figures
            for figNum = [f3C1 f3C3 f3C4]

                fig = figure(figNum);
                ax = gca;

                % Same colors
                colormap(ax, refCMap);

                % Same color scaling
                clim(ax, refCLim);

                % Get this figure's colorbar
                cb = colorbar(ax);

                % Same ticks and labels
                cb.Ticks = refCb.Ticks;
                cb.TickLabels = refCb.TickLabels;
                cb.TickLabelInterpreter = refCb.TickLabelInterpreter;
                cb.Label.String = refCb.Label.String;

            end


        end

        function obj = makeFig3D(obj,f3D1,f3D2,f3D3)

            Model_chg =obj.ModelKoffSig;
            % vNCells = round(logspace(2, 3, 10));
            vNCells = [100, 300, 600, 1000];

            nExperiments = 4;
            experimentFIMs = cell(length(vNCells), nExperiments);
            EOpt = zeros(length(vNCells), nExperiments);
            DOpt = zeros(length(vNCells), nExperiments);

            DsOptIdx = 4;
            DsOpt = zeros(length(vNCells), nExperiments);

            for a = 1:length(vNCells)
                Ncells = vNCells(a)

                % exp 1 - steady states
                exp1NCells = zeros(size(obj.FIM));
                exp1NCells(5,5, 1) = Ncells/2;
                exp1NCells(1,1, 1) = Ncells/2;

                % exp 2 - steady state plus one during transition
                exp2NCells = zeros(size(obj.FIM));
                exp2NCells(1,5, [1,6,31]) = Ncells/3;

                % exp 3 - multiple of exp 2
                exp3NCells = zeros(size(obj.FIM));
                exp3NCells(1, 5, [1,6,31]) = Ncells/6;
                exp3NCells(5, 3, [1,6,31]) = Ncells/6;

                experiments = {exp1NCells, exp2NCells, exp3NCells};

                for b = 1:length(experiments)
                    % compute fim for each experiment
                    cells4experiment = experiments{b};
                    runningFIM = zeros(size(obj.FIM{1}));
                    for i = 1:size(obj.FIM,1)
                        for j = 1:size(obj.FIM,2)
                            for k = 1:size(obj.FIM,3)
                                runningFIM = runningFIM + ...
                                    cells4experiment(i,j,k) .* obj.FIM{i,j,k};
                            end
                        end
                    end
                    experimentFIMs{a,b} = runningFIM;

                    % compute optimality
                    [V,D] = eig(runningFIM);
                    EOpt(a,b) = 1/min(diag(D));
                    DOpt(a,b) = 1/det(runningFIM);
                    % DsOpt(a,b) = runningFIM(DsOptIdx, DsOptIdx);
                    runningFIMInv = inv(runningFIM);
                    DsOpt(a,b) = runningFIMInv(DsOptIdx, DsOptIdx);

                end
                % compute optimal optimality
                % D opt
                OptExperiment = Model_chg.optimizeCellCounts(obj.allFims,Ncells,'D-cov');
                FIM_Opt = 0;
                for i = 1:length(OptExperiment)
                    FIM_Opt = FIM_Opt + OptExperiment(i)*obj.allFims{i};
                end
                DOpt(a,nExperiments) = 1/det(FIM_Opt);

                % E opt
                OptExperiment = Model_chg.optimizeCellCounts(obj.allFims,Ncells,'E-cov');
                FIM_Opt = 0;
                for i = 1:length(OptExperiment)
                    FIM_Opt = FIM_Opt + OptExperiment(i)*obj.allFims{i};
                end
                [V,D] = eig(FIM_Opt);
                EOpt(a,nExperiments) = 1/min(diag(D));

                % Ds opt
                OptExperiment = Model_chg.optimizeCellCounts(obj.allFims,Ncells,sprintf('D-cov-sub%d',DsOptIdx));
                FIM_Opt = 0;
                for i = 1:length(OptExperiment)
                    FIM_Opt = FIM_Opt + OptExperiment(i)*obj.allFims{i};
                end
                FIM_OptInv = inv(FIM_Opt);
                DsOpt(a,nExperiments) = FIM_OptInv(DsOptIdx,DsOptIdx);

            end

            %% Plot Optimality Criteria
            % Cell counts
            x = vNCells(:);

            % Matrices of performance metrics
            data = {DOpt, EOpt, DsOpt};
            names = {'DOpt', 'EOpt', 'DsOpt'};

            figs = [f3D1,f3D2,f3D3];
            for k = 1:length(data)

                figure(figs(k));
                clf

                Y = data{k};

                % Plot all experiments
                plot(x, Y, '-', 'LineWidth', 2);
                hold on;

                % ---- Fit lines in log-log space ----
                % Each column is an experiment
                coeff = zeros(2, size(Y,2));

                for j = 1:size(Y,2)
                    coeff(:,j) = polyfit(log10(x), log10(Y(:,j)), 1);
                end

                % Smooth x values for fitted lines
                xfit = logspace(log10(min(x)), log10(max(x)), 200);

                % Plot fitted lines
                for j = 1:size(Y,2)
                    yfit = 10.^(polyval(coeff(:,j), log10(xfit)));

                    % plot(xfit, yfit, '--', 'LineWidth', 1.5);
                end

                % ---- Optimal design at 300 cells ----
                xOpt = 150;

                % Evaluate optimal-design fit at 300 cells
                yOpt300 = 10.^(polyval(coeff(:,end), log10(xOpt)));

                % Horizontal line through optimality at 300 cells
                yline(yOpt300, 'k-', 'LineWidth', 2, ...
                    'DisplayName', 'Optimal @ 300 cells');

                % Mark optimal point
                % plot(xOpt, yOpt300, 'kp', ...
                %     'MarkerSize', 12, ...
                %     'MarkerFaceColor', 'k');
                xOpt
                xline(xOpt, 'k-', 'LineWidth', 2)

                % ---- Find intersections with other experiments ----
                for j = 1:size(Y,2)-1

                    % log10(y) = m*log10(x) + b
                    m = coeff(1,j);
                    b = coeff(2,j);

                    % Solve:
                    % log10(yOpt300) = m*log10(x) + b
                    logxIntersect = (log10(yOpt300) - b) / m;
                    xIntersect = 10^logxIntersect

                    % Only display intersections within plotted range
                    if xIntersect >= min(xfit) && xIntersect <= max(xfit)

                        % Plot intersection
                        % plot(xIntersect, yOpt300, 'rx', ...
                        %     'MarkerSize', 12, ...
                        %     'LineWidth', 2);
                        xline(xIntersect, 'r-', 'LineWidth', 2)

                        % Annotate number of cells
                        % text(xIntersect, yOpt300, ...
                        %     sprintf('  %.0f cells', xIntersect), ...
                        %     'FontSize', 10, ...
                        %     'FontWeight', 'bold', ...
                        %     'VerticalAlignment', 'bottom');
                    end
                end

                % Log axes
                set(gca, 'XScale', 'log', 'YScale', 'log');

                xlabel('Number of Cells');
                ylabel(names{k});
                title(names{k});

                grid on;
                legend('Location', 'best');
                hold off;

                ax = gca;
                ax.LineWidth = 1.5;
            end

        end
        function exportFigs3(obj,fignums,outputFolder,opts)
            arguments
                obj
                fignums
                outputFolder = 'AnnualReview_Figures';
                opts.fileType = 'svg';
            end

            if ~exist(outputFolder, 'dir')
                mkdir(outputFolder);
            end

            % Overall paper canvas
            fullWidth = 6.33;
            fullHeight = 6.33;

            % 4 x 3 grid
            plotWidth = fullWidth / 4;
            plotHeight = fullHeight / 4;

            for figNum = fignums

                fig = figure(figNum);

                % Remove figure-level title
                sgtitle(fig, '');

                % Find all axes
                axesList = findall(fig, 'Type', 'Axes');

                for i = 1:length(axesList)

                    ax = axesList(i);

                    % Remove title
                    ax.Title.String = '';
                    ax.Title.Visible = 'off';

                    % Remove axis labels
                    ax.XLabel.String = '';
                    ax.XLabel.Visible = 'off';

                    ax.YLabel.String = '';
                    ax.YLabel.Visible = 'off';

                    % Remove ticks and tick labels
                    % ax.XTick = [];
                    % ax.YTick = [];

                    ax.XTickLabel = [];
                    ax.YTickLabel = [];

                    % Remove tick marks
                    % ax.TickLength = [0 0];

                end

                % Remove legends
                legends = findall(fig, 'Type', 'Legend');

                if ~isempty(legends)
                    delete(legends);
                end

                % Remove colorbar labels/ticks
                colorbars = findall(fig, 'Type', 'ColorBar');

                for i = 1:length(colorbars)

                    cb = colorbars(i);

                    cb.TickLabels = [];
                    cb.Label.String = '';

                end

                % Set physical dimensions for 4 x 3 grid
                fig.Units = 'inches';
                fig.Position(3:4) = [plotWidth plotHeight];

                % Export
                fileName = sprintf('figure%d.%s', figNum.Number, opts.fileType);

                exportgraphics(fig, ...
                    fullfile(outputFolder, fileName), ...
                    'ContentType', 'vector');

            end

            disp('Figures exported successfully.');

        end

        %% Figure 4
        function obj = makeFig4ABC(obj,f4A,f4B,f4C)
            arguments
                obj
                f4A
                f4B
                f4C
            end
            Model_chg = obj.ModelKoffSig;

            % Solve and plot FSP without PDO effect.
            Model_chg.fspOptions.bounds = [];
            Model_chg.fspOptions.stateSpace = [];
            Model_chg = Model_chg.solve(solver='fsp');
            Model_chg.plotFSP(figureNums=f4A,plotType='marginals',indTimes=length(Model_chg.tSpan),speciesNames='mRNA',Colors={'k'})
            Model_chg.plotFSP(figureNums=f4C,plotType='marginals',indTimes=length(Model_chg.tSpan),speciesNames='mRNA',Colors={'k'}) %  lineProps={'LineWidth',2,'LineStyle','--'}

            % Create and plot Binomial PDO
            dropOut = 0.6; % fraction dropout
            obj.Model_BinomialPDO = Model_chg;
            obj.Model_BinomialPDO.pdoOptions.type = 'Binomial';
            obj.Model_BinomialPDO.pdoOptions.unobservedSpecies = 'gON';
            obj.Model_BinomialPDO.pdoOptions.props.CaptureProbabilityS1 = 0;    % Gene State is not measured
            obj.Model_BinomialPDO.pdoOptions.props.CaptureProbabilityS2 = 1-dropOut; % 95% dropout from RNA
            [~,obj.Model_BinomialPDO] = obj.Model_BinomialPDO.generatePDO(...
                showPlot=true, Title='Binomial PDO');
            fPDO = gcf;
            clf;
            copyobj(allchild(fPDO), f4B);
            close(fPDO);

            %Plot distributons with effect of PDO
            figure(f4C)
            hold on
            obj.Model_BinomialPDO.plotFSP(figureNums=f4C,plotType='marginals',indTimes=length(obj.Model_BinomialPDO.tSpan),...
                speciesNames='mRNA',includePDO=true,Colors={'r'})
        end

        function obj = prepareFig4DEF(obj,opts)
            arguments
                obj
                opts.nMLE = 20;
            end

            % First, generate the MLE scatter plot and FIM overlay (same as above).
            nMLE = opts.nMLE;

            MLE_noDistortion = obj.Model_BinomialPDO.estimateMLEspread(nCells=obj.nCellsInExperiment,...
                observableSpecies={'mRNA'},nMLE=nMLE,simsSaveFile='BurstFIMSimsPDO.csv',...
                freePars=obj.freeParsFig4,restart=true,useDistortions=false,correctDistortions=false,...
                nIter = 500);

            % Next, find MLE estimates WITHOUT correcting for the distortion.
            MLE_PDO_Uncorrected = obj.Model_BinomialPDO.estimateMLEspread(nCells=obj.nCellsInExperiment,...
                observableSpecies={'mRNA'},nMLE=nMLE,simsSaveFile='BurstFIMSimsPDO.csv',...
                freePars=obj.freeParsFig4,restart=false,useDistortions=true,correctDistortions=false,...
                nIter = 500);

            % Next, find MLE estimates with correcting for the distortion.
            % Model_chg = Model_chg.solve(solver='fsp');
            MLE_PDO_Corrected = obj.Model_BinomialPDO.estimateMLEspread(nCells=obj.nCellsInExperiment,...
                observableSpecies={'mRNA'},nMLE=nMLE,simsSaveFile='BurstFIMSimsPDO.csv',...
                freePars=obj.freeParsFig4,restart=false,useDistortions=true,correctDistortions=true,...
                nIter = 500);

            save('MLEforDistortions.mat', 'MLE_noDistortion', 'MLE_PDO_Uncorrected', 'MLE_PDO_Corrected')

        end

        function makeFigs4DEF(obj)
            load('MLEforDistortions.mat')

            %% Plot the spread of the mle and FIM estiamte
            % plot unaltered spread
            f1 = figure(304);

            fTrash = figure(350);

            Model_chg = obj.ModelKoffSig;

            FIMs = Model_chg.computeFIM(freePars=obj.freeParsFig4,scale='log');
            FIMTotal = Model_chg.totalFim(FIMs,obj.nCellsInExperiment);
            FIM = FIMTotal{1};

            MLElog = MLE_noDistortion.mhSamples/log(10);

            Model_chg.plotFIMResults(FIM^(-1)/log(10)^2, 'log',...
                Model_chg.parameters(1:5,1),...
                [Model_chg.parameters{1:5,2}],...
                PlotEllipses=true, ...
                EllipseFigure=f1,...
                EllipseLevel=0.95,...
                Colors=struct('EllipseColors',[0,0,0],'CenterSquare',[0,0,0]),...
                EllipsePairs=[1,2],...
                FigureHandle=fTrash,...
                LogThreshold=-4,...
                HeatMapType='invfim',...
                MatrixType='invfim');

            hold on

            scatter(MLElog(:,2),MLElog(:,1),10,[0.5 0.5 0.5],'filled');
            obj.plotMLEEllipse(MLElog(:,2),MLElog(:,1),0.95);


            f1 = figure(305);

            Model_chg.plotFIMResults(FIM^(-1)/log(10)^2, 'log',...
                Model_chg.parameters(1:5,1),...
                [Model_chg.parameters{1:5,2}],...
                PlotEllipses=true, ...
                EllipseLevel=0.95,...
                EllipseFigure=f1,...
                Colors=struct('EllipseColors',[0,0,0],'CenterSquare',[0,0,0]),...
                EllipsePairs=[3,4],...
                FigureHandle=fTrash,...
                LogThreshold=-4,...
                HeatMapType='invfim',...
                MatrixType='invfim');

            hold on

            scatter(MLElog(:,4),MLElog(:,3),10,[0.5 0.5 0.5],'filled');

            obj.plotMLEEllipse(MLElog(:,4),MLElog(:,3),0.95);

            % plot distorted

            f1 = figure(306);

            FIMs = Model_chg.computeFIM(freePars=obj.freeParsFig4,scale='log');
            FIMTotal = Model_chg.totalFim(FIMs,obj.nCellsInExperiment);
            FIM = FIMTotal{1};

            MLElog = MLE_PDO_Uncorrected.mhSamples/log(10);

            obj.Model_BinomialPDO.plotFIMResults(FIM^(-1)/log(10)^2, 'log',...
                Model_chg.parameters(1:5,1),...
                [Model_chg.parameters{1:5,2}],...
                PlotEllipses=true, ...
                EllipseFigure=f1,...
                EllipseLevel=0.95,...
                Colors=struct('EllipseColors',[0,0,0],'CenterSquare',[0,0,0]),...
                EllipsePairs=[1,2],...
                FigureHandle=fTrash,...
                LogThreshold=-4,...
                HeatMapType='invfim',...
                MatrixType='invfim');

            hold on

            scatter(MLElog(:,2),MLElog(:,1),10,[0.5 0.5 0.5],'filled');
            obj.plotMLEEllipse(MLElog(:,2),MLElog(:,1),0.95);


            f1 = figure(307);

            obj.Model_BinomialPDO.plotFIMResults(FIM^(-1)/log(10)^2, 'log',...
                Model_chg.parameters(1:5,1),...
                [Model_chg.parameters{1:5,2}],...
                PlotEllipses=true, ...
                EllipseFigure=f1,...
                EllipseLevel=0.95,...
                Colors=struct('EllipseColors',[0,0,0],'CenterSquare',[0,0,0]),...
                EllipsePairs=[3,4],...
                FigureHandle=fTrash,...
                LogThreshold=-4,...
                HeatMapType='invfim',...
                MatrixType='invfim');

            hold on

            scatter(MLElog(:,4),MLElog(:,3),10,[0.5 0.5 0.5],'filled');
            obj.plotMLEEllipse(MLElog(:,4),MLElog(:,3),0.95);

            % plot corrected distorted

            f1 = figure(308);

            FIMs = obj.Model_BinomialPDO.computeFIM(freePars=obj.freeParsFig4,scale='log');
            FIMTotal = obj.Model_BinomialPDO.totalFim(FIMs,obj.nCellsInExperiment);
            FIM = FIMTotal{1};

            MLElog = MLE_PDO_Corrected.mhSamples/log(10);

            obj.Model_BinomialPDO.plotFIMResults(FIM^(-1)/log(10)^2, 'log',...
                Model_chg.parameters(1:5,1),...
                [Model_chg.parameters{1:5,2}],...
                PlotEllipses=true, ...
                EllipseFigure=f1,...
                EllipseLevel=0.95,...
                Colors=struct('EllipseColors',[0,0,0],'CenterSquare',[0,0,0]),...
                EllipsePairs=[1,2],...
                FigureHandle=fTrash,...
                LogThreshold=-4,...
                HeatMapType='invfim',...
                MatrixType='invfim');

            hold on

            scatter(MLElog(:,2),MLElog(:,1),10,[0.5 0.5 0.5],'filled');
            obj.plotMLEEllipse(MLElog(:,2),MLElog(:,1),0.95);


            f1 = figure(309);

            Model_chg.plotFIMResults(FIM^(-1)/log(10)^2, 'log',...
                Model_chg.parameters(1:5,1),...
                [Model_chg.parameters{1:5,2}],...
                PlotEllipses=true, ...
                EllipseLevel=0.95,...
                EllipseFigure=f1,...
                Colors=struct('EllipseColors',[0,0,0],'CenterSquare',[0,0,0]),...
                EllipsePairs=[3,4],...
                FigureHandle=fTrash,...
                LogThreshold=-4,...
                HeatMapType='invfim',...
                MatrixType='invfim');

            hold on

            scatter(MLElog(:,4),MLElog(:,3),10,[0.5 0.5 0.5],'filled');
            obj.plotMLEEllipse(MLElog(:,4),MLElog(:,3),0.95);
        end
        function makeFig4GHI(obj)
            % In this section, we compute the FIM for different dropout fractions.  The
            % current analysis only allows for a single define experiment (i.e., the
            % change from a pre-specified S0 to a pre-specified S1). The experiment
            % design option is to decide on the time points at which to take the
            % observations and how masny cells to observe at each time point.           
            N = 50;
            vDropOut = linspace(0,0.98,N);
            OptExptVsDropOut = zeros(50,length(Model_chg.tSpan));
            ModelPDO = obj.ModelKoffSig;
            ModelPDO = ModelPDO.solve(solver='fspsens');
            ModelPDO.pdoOptions.type = 'Binomial';
            ModelPDO.pdoOptions.unobservedSpecies = 'gON';
            TotalFim = cell(N,1);
            detFIMOrig = zeros(N,1);
            detFIMOpt = zeros(N,1);
            nCellsOrig = zeros(N,1);
            nCellsOpt = zeros(N,1);
            NCellsTotal = sum(obj.nCellsInExperiment);
            for i = 1:N
                dropOut = vDropOut(i);
                ModelPDO.pdoOptions.props.CaptureProbabilityS1 = 0;    % Gene State is not measured
                ModelPDO.pdoOptions.props.CaptureProbabilityS2 = 1-dropOut; % 95% dropout from RNA
                [~,ModelPDO] = ModelPDO.generatePDO;
                FIMs = ModelPDO.computeFIM(scale='log',freePars=obj.freeParsFig4,...
                    observed={'mRNA'});
                OptExperiment(i,:) = ModelPDO.optimizeCellCounts(FIMs,NCellsTotal,'D-opt');
                TotalFimOrig(i,1) = ModelPDO.totalFim(FIMs,obj.nCellsInExperiment);
                TotalFimOpt(i,1) = ModelPDO.totalFim(FIMs,OptExperiment(i,:));
                detFIMOrig(i) = det(TotalFimOrig{i,1});
                detFIMOpt(i) = det(TotalFimOpt{i,1});
                nCellsOrig(i) = NCellsTotal*(detFIMOrig(1)/detFIMOrig(i))^(1/4);
                nCellsOpt(i) = NCellsTotal*(detFIMOrig(1)/detFIMOpt(i))^(1/4);
            end

            % Plot the determinant of the inverse FIM versus the drop out rate
            figure(310); clf;
            plot(vDropOut,1./detFIMOrig,'b',vDropOut,1./detFIMOpt,'r--','linewidth',3)
            set(gca,'yscale','log')
            xlabel('Drop Out Fraction')
            ylabel('Det(FIM^{-1})')

            % Plot the nmber of cells that need to be measured to achieve the same
            % information (same expected determinant of FIM) as was achieved when we
            % did the original experiment design with 600 cells.
            figure(311); clf;
            plot(vDropOut,nCellsOrig,'b',vDropOut,nCellsOpt,'r--','linewidth',3)
            set(gca,'yscale','log')
            xlabel('Drop Out Fraction')
            ylabel('Required Number of Cells')

            figure(312); clf;
            % Plot the optimal experiment design versus the dropout rate, constrained
            % to have the same original number of cells (600).  In this plot, the
            % colors will represent the fraction of cells that are measure at each time
            % point.
            pcolor(vDropOut,[ModelPDO.tSpan,ModelPDO.tSpan(end)+(ModelPDO.tSpan(end)-ModelPDO.tSpan(end-1))],[OptExperiment,zeros(N,1)]'/600)
            % set(gca,'yscale','log')
            ylabel('Measurement Time')
            xlabel('Drop Out Fraction')
            c = colorbar;
            c.Label.String = 'Fraction of Cells'

        end

    end
    methods (Static)
        %% Functions
        function plotHeatmap(M,rowNames,colNames,titleText,logThreshold)
            % plotHeatmap
            %
            % Standalone heatmap using the same log transformation and colormap
            % as plotFIMResults.
            %
            % Infinite values are displayed at the maximum color scale and labeled
            % with +inf / -inf rather than +/-.
            %
            %   M             : matrix to display
            %   rowNames      : row names
            %   colNames      : column names
            %   titleText     : optional title
            %   logThreshold  : optional log10 threshold, default = -5
            %
            % Example:
            %
            %   plotHeatmap(M,params,params,'FIM',-4)

            if nargin < 4
                titleText = '';
            end

            if nargin < 5
                logThreshold = -5;
            end

            % -------------------------------------------------------------
            % Check dimensions
            % -------------------------------------------------------------

            if size(M,1) ~= numel(rowNames)
                error('Number of row names must equal number of rows.');
            end

            if size(M,2) ~= numel(colNames)
                error('Number of column names must equal number of columns.');
            end

            rowNames = cellstr(rowNames);
            colNames = cellstr(colNames);

            % -------------------------------------------------------------
            % Log transformation
            % -------------------------------------------------------------

            threshold = 10^logThreshold;

            isPosInf = isinf(M) & M > 0;
            isNegInf = isinf(M) & M < 0;
            isInf    = isPosInf | isNegInf;

            % -------------------------------------------------------------
            % Transform finite values only
            % -------------------------------------------------------------

            Mfinite = M;
            Mfinite(isInf) = NaN;

            posValues = ...
                1*(Mfinite >= threshold) + ...
                -1*(Mfinite <= -threshold);

            logMag = ...
                max(0,log10(abs(Mfinite))-logThreshold) .* posValues;

            fimDisp = logMag;

            % -------------------------------------------------------------
            % Determine maximum finite display value
            % -------------------------------------------------------------

            finiteDisp = fimDisp(isfinite(fimDisp));

            if isempty(finiteDisp)

                % Entire matrix is infinite
                x1 = 1;

            else

                x1 = max(abs(finiteDisp),[],'all');

                if x1 == 0
                    x1 = 1;
                end

            end

            % -------------------------------------------------------------
            % Make the color scale symmetric and integer-bounded
            %
            % This is important because the colorbar ticks and the colormap
            % must use exactly the same endpoints.
            % -------------------------------------------------------------

            colorMax = ceil(x1);

            if colorMax == 0
                colorMax = 1;
            end

            % -------------------------------------------------------------
            % Assign Inf values directly to the color-scale endpoints
            %
            % +Inf -> +colorMax -> darkest red
            % -Inf -> -colorMax -> darkest blue
            % -------------------------------------------------------------

            fimDisp(isPosInf) = colorMax;
            fimDisp(isNegInf) = -colorMax;

            % -------------------------------------------------------------
            % Color range
            % -------------------------------------------------------------

            rangeColors = [-colorMax, 0, colorMax];

            % -------------------------------------------------------------
            % Colorbar ticks
            % -------------------------------------------------------------

            cbTicks = [ ...
                floor(rangeColors(1)):0, ...
                1:ceil(rangeColors(3))];

            % Remove any duplicate zero
            cbTicks = unique(cbTicks,'stable');

            % -------------------------------------------------------------
            % Colorbar labels
            % -------------------------------------------------------------

            if logThreshold < 0

                cbTickLabels = cell(size(cbTicks));

                for k = 1:numel(cbTicks)

                    v = cbTicks(k);

                    if v < 0

                        cbTickLabels{k} = sprintf( ...
                            '$-10^{%g}$', ...
                            -v + logThreshold);

                    elseif v == 0

                        cbTickLabels{k} = ...
                            ['$\pm 10^{',num2str(logThreshold),'}$'];

                    else

                        cbTickLabels{k} = sprintf( ...
                            '$10^{%g}$', ...
                            v + logThreshold);

                    end

                end

            elseif logThreshold > 0

                cbTickLabels = cell(size(cbTicks));

                for k = 1:numel(cbTicks)

                    v = cbTicks(k);

                    if v < 0

                        cbTickLabels{k} = sprintf( ...
                            '$-10^{%g}$', ...
                            -v + logThreshold);

                    elseif v == 0

                        cbTickLabels{k} = ...
                            ['$\pm 10^{',num2str(logThreshold),'}$'];

                    else

                        cbTickLabels{k} = sprintf( ...
                            '$10^{%g}$', ...
                            v + logThreshold);

                    end

                end

            else

                cbTickLabels = cell(size(cbTicks));

                for k = 1:numel(cbTicks)

                    v = cbTicks(k);

                    if v < 0

                        cbTickLabels{k} = sprintf( ...
                            '$-10^{%g}$', ...
                            -v);

                    elseif v == 0

                        cbTickLabels{k} = '$0$';

                    else

                        cbTickLabels{k} = sprintf( ...
                            '$10^{%g}$', ...
                            v);

                    end

                end

            end

            % -------------------------------------------------------------
            % Colormap
            %
            % The endpoints are now EXACTLY:
            %
            %   -colorMax -> dark blue
            %   0         -> white
            %   +colorMax -> dark red
            %
            % -------------------------------------------------------------

            cmap = Forman2027.blueWhiteFlatRed( ...
                -colorMax, ...
                0, ...
                0, ...
                colorMax);

            % -------------------------------------------------------------
            % Plot
            % -------------------------------------------------------------

            imagesc(fimDisp);

            ax = gca;

            axis square;

            colormap(ax,cmap);

            % IMPORTANT:
            % Explicitly use the same endpoints as the colormap.
            clim([-colorMax,colorMax]);

            % -------------------------------------------------------------
            % Colorbar
            % -------------------------------------------------------------

            cb = colorbar;

            cb.Ticks = cbTicks;
            cb.TickLabels = cbTickLabels;
            cb.TickLabelInterpreter = 'latex';

            cb.Label.String = 'Value';

            % -------------------------------------------------------------
            % Axis formatting
            % -------------------------------------------------------------

            ax.XTick = 1:numel(colNames);
            ax.YTick = 1:numel(rowNames);

            ax.XTickLabel = colNames;
            ax.YTickLabel = rowNames;

            ax.FontSize = 11;
            ax.LineWidth = 0.5;
            ax.TickDir = 'out';
            ax.Box = 'on';

            xlabel('Parameter');
            ylabel('Parameter');

            % -------------------------------------------------------------
            % Cell boundaries
            % -------------------------------------------------------------

            hold on;

            for x = 0.5:1:size(M,2)+0.5

                plot([x x], ...
                    [0.5 size(M,1)+0.5], ...
                    'Color',[0.5 0.5 0.5], ...
                    'LineWidth',0.5);

            end

            for y = 0.5:1:size(M,1)+0.5

                plot([0.5 size(M,2)+0.5], ...
                    [y y], ...
                    'Color',[0.5 0.5 0.5], ...
                    'LineWidth',0.5);

            end

            % -------------------------------------------------------------
            % Overlay symbols
            % -------------------------------------------------------------

            for i = 1:size(M,1)

                for j = 1:size(M,2)

                    % -----------------------------------------------------
                    % Infinite values
                    % -----------------------------------------------------

                    if isPosInf(i,j)

                        txt = '+\infty';
                        fontSize = 5;

                    elseif isNegInf(i,j)

                        txt = '-\infty';
                        fontSize = 5;

                        % -----------------------------------------------------
                        % Normal finite values
                        % -----------------------------------------------------

                    elseif M(i,j) >= threshold

                        txt = '+';
                        fontSize = 11;

                    elseif M(i,j) <= -threshold

                        txt = '-';
                        fontSize = 11;

                    else

                        continue;

                    end

                    text(j,i,txt, ...
                        'HorizontalAlignment','center', ...
                        'VerticalAlignment','middle', ...
                        'FontWeight','bold', ...
                        'FontSize',fontSize, ...
                        'Color','k', ...
                        'Interpreter','tex');

                end

            end

            hold off;

            % -------------------------------------------------------------
            % Title
            % -------------------------------------------------------------

            if ~isempty(titleText)

                title(titleText, ...
                    'FontSize',14, ...
                    'FontWeight','normal');

            end

        end

        % =============================================================
        % Blue -> white -> red
        % =============================================================

        function cmap = blueWhiteFlatRed(x1,x2,x3,x4,n)

            if nargin < 5
                n = 256;
            end

            xs = linspace(x1,x4,n);

            blue  = [0 0 0.6];
            white = [1 1 1];
            red   = [0.6 0 0];

            cmap = zeros(n,3);

            for i = 1:n

                x = xs(i);

                if x <= x2

                    t = (x-x1)/(x2-x1);

                    cmap(i,:) = ...
                        (1-t)*blue + t*white;

                elseif x <= x3

                    cmap(i,:) = white;

                else

                    t = (x-x3)/(x4-x3);

                    cmap(i,:) = ...
                        (1-t)*white + t*red;

                end


            end

        end

        function h = plotMLEEllipse(x,y,confidence)

            % Remove invalid samples
            valid = isfinite(x) & isfinite(y);
            x = x(valid);
            y = y(valid);

            % Mean of MLE samples
            mu = [mean(x), mean(y)];

            % Empirical covariance
            C = cov([x y]);

            % Eigenvectors/eigenvalues of covariance matrix
            [V,D] = eig(C);

            % Sort eigenvalues from largest to smallest
            [lambda,idx] = sort(diag(D),'descend');
            V = V(:,idx);

            % Chi-square scaling for a 2D confidence ellipse
            scale = sqrt(chi2inv(confidence,2));

            % Parametric ellipse
            theta = linspace(0,2*pi,300);

            ellipse = V * diag(sqrt(lambda)) * scale * ...
                [cos(theta); sin(theta)];

            ellipse(1,:) = ellipse(1,:) + mu(1);
            ellipse(2,:) = ellipse(2,:) + mu(2);

            % Plot ellipse
            h = plot(ellipse(1,:),ellipse(2,:),...
                'c-',...
                'LineWidth',2);

        end
    end

end
