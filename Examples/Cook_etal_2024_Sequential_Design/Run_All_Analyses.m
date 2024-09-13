% This script re-generated all figures associated with the manuscript:
%    J Cook, E Ron, D Svetlov, LU Aguilera, B Munsky, "Sequential design of
%    single-cell experiments to identify discrete stochastic models for
%    gene expression."
%
%%
close all
clear all
addpath(genpath('../../src'));
addpath('codes')

%% Specific which figures to generate and then run script.
Figure_to_Generate = 1;
switch Figure_to_Generate
    case 1  % Figure 1.
        %%   Generate (or load saved) results
        [TestCasesBFBD,finalExperimentDesignBFBD] = sequentialExptDesignBatchRunner('Poisson',1,1,true,false,[4]);
        [TestCasesRD,finalExperimentDesignRD] = sequentialExptDesignBatchRunner('Poisson',2,1,true,false,[4]);
        close(1,2,3,4,103)
        %%   Reformat Figure 1B Results (Convergence of |COV|).
        f = figure(101);
        set(f,'Position',[ 616   659   429   259])
        set(gca,'ylim',[1e-9,1e-4],'xlim',[1.5,8],'XTick',[2:8])
        grid on

        legs = {'BFBD: MHA','BFBD: FIM$_{\rm Est}^{-1}$','BFBD: FIM$_{\rm True}^{-1}$',...
            'RD: MHA','RD: FIM$_{\rm Est}^{-1}$','RD: FIM$_{\rm True}^{-1}$'};
        legend(legs,'Interpreter','latex')
        xlabel('Experiment Round','Interpreter','latex')
        ylabel('Determinant of Covariance','Interpreter','latex')
        %%   Print |COV| values at specific rounds.
        clc
        disp(['|COV| after 8 rounds for RD = ',num2str(TestCasesRD.detCov(8))])
        disp(['|COV| after 5 rounds for BFBD = ',num2str(TestCasesBFBD.detCov(5))])

        %%   Figure 1D, Top (Optimized Experiment Design)
        f2 = figure;
        set(f2,'Position',[939   687   404   238])
        colormap('sky')
        sz = size(finalExperimentDesignBFBD);
        pCol = finalExperimentDesignBFBD; pCol(end+1,end+1) = 0;
        pcolor([-0.25:0.5:10.25],[0.5:10.5],pCol)
        cb = colorbar;
        cb.Label.FontSize = 16;
        cb.Label.String = 'Number of Cells';
        cb.Label.Interpreter = 'latex';
        set(gca,'FontSize',16,'TickLabelInterpreter','latex');
        xlabel('Measurement Time (min)','Interpreter','latex')
        ylabel('Input ($\mu$M)','Interpreter','latex')

        %%   Figure 1D, Bottom (Random Experiment Design)
        f3 = figure;
        set(f3,'Position',[939   687   404   238])
        colormap('sky')
        sz = size(finalExperimentDesignRD);
        pCol = finalExperimentDesignRD; pCol(end+1,end+1) = 0;
        pcolor([-0.25:0.5:10.25],[0.5:10.5],pCol)
        cb = colorbar;
        cb.Label.String = 'Number of Cells';
        cb.Label.Interpreter = 'latex';
        cb.Label.FontSize = 16;
        set(gca,'FontSize',16,'TickLabelInterpreter','latex');
        xlabel('Measurement Time (min)','Interpreter','latex')
        ylabel('Input ($\mu$M)','Interpreter','latex')


        %%   Figure 1C (Scatter plots MHA, Round 4)
        fs = [122,130];
        for i = fs
            f122 = figure(i);

            set(f122,'Position',[1098         576         560         420])
            subplot(2,2,1);
            lx1 = [0.9 1.25];
            lx2 = [-0.6 -0.35];
            lx3 = [0.6 0.95];
            set(gca,'FontSize',16,'TickLabelInterpreter','latex','xlim',lx2,'ylim',lx1);
            ylabel('$\log_{10} (k_r)$','Interpreter','latex')
            xlabel('$\log_{10} (\gamma)$','Interpreter','latex')

            subplot(2,2,2);
            set(gca,'FontSize',16,'TickLabelInterpreter','latex','ylim',lx1,'xlim',lx3);
            xlabel('$\log_{10} (k_D)$','Interpreter','latex')

            subplot(2,2,4);
            set(gca,'FontSize',16,'TickLabelInterpreter','latex','ylim',lx2,'xlim',lx3);
            ylabel('$\log_{10} (\gamma)$','Interpreter','latex')
            xlabel('$\log_{10} (k_D)$','Interpreter','latex')
        end

        %%
    case 2  % Figure 2
        clear all
        close all

        %%   Generate or load data
        [TestCasesBFBD,finalExperimentDesignBFBD] = sequentialExptDesignBatchRunner('Burst',1,2,true,false,[4]);
        [TestCasesRD,finalExperimentDesignRD] = sequentialExptDesignBatchRunner('Burst',2,2,true,false,[4]);
        close(1,2,3,4,503)
        %%   Figure 2B (Convergence of |COV|)
        f = figure(501);
        set(f,'Position',[ 616   659   429   259])
        set(gca,'xlim',[2,8],'XTick',[2:8],'YLim',[1e-13,1e-8])
        grid on

        legs = {'BFBD: MHA','BFBD: FIM$_{\rm Est}^{-1}$','BFBD: FIM$_{\rm True}^{-1}$',...
            'RD: MHA','RD: FIM$_{\rm Est}^{-1}$','RD: FIM$_{\rm True}^{-1}$'};
        legend(legs,'Interpreter','latex')
        xlabel('Experiment round','Interpreter','latex')
        ylabel('Determinant of Covariance','Interpreter','latex')

        %%   Figure 2D, Top (Optimized Experiment Design)
        f2 = figure;
        set(f2,'Position',[939   687   404   238])
        colormap('sky')
        sz = size(finalExperimentDesignBFBD);
        pCol = finalExperimentDesignBFBD; pCol(end+1,end+1) = 0;
        pcolor([-0.5:1:30.5],[0.5:10.5],pCol)
        cb = colorbar;
        cb.Label.FontSize = 16;
        cb.Label.String = 'Number of Cells';
        cb.Label.Interpreter = 'latex';
        set(gca,'FontSize',16,'TickLabelInterpreter','latex');
        xlabel('Measurement Time (min)','Interpreter','latex')
        ylabel('Input ($\mu$M)','Interpreter','latex')

        %%   Figure 1D, Bottom (Random Experiment Design)
        f3 = figure;
        set(f3,'Position',[939   687   404   238])
        colormap('sky')
        sz = size(finalExperimentDesignRD);
        pCol = finalExperimentDesignRD; pCol(end+1,end+1) = 0;
        pcolor([-0.5:1:30.5],[0.5:10.5],pCol)
        cb = colorbar;
        cb.Label.FontSize = 16;
        cb.Label.String = 'Number of Cells';
        cb.Label.Interpreter = 'latex';
        set(gca,'FontSize',16,'TickLabelInterpreter','latex');
        xlabel('Measurement Time (min)','Interpreter','latex')
        ylabel('Input ($\mu$M)','Interpreter','latex')

        %%   Figure 2C (Scatter plots MHA, Round 4)
        fs = [522,530];
        for i = fs
            f122 = figure(i);

            set(f122,'Position',[389         143        1155         852])
            subplot(5,5,9);
            lx1 = [-1.2 -0.7];
            lx2 = [-0.9 -0.3];
            lx3 = [0.8 1.2];
            lx4 = [-0.7 -0.35];
            lx5 = [0 1];
            set(gca,'FontSize',16,'TickLabelInterpreter','latex','xlim',lx5,'ylim',lx2);
            ylabel('$\log_{10} (k_{\rm off}^{0})$','Interpreter','latex')
            xlabel('$\log_{10} (k_{\rm D})$','Interpreter','latex')

            subplot(5,5,13);
            set(gca,'FontSize',16,'TickLabelInterpreter','latex','ylim',lx3,'xlim',lx4);
            ylabel('$\log_{10} (k_{\rm r})$','Interpreter','latex')
            xlabel('$\log_{10} (\gamma_{\rm r})$','Interpreter','latex')

        end

        %%
    case 3  % Figure 3
        clear all
        close all
        %%   Generate or load data
        [TestCasesBFBD,finalExperimentDesignBFBD] = sequentialExptDesignBatchRunner('GR',1,2,true,false,[4]);
        [TestCasesRD,finalExperimentDesignRD] = sequentialExptDesignBatchRunner('GR',2,2,true,false,[4]);
        
        %%   Figure 3B (Convergence of |COV|)
        f = figure(801);
        set(f,'Position',[ 616   659   429   259])
        set(gca,'xlim',[1,6],'XTick',[2:6],'YLim',10.^[-8.5,-6])
        grid on
        legs = {'BFBD: MHA','BFBD: FIM$_{\rm Pre}^{-1}$','BFBD: FIM$_{\rm Post}^{-1}$',...
            'RD: MHA','RD: FIM$_{\rm Pre}^{-1}$','RD: FIM$_{\rm Post}^{-1}$'};

        Ch = get(gca,'Children');
        set(Ch(5),'Marker','o','MarkerSize',20);
        set(Ch(2),'Marker','o','MarkerSize',20);

        legend(legs,'Interpreter','latex','FontSize',13)
        xlabel('Experiment round','Interpreter','latex')
        ylabel('Determinant of Covariance','Interpreter','latex')

        %%   Figure 2E,(Optimized Experiment Design)
        f2 = figure;
        set(f2,'Position',[939   687   404   238])
        colormap('sky')
        sz = size(finalExperimentDesignBFBD);
        pCol = finalExperimentDesignBFBD; pCol(end+1,end+1) = 0;
        pcolor([-0:1:6],[0:3],pCol)
        cb = colorbar;
        cb.Label.FontSize = 16;
        cb.Label.String = 'Number of Cells';
        cb.Label.Interpreter = 'latex';
        set(gca,'FontSize',16,'TickLabelInterpreter','latex',...
            'ytick',[0.5:1:3.5],'YTickLabel',[1,10,100],...
            'xtick',[0.5:1:5.5],'XTickLabel',[0,10,30,50,75,180]);
        xlabel('Measurement Time (min)','Interpreter','latex')
        ylabel('Input ($\mu$M)','Interpreter','latex')

        %%   Figure 3C (Scatter plots MHA, Round 4)
        fs = [822];
        for i = fs
            f122 = figure(i);

            set(f122,'Position',[389         143        1155         852])
            lx{1} = [-2.4 -2.1];
            lx{2} = [1.1, 1.25];
            lx{3} = [-2 -1.7];
            lx{4} = [-1.95 -1.75];
            lx{5} = [-2.35 -2.2];
            lx{6} = [0.9 1.1];

            labs = {'k_{\rm cn}^0';'k_{\rm cn}^{1}';'k_{\rm nc}';'k';'\gamma_{\rm nuc}';'M_{\rm Dex}'};

            for ii = 1:5
                for jj = ii:5
                    k = (ii-1)*5+jj;
                    subplot(5,5,k);
                    set(gca,'FontSize',16,'TickLabelInterpreter','latex','xlim',lx{jj+1},'ylim',lx{ii});
                    ylabel('$\log_{10} (k_{\rm off}^{0})$','Interpreter','latex')
                    ylabel(['$',labs{ii},'$'],'Interpreter','latex')
                    xlabel(['$',labs{jj+1},'$'],'Interpreter','latex')
                end
            end
        end


        %%   Figure 3F (Make Predictions)
        %%       Solve Predictions Model
        ModelTrue = SSIT;
        ModelTrue.species = {'cytGR';'nucGR'};
        ModelTrue.initialCondition = [20;1];
        ModelTrue.fspOptions.bounds = [0,0,30,30];
        ModelTrue.fspOptions.verbose = false;
        ModelTrue.fspOptions.fspIntegratorAbsTol = 1e-10;
        ModelTrue.propensityFunctions = {'kcn0*(1 + (t>0)*kcn1*IDex/(MDex+IDex)) * cytGR';...
            'knc*nucGR'; 'kg1';'gnuc*nucGR'};
        ModelTrue.stoichiometry = [-1,1,1,0;...
            1,-1,0,-1];
        ModelTrue.customConstraintFuns = {'x1+x2'};

        ModelTrue.parameters = {'kcn0',1;'kcn1',1;'knc',1;...
            'kg1',1;'gnuc',1;'MDex',1};

        % Update to fitted parameters.
        round4Predictions = 6;
        ModelTrue.parameters(:,2) = num2cell(TestCasesBFBD.pars(round4Predictions,:));
        ModelTrue.inputExpressions = {'IDex','10'};

        ModelTrue.fspOptions.initApproxSS = true;
        ModelTrue.fspOptions.usePiecewiseFSP = true;
        ModelTrue.fspOptions.constantJacobian = true;
        ModelTrue.fspOptions.constantJacobianTime = 0.1;
        ModelTrue = ModelTrue.formPropensitiesGeneral(['Predictor'],true);

        for i=1:3
            ModelTrue.fspOptions.fspTol = 1e-8;
            [fspSoln,ModelTrue.fspOptions.bounds] = ModelTrue.solve;
        end

        ModelTrue = ModelTrue.loadData("ExampleData/Data.csv",...
            {'nucGR','normgrnuc';'cytGR','normgrcyt'},...
            {'Dex_Conc','10'});

        %%       Make  figures.
        close(1,2,3,4,803,830)
        ModelTrue.makeFitPlot([],1);
        %%       Reformat Fig 3F,Bottom (Nuclear Distributions)
        origFigs = [2];
        titles = {'Nuc GR 10nM'};
        newFigs = [1022];

        for idex = 1
            oldFig = figure(origFigs(idex));
            newFig = figure(newFigs(idex));clf
            set(newFig,'Position',[1000        1039        1361         199])

            set(newFig,'Name',titles{idex})
            % Get handles of all subplots in the original figure
            subplotHandles = findobj(oldFig, 'type', 'axes');

            % Iterate over each subplot and copy its contents to the new figure
            orderset = [7,6,5,4,3,2,1];
            for i2 = 1:length(orderset)
                i = orderset(i2);
                subplotHandle = subplotHandles(i);

                % Create new subplot in the new figure
                ax = subplot(1, 7, i2, 'Parent', newFig);

                % Copy contents of the original subplot to the new subplot
                copyobj(allchild(subplotHandle), ax);

                if i2~=1
                    ylabel('')
                    set(gca,'yticklabels',[])
                end

                h = gca;
                set(h,'Children',h.Children([1,2,3,4]))

                grid on

                h.Children(3).Visible = 'off';
                h.Children(4).Visible = 'off';

                h.Children(1).LineWidth = 4;
                h.Children(2).LineWidth = 4;

                % Adjust concentration scale for Nuclear GR.
                ratio = 0.520900129552663; %Fit in "ProcessEric_Data_Feb6".
                for ich=1:2
                    h.Children(ich).XData = h.Children(ich).XData/ratio;
                end

                set(gca,'xlim',[0,40],'ylim',[0,0.3],'FontSize',16)

            end

        end

        %%       Reformat Fig 3F,Top (Cytoplasmic Distributions
        origFigs = [2];
        titles = {'Cyt GR 10nM'};
        newFigs = [1032];

        for idex = 1:1
            oldFig = figure(origFigs(idex));
            newFig = figure(newFigs(idex));clf
            set(newFig,'Position',[1000        1039        1361         199])

            set(newFig,'Name',titles{idex})
            % Get handles of all subplots in the original figure
            subplotHandles = findobj(oldFig, 'type', 'axes');

            % Iterate over each subplot and copy its contents to the new figure
            orderset = [7,6,5,4,3,2,1];
            for i2 = 1:length(orderset)
                i = orderset(i2);
                subplotHandle = subplotHandles(i);

                % Create new subplot in the new figure
                ax = subplot(1, 7, i2, 'Parent', newFig);

                % Copy contents of the original subplot to the new subplot
                copyobj(allchild(subplotHandle), ax);

                if i2~=1
                    ylabel('')
                    set(gca,'yticklabels',[])
                end

                h = gca;
                set(h,'Children',h.Children([1,2,3,4]))

                grid on

                h.Children(3).LineWidth = 4;
                h.Children(4).LineWidth = 4;

                h.Children(1).Visible = 'off';
                h.Children(2).Visible = 'off';

                set(gca,'xlim',[0,20],'ylim',[0,0.5],'FontSize',16)

            end
        end
        close(1,2,3,4)

end

