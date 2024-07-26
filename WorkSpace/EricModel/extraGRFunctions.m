%% Extra Functions
classdef extraGRFunctions
    properties
        
    end
    methods
        function makeGRPlots(combinedModel,GRpars)
        combinedGRModel = combinedModel.updateModels(GRpars,false);
        nMods = length(combinedGRModel.SSITModels);
        ModelGroup = cell(nMods,1);
        for i=1:nMods
            %  Update parameters in original models.
            ModelGroup{i} = combinedGRModel.SSITModels{i};
            ModelGroup{i}.tSpan = sort(unique([ModelGroup{i}.tSpan,linspace(0,180,30)]));
            ModelGroup{i}.makeFitPlot([],1,[],true,'STD')
        end
        end
        function plotODEresults(extendedMod,soln,modeWithGRData,fignum)
        arguments
            extendedMod
            soln
            modeWithGRData
            fignum = 1;
        end
        figure(fignum); clf;
        % Plot GR levels vs. Time
        subplot(2,1,1)
        plot(extendedMod.tSpan,soln.ode(:,3:4),'--','LineWidth',2);hold on
        plot(modeWithGRData.dataSet.times,modeWithGRData.dataSet.mean,'s','MarkerSize',16,'MarkerFaceColor','k','LineWidth',3)
        legend(extendedMod.species(3:4))
        set(gca,'xlim',[-10,200],'ylim',[0,12],'fontsize',16)
        ylabel('GR Concentrations (UA)')
        legend({'Cyt-Model','Nuc-Model','Cyt-Data','Nuc-Data'})
        title('GR')

        % Plot DUSP1 levels vs. Time
        subplot(2,1,2)
        plot(extendedMod.tSpan,soln.ode(:,5:6),'--','LineWidth',2);hold on
        plot(extendedMod.dataSet.times,extendedMod.dataSet.mean,'s','MarkerSize',16,'MarkerFaceColor','k','LineWidth',3)
        legend(extendedMod.species(3:4))
        set(gca,'xlim',[-10,200],'ylim',[0,160],'fontsize',16)
        ylabel('GR Concentrations (UA)')
        legend({'Nuc-Model','Cyt-Model','Nuc-Data','Cyt-Data'})
        title('DUSP1')
        end

        function makeCytDistPlots(ssaSoln_100,extendedMod,fignum,timeIndsMod,timeIndsDat,speciesIndMod,speciesIndDat)
        arguments
            ssaSoln_100
            extendedMod
            fignum = 1;
            timeIndsMod = [];
            timeIndsDat = [];
            speciesIndMod = 1;
            speciesIndDat = 1;
        end
        figure(fignum); clf;
        if isempty(timeIndsMod)
            timeIndsMod = [1:length(ssaSoln_100.T_array)];
        end
        if isempty(timeIndsDat)
            timeIndsDat = [1:length(extendedMod.dataSet.times)];
        end

        if length(timeIndsMod)~=length(timeIndsDat)
            error('Length of data and model time points must be equal')
        end

        Nrows = ceil(sqrt(length(timeIndsMod)));
        for i = 1:length(timeIndsMod)
            subplot(Nrows,Nrows,i)

            % Add SSA to histogram plot
            M = squeeze(ssaSoln_100.trajs(speciesIndMod,timeIndsMod(i),:));
            H = histogram(M,'Normalization','pdf');
            hold on

            % Add data to histogram plot
            dMat = double(extendedMod.dataSet.app.DataLoadingAndFittingTabOutputs.dataTensor(timeIndsDat(i),:,:));
            N = sum(dMat,"all");
            if speciesIndDat==1
                PD = [0;sum(dMat,2)/N];
            else
                PD = [0;sum(dMat,1)'/N];
            end
    
            binEdges = round(H.BinEdges);
            nBins = length(binEdges)-1;
            PDbinned = zeros(nBins,1);
            binwidth = binEdges(2)-binEdges(1);
            for j = 1:nBins
                PDbinned(j) = sum(PD(binEdges(j)+1:binEdges(j+1)));
            end

            PDbinned = [PDbinned;1-sum(PDbinned)];
            stairs(binEdges,PDbinned/binwidth,'linewidth',2)

            set(gca,'FontSize',15,'ylim',[0,0.03])
    
        end
        end
    end
end