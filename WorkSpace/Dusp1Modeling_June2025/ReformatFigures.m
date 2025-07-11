splitReps = true;

%% To generate GR Nuc plots 
origFigs = [3,7,11];
titles = {'Nuc GR 1nM','Nuc GR 10nM','Nuc GR 100nM'};
newFigs = [1001,1002,1003];

for idex = 1:3
    oldFig = figure(origFigs(idex));
    newFig = figure(newFigs(idex));clf
    set(gcf,'Position',[1000 770   375  247])
    
    copyobj(get(oldFig,'children'), newFig);
    newFig = figure(newFig);
    set(newFig,'Name',titles{idex})
    h = gca;
    grid on
    % ch = get(gca,'Children');
    set(h,'Children',h.Children([1,3,4,2,5,6]))
    h.Children(4).Visible = 'off';
    h.Children(5).Visible = 'off';
    h.Children(6).Visible = 'off';

    h.Children(1).LineWidth = 4;
    h.Children(2).LineWidth = 4;

    % Adjust concentration scale for Nuclear GR.
    ratio = 0.5514; %Fit in "ProcessEric_Data_Feb6".
    for ich=1:3
        h.Children(ich).YData = h.Children(ich).YData/ratio;
    end
    h.Children(1).YNegativeDelta = h.Children(1).YNegativeDelta/ratio;
    h.Children(1).YPositiveDelta = h.Children(1).YPositiveDelta/ratio;

    legend('','','','Model \mu \pm \sigma','Model \mu','Data \mu \pm \sigma ')
    set(gca,'fontsize',16,'ylim',[0,50])
    xlabel('')
    title('')
    ylabel('')
end

%% To generate GR Cyt plots 
origFigs = [3,7,11];
titles = {'Cyt GR 1nM','Cyt GR 10nM','Cyt GR 100nM'};
newFigs = [1011,1012,1013];

for idex = 1:3
    oldFig = figure(origFigs(idex));
    newFig = figure(newFigs(idex));clf
    set(gcf,'Position',[1000 770   375  247])
    
    copyobj(get(oldFig,'children'), newFig);
    newFig = figure(newFig);
    set(newFig,'Name',titles{idex})
    h = gca;
    grid on
    % ch = get(gca,'Children');
    set(h,'Children',h.Children([2,5,6,1,3,4]))
    h.Children(4).Visible = 'off';
    h.Children(5).Visible = 'off';
    h.Children(6).Visible = 'off';

    h.Children(1).LineWidth = 4;
    h.Children(2).LineWidth = 4;

    legend('','','','Model \mu \pm \sigma','Model \mu','Data \mu \pm \sigma ')
    set(gca,'fontsize',16,'ylim',[0,25])
    xlabel('')
    title('')
    ylabel('')

end


%% Nuclear Distributions
origFigs = [2,6,10];
titles = {'Nuc GR 1nM','Nuc GR 10nM','Nuc GR 100nM'};
newFigs = [1021,1022,1023];

for idex = 1:3
    oldFig = figure(origFigs(idex));
    newFig = figure(newFigs(idex));clf
    set(newFig,'Position',[1000        1039        1361         199])

    set(newFig,'Name',titles{idex})
    % Get handles of all subplots in the original figure
    subplotHandles = findobj(oldFig, 'type', 'axes');

    % Iterate over each subplot and copy its contents to the new figure
    orderset = [7,6,5,4,3,1];
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
        set(h,'Children',h.Children([1:length(h.Children)]))

        grid on

        if splitReps
            for ih = [1,2,5,6,9,10]+2
                h.Children(ih).Visible = 'off';
            end
            h.Children(1)
            for ih = [1,2,5,6,9,10]
                h.Children(ih).Visible = 'on';
                h.Children(ih).LineWidth = 2;
            end
            h.Children(1).Color = 'k';
            h.Children(1).LineWidth = 4;
        else
            h.Children(1).LineWidth = 4;
            h.Children(2).LineWidth = 4;

            h.Children(3).Visible = 'off';
            h.Children(4).Visible = 'off';
        end

        % Adjust concentration scale for Nuclear GR.
        ratio = 0.5514; %Fit in "ProcessEric_Data_Feb6".
        for ich=1:length(h.Children)
            h.Children(ich).XData = h.Children(ich).XData/ratio;
        end

        set(gca,'xlim',[0,70],'ylim',[0,0.3],'FontSize',16)
        % h.Children(6).Visible = 'off'

        %text(10,0.28,[subplotHandle.Title.String(5:end),' min'],'Fontsize',16)
        ytxt = 0.25; xtxt = 20;
        if idex==3||i2==1
            N = combinedGRModel.SSITModels{3}.dataSet.nCells(i2);          
        else
            N = combinedGRModel.SSITModels{idex}.dataSet.nCells(i2-1);
        end
        t = text(xtxt,ytxt,['N = ',num2str(N)]); t.FontSize = 15;

    end

end

%% Cytoplasmic Distributions
origFigs = [2,6,10];
titles = {'Cyt GR 1nM','Cyt GR 10nM','Cyt GR 100nM'};
newFigs = [1031,1032,1033];

for idex = 1:3
    oldFig = figure(origFigs(idex));
    newFig = figure(newFigs(idex));clf
    set(newFig,'Position',[1000        1039        1361         199])

    set(newFig,'Name',titles{idex})
    % Get handles of all subplots in the original figure
    subplotHandles = findobj(oldFig, 'type', 'axes');

    % Iterate over each subplot and copy its contents to the new figure
    orderset = [7,6,5,4,3,1];
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
        set(h,'Children',h.Children([1:length(h.Children)]))

        grid on

        if splitReps
            for ih = [1,2,5,6,9,10]
                h.Children(ih).Visible = 'off';
            end
            h.Children(1)
            for ih = [1,2,5,6,9,10]+2
                h.Children(ih).Visible = 'on';
                h.Children(ih).LineWidth = 2;
            end
            h.Children(3).Color = 'k';
            h.Children(3).LineWidth = 4;
        else
            h.Children(3).LineWidth = 4;
            h.Children(4).LineWidth = 4;

            h.Children(1).Visible = 'off';
            h.Children(2).Visible = 'off';
        end

        set(gca,'xlim',[0,30],'ylim',[0,0.4],'FontSize',16)
        % h.Children(6).Visible = 'off'

        %text(6,0.42,[subplotHandle.Title.String(5:end),' min'],'Fontsize',16)
        ytxt = 0.35; xtxt = 5;
        if idex==3||i2==1
            N = combinedGRModel.SSITModels{3}.dataSet.nCells(i2);          
        else
            N = combinedGRModel.SSITModels{idex}.dataSet.nCells(i2-1);
        end
        t = text(xtxt,ytxt,['N = ',num2str(N)]); t.FontSize = 15;

    end
end

%% To generate Dusp1 means and variance plots 
origFigs = [201,301,302,303];
titles = {'Dusp1 100nM','Dusp1 GR 10nM','Nuc GR 1nM','Nuc GR 0.3nM'};
newFigs = [1041,1042,1043,1044];

for idex = 1:4
    oldFig = figure(origFigs(idex));
    newFig = figure(newFigs(idex));clf
    set(gcf,'Position',[1000 991   496  247])
    
    copyobj(get(oldFig,'children'), newFig);
    newFig = figure(newFig);
    set(newFig,'Name',titles{idex})
    h = gca;
    grid on
    legend('Model \mu \pm \sigma','Model \mu','Data \mu \pm \sigma ')
    set(gca,'fontsize',16,'ylim',[0,100])
    xlabel('')
    title('')
    ylabel('')
end

%% Dusp1 Distributions
origFigs = [221,321,322,323];
titles = {'Dusp1 100nM','Dusp1 10nM','Dusp1 1nM','Dusp1 0.3nM'};
newFigs = [1051,1052,1053,1054];

oldFig100 = figure(221);
subplotHandles100 = findobj(oldFig100, 'type', 'axes');

for idex = 1:4
    oldFig = figure(origFigs(idex));
    newFig = figure(newFigs(idex));clf
    set(newFig,'Position',[50        139        1361         199])

    set(newFig,'Name',titles{idex})
    % Get handles of all subplots in the original figure
    subplotHandles = findobj(oldFig, 'type', 'axes');

    % Iterate over each subplot and copy its contents to the new figure
    if idex==1        
        orderset = [12,10,9,8,6,4,1];
    else
        orderset = [NaN,6,5,4,3,2,1];
    end
    for i2 = 1:length(orderset)
        i = orderset(i2);
        if i2==1
            subplotHandle = subplotHandles100(12);
        else
            subplotHandle = subplotHandles(i);
        end

        % Create new subplot in the new figure
        ax = subplot(1, 7, i2, 'Parent', newFig);

        % Copy contents of the original subplot to the new subplot
        copyobj(allchild(subplotHandle), ax);

        if i2~=1
            ylabel('')
            set(gca,'yticklabels',[])
        end

        h = gca;
        set(h,'Children',h.Children([1,2]))

        grid on
        ax.XMinorGrid = 'on';

        h.Children(1).LineWidth = 4;
        h.Children(2).LineWidth = 4;

        set(gca,'xlim',[0,200],'ylim',[0,0.2],'FontSize',16)
        % h.Children(6).Visible = 'off'

        %text(60,0.17,[subplotHandle.Title.String(5:end),' min'],'Fontsize',16)

    end
end

%% Make pot of experiment design.
f2 = figure;
set(f2,'Position',[939   687   404   238])
colormap('sky')
pCol = [ModelGRfit{3}.dataSet.nCells(1),ModelGRfit{1}.dataSet.nCells([1:6])';
    ModelGRfit{3}.dataSet.nCells(1),ModelGRfit{2}.dataSet.nCells([1:6])';
    ModelGRfit{3}.dataSet.nCells([1:7])'];
% sz = size(finalExperimentDesignBFBD);
pCol(end+1,end+1) = 0;
pCol(end,end) = 450;
pcolor([-0:1:7],[0:3],pCol)
cb = colorbar;
cb.Label.FontSize = 16;
cb.Label.String = 'Number of Cells';
cb.Label.Interpreter = 'latex';
set(gca,'FontSize',16,'TickLabelInterpreter','latex',...
    'ytick',[0.5:1:3.5],'YTickLabel',[1,10,100],...
    'xtick',[0.5:1:7.5],'XTickLabel',[0,10,30,50,75,120,180]);
xlabel('Measurement Time (min)','Interpreter','latex')
ylabel('Input ($\mu$M)','Interpreter','latex')


%% Reformatting for the titration plot
set(gca,'xtick',10.^[-3:4],'ylim',[-5,150])