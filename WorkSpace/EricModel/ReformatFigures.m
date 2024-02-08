%% To generate GR Nuc plots 
origFigs = [3,7,11];
titles = {'Nuc GR 1nM','Nuc GR 10nM','Nuc GR 100nM'};
newFigs = [1001,1002,1003];

for idex = 1:3
    oldFig = figure(origFigs(idex));
    newFig = figure(newFigs(idex));clf
    set(gcf,'Position',[1000 991   496  247])
    
    copyobj(get(oldFig,'children'), newFig);
    newFig = figure(newFig);
    set(newFig,'Name',titles{idex})
    h = gca;
    % ch = get(gca,'Children');
    set(h,'Children',h.Children([1,3,4,2,5,6]))
    h.Children(4).Visible = 'off'
    h.Children(5).Visible = 'off'
    h.Children(6).Visible = 'off'
    legend('','','','Model \mu \pm \sigma','Model \mu','Data \mu \pm \sigma ')
    set(gca,'fontsize',16)
    xlabel('')
    title('')
end

%% To generate GR Cyt plots 
origFigs = [3,7,11];
titles = {'Cyt GR 1nM','Cyt GR 10nM','Cyt GR 100nM'};
newFigs = [1011,1012,1013];

for idex = 1:3
    oldFig = figure(origFigs(idex));
    newFig = figure(newFigs(idex));clf
    set(gcf,'Position',[1000 991   496  247])
    
    copyobj(get(oldFig,'children'), newFig);
    newFig = figure(newFig);
    set(newFig,'Name',titles{idex})
    h = gca;
    % ch = get(gca,'Children');
    set(h,'Children',h.Children([2,5,6,1,3,4]))
    h.Children(4).Visible = 'off'
    h.Children(5).Visible = 'off'
    h.Children(6).Visible = 'off'
    legend('','','','Model \mu \pm \sigma','Model \mu','Data \mu \pm \sigma ')
    set(gca,'fontsize',16)
    xlabel('')
    title('')
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

        h.Children(3).Visible = 'off'
        h.Children(4).Visible = 'off'

        set(gca,'xlim',[0,20],'ylim',[0,0.3],'FontSize',16)
        % h.Children(6).Visible = 'off'
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

        h.Children(1).Visible = 'off'
        h.Children(2).Visible = 'off'

        set(gca,'xlim',[0,20],'ylim',[0,0.5],'FontSize',16)
        % h.Children(6).Visible = 'off'
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
    legend('Model \mu \pm \sigma','Model \mu','Data \mu \pm \sigma ')
    set(gca,'fontsize',16)
    xlabel('')
    title('')
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
    set(newFig,'Position',[1000        1039        1361         199])

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

        % h.Children(1).Visible = 'off'
        % h.Children(2).Visible = 'off'

        set(gca,'xlim',[0,200],'ylim',[0,0.2],'FontSize',16)
        % h.Children(6).Visible = 'off'
    end
end
