
%% To generate GR plots 
f1 = figure(3);
newFig = figure(105);clf
set(gcf,'Position',[1000 991   496  247])
clf(newFig)
copyobj(get(f1,'children'), newFig);

fn = figure(newFig);
h = gca;
% ch = get(gca,'Children');
set(h,'Children',h.Children([1,3,4,2,5,6]))
h.Children(4).Visible = 'off'
h.Children(5).Visible = 'off'
h.Children(6).Visible = 'off'
leg = legend('','','','Model \mu \pm \sigma','Model \mu','Data \mu \pm \sigma ')
set(gca,'fontsize',16)
xlabel('')
title('')
set(gcf,'Name','Nuc GR Vs Time')

%%
newFig = figure(106);clf
set(gcf,'Position',[1000 991   496  247])
clf(newFig)
copyobj(get(f1,'children'), newFig);

fn = figure(newFig);
get(newFig,'Name')
h = gca;
% ch = get(gca,'Children');
set(h,'Children',h.Children([2,5,6,1,3,4]))
h.Children(4).Visible = 'off'
h.Children(5).Visible = 'off'
h.Children(6).Visible = 'off'
leg = legend('','','','Model \mu \pm \sigma','Model \mu','Data \mu \pm \sigma ')
set(gca,'fontsize',16)
xlabel('')
set(gcf,'Name','Cyt GR Vs Time')

%% Cytoplasmic
f2 = figure(2);

% Create a new figure
f2new = figure(107); clf
set(f2new,'Name','Nuclear GR Distributuions')
set(f2new,'Position',[1000        1039        1361         199])
% Get handles of all subplots in the original figure
subplotHandles = findobj(f2, 'type', 'axes');

% Iterate over each subplot and copy its contents to the new figure
orderset = [7,6,5,4,3,2,1];
for i2 = 1:length(orderset)
    i = orderset(i2);
    subplotHandle = subplotHandles(i);
    
    % Create new subplot in the new figure
    ax = subplot(1, 7, i2, 'Parent', f2new);
    
    % Copy contents of the original subplot to the new subplot
    copyobj(allchild(subplotHandle), ax);

    if i2~=1
        ylabel('')
        set(gca,'yticklabels',[])
    end

    h = gca;
    set(h,'Children',h.Children([1,2,3,4]))
    
    h.Children(3).Visible = 'off'
    h.Children(4).Visible = 'off'

    set(gca,'xlim',[0,20],'ylim',[0,0.25],'FontSize',16)
    % h.Children(6).Visible = 'off'


end

%% Nuclear
f2 = figure(2);

% Create a new figure
f2new = figure(108); clf
set(f2new,'Name','Cytoplasmic GR Distributuions')
set(f2new,'Position',[1000        1039        1361         199])
% Get handles of all subplots in the original figure
subplotHandles = findobj(f2, 'type', 'axes');

% Iterate over each subplot and copy its contents to the new figure
orderset = [7,6,5,4,3,2,1];
for i2 = 1:length(orderset)
    i = orderset(i2);
    subplotHandle = subplotHandles(i);
    
    % Create new subplot in the new figure
    ax = subplot(1, 7, i2, 'Parent', f2new);
    
    % Copy contents of the original subplot to the new subplot
    copyobj(allchild(subplotHandle), ax);

    if i2~=1
        ylabel('')
        set(gca,'yticklabels',[])
    end

    h = gca;
    set(h,'Children',h.Children([1,2,3,4]))
    
    h.Children(1).Visible = 'off'
    h.Children(2).Visible = 'off'

    set(gca,'xlim',[0,20],'ylim',[0,0.5],'FontSize',16)
    % h.Children(6).Visible = 'off'


end
