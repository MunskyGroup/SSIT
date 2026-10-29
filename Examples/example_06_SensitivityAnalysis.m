%% SSIT/Examples/example_06_SensitivityAnalysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 3.2.6: Sensitivity analysis 
%   * Compute model sensitivities to small changes in the model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminaries
% Use the models from example_01_CreateSSITModels and  
% computed FSP solutions from example_04_SolveSSITModels_FSP

% clear
% close all

% example_01_CreateSSITModels  
% example_04_SolveSSITModels_FSP

%% Load pre-computed FSP solutions:
load('ExampleSaveFiles/example_4_SolveSSITModels_FSP.mat')

% View model summaries
Model.summarizeModel
STL1.summarizeModel
STL1_4state.summarizeModel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(1): Solve sensitivities of the bursting gene model
%  from example_01_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Solve FSP sensitivities
% Solve the sensitivity problem:
Model = Model.solve(solver='fspSens');

% Plot the results from the sensitivity analysis:
Model.plotFSP(speciesNames=Model.species(3), plotType='sens',...
   indTimes=40, lineProps={'linewidth',3}, AxisLabelSize=12,...
   TickLabelSize=12, XLim=[0,10], TitleFontSize=22,...
   Title="Bursting Gene", Colors=[0.93,0.69,0.13])
shareAxisLabels(gcf, 14)  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(2): Solve sensitivities of the time-varying STL1 yeast model
%  from example_01_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Solve FSP sensitivities
% Solve the sensitivity problem: 
STL1 = STL1.solve(solver='fspSens'); 

% Plot the results from the sensitivity analysis:
STL1.plotFSP(speciesNames=STL1.species(3), plotType='sens',...
    indTimes=40, lineProps={'linewidth',3}, AxisLabelSize=12,...
    TickLabelSize=12, XLim=[0,10], Title="STL1",...
    TitleFontSize=22, Colors=[0.93,0.69,0.13])
shareAxisLabels(gcf, 14)  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex(3): Solve sensitivities of the 4-state time-varying STL1 yeast model
%  from example_01_CreateSSITModels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve FSP sensitivities
% Solve the sensitivity problem:
STL1_4state = STL1_4state.solve(solver='fspSens');

% Plot the results from the sensitivity analysis
STL1_4state.plotFSP(speciesNames=STL1_4state.species(5),...
    plotType='sens', indTimes=51, lineProps={'linewidth',3},...
    Colors=[0.23,0.67,0.2], AxisLabelSize=15, TickLabelSize=12,...
    XLim=[0,100], Title="4-state STL1 (t=25 min.)", TitleFontSize=24)
shareAxisLabels(gcf, 17)  


%% Save models & sensitivities
saveNames = unique({ ...
    'Model'
    'STL1'
    'STL1_4state'
    });
    
save('ExampleSaveFiles/example_6_SensitivityAnalysis',saveNames{:})

%% Local helper functions
function shareAxisLabels(fig, fontSize)
% Replace the per-subplot x/y labels with a single shared label per figure.
ax = findobj(fig, 'Type', 'axes');
if isempty(ax), return; end
xl = ax(1).XLabel.String;
yl = ax(1).YLabel.String;
for k = 1:numel(ax)
    xlabel(ax(k), '');
    ylabel(ax(k), '');
end
% Invisible axes spanning all subplots, used only to hold the shared labels
% (works in all MATLAB versions; supxlabel/supylabel need R2020b+).
pos = vertcat(ax.Position);
left = min(pos(:,1));  bottom = min(pos(:,2));
right = max(pos(:,1) + pos(:,3));  top = max(pos(:,2) + pos(:,4));
bigAx = axes(fig, 'Position', [left, bottom, right-left, top-bottom], ...
    'Visible', 'off', 'HitTest', 'off', 'HandleVisibility', 'off');
bigAx.XLabel.Visible = 'on';
bigAx.YLabel.Visible = 'on';
xlabel(bigAx, xl, 'FontSize', fontSize);
ylabel(bigAx, yl, 'FontSize', fontSize);
% Nudge labels outward so they clear the subplot tick labels
bigAx.XLabel.Units = 'normalized';  bigAx.XLabel.Position(2) = -0.08;
bigAx.YLabel.Units = 'normalized';  bigAx.YLabel.Position(1) = -0.06;
end
