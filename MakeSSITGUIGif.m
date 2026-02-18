% makeGUIgif
install;  
stubDir = tempname;
mkdir(stubDir);

% --- Write stub questdlg to that folder ---
fid = fopen(fullfile(stubDir, 'questdlg.m'), 'w');
fprintf(fid, 'function sel = questdlg(varargin)\n');
fprintf(fid, 'sel = ''Yes - Overwrite'';\n');  % always return "Yes"
fprintf(fid, 'end\n');
fclose(fid);

% --- Put stub in front of MATLAB path ---
origPath = addpath(stubDir);
cleanup = onCleanup(@() path(origPath)); % restore later

close all force
A = SSITGUI;

gifFile = 'SSITAnimation.gif';

% --- capture current figure ---
frame = getframe(A.UIFigure);
im = frame2im(frame);
[Amovie,map] = rgb2ind(im,256);

% --- first frame: create GIF ---
imwrite(Amovie,map,gifFile,'gif', ...
    'LoopCount',inf, ...
    'DelayTime',0.2);

A.TabGroup.SelectedTab = A.ModelLoadingandBuildingTab;
A = loadModelBP(A, [], 'tests/test_data/GRDusp1ModelTestLibrary.mat');
A.ChooseSSITModel.Value = 'ModelDUSP1_100nM';
A = ChooseSSITModelValue(A);

updateGif(A,gifFile);

A.TabGroup.SelectedTab = A.StochasticSimulationTab;
A.PrintTimesEditField.Value = '[0:10:180]';
A.PrintTimesEditField.ValueChangedFcn(A,[]);

updateGif(A,gifFile);

A.SsaRunButton.ButtonPushedFcn(A, []);
updateGif(A,gifFile);

A.SpeciestoShowListBox.Value = 'rna';
A.SSAYmaxEditField.Value = 100;
A.SSAXmaxEditField.Value = 180;
A.SSA_hist_XmaxEditField.Value = 100;
A.SSA_hist_YmaxEditField.Value = 0.4;
A.SSA_hist_XmaxEditField.ValueChangedFcn(A,[]);

A.SpeciestoShowListBox.ValueChangedFcn(A,[]);
updateGif(A,gifFile);
%%
for x = 0:10:180
    A.SsaTimeSlider.Value = x;
    A.SsaTimeSlider.ValueChangedFcn(A,[]);
    updateGif(A,gifFile);
end


function gifFile = updateGif(A,gifFile)
% --- capture again and append ---
frame = getframe(A.UIFigure);
im = frame2im(frame);
[Amovie,map] = rgb2ind(im,256);

imwrite(Amovie,map,gifFile,'gif', ...
    'WriteMode','append', ...
    'DelayTime',0.2);
end


