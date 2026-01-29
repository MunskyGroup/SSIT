function exportSsaData(app,filename)
% This function serves to export the SSA Trajectory data into an Excel
% sheet. Use READCELL to read these sheets.

if nargin<2
    % Open a dialog box to ask for filename
    [filename,pathname] = uiputfile('*.xlsx', 'Save SSA data as:');
end

numSamples = size(app.SSITModel.Solutions.trajs,3);
numTimes = size(app.SSITModel.Solutions.trajs,2);
numSpecies = size(app.SSITModel.Solutions.trajs,1);
%preallocate
ssaDat = zeros(numSamples*numTimes+1,numSpecies);
%load in the data from app
Trajectories = reshape(app.SSITModel.Solutions.trajs,...
    [numSpecies,numSamples*numTimes]);
Trajectories = Trajectories';
ssaDat(2:end,2:2+numSpecies-1) = Trajectories;
% Add distorted values
% if isfield(app.)
TrajectoriesDistorted = reshape(app.SSITModel.Solutions.trajsDistorted,...
    [numSpecies,numSamples*numTimes]);
TrajectoriesDistorted = TrajectoriesDistorted';
ssaDat(2:end,2+numSpecies:2+2*numSpecies-1) = TrajectoriesDistorted;
%construct time vector
ini_time = eval(app.PrintTimesEditField.Value); %initialize while loop
i = 1;
while i < numSamples
    if i == 1
        time = horzcat(ini_time,ini_time);
        i = i+1;
    else
        time = horzcat(time,ini_time);
        i = i+1;
    end
end
%load time vector into data matrix to be exported
ssaDat(2:end,1) = time;
ssaDat = num2cell(ssaDat); %convert to cell to populate labels
specDistortLAbel = cellfun(@(s)[s '_distorted'], app.SSITModel.species', 'UniformOutput', false);
ssaDat(1,:) = ['time',app.SSITModel.species',specDistortLAbel]; % populate labels
writecell(ssaDat,[pathname,filename])