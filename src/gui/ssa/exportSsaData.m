function exportSsaData(app,filename)
% This function serves to export the SSA Trajectory data into an Excel
% sheet. Use READCELL to read these sheets.

if nargin<2
    % Open a dialog box to ask for filename
    [filename,pathname] = uiputfile('*.xlsx', 'Save SSA data as:');
end

numSamples = app.SsaNumSimField.Value;
numTimes = length(eval(app.PrintTimesEditField.Value));
%preallocate
ssaDat = zeros(numSamples*numTimes+1,4);
%load in the data from app
Trajectories = reshape(app.StochasticSimulationTabOutputs.samples,...
    [3,numSamples*numTimes]);
Trajectories = Trajectories';
ssaDat(2:end,2:4) = Trajectories;
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
ssaDat(1,:) = {'Time','x1','x2','x3'}; % populate labels
writecell(ssaDat,[pathname,filename])







