function runBackgroundFit(data_file,results_file,Merged)
arguments
    data_file
    results_file
    Merged = false;
end
% this script demonstrates how to create a background Matlab job from
% within Matlab.

command1 = '';

%  ML = '/Applications/MATLAB_R2020b.app/bin/./matlab ';
tmp = which('matlab');
tmpb = find(tmp=='/'|tmp=='\');
ML = [tmp(1:tmpb(3)),'bin/./matlab '];
ARGS = '-nodesktop ';


% For BASH, we need to replace spaces and '(' or ')' in directory names
% with '\ ', '\(' or ')'.
results_file_bash = results_file;
ss = {' ','(',')'};
for j=1:3
    k = strfind(results_file_bash,ss{j});
    for i = length(k):-1:1
        results_file_bash = [results_file_bash(1:k(i)-1),'\',results_file_bash(k(i):end)];
    end
%     k = strfind(data_file,ss{j});
%     for i = length(k):-1:1
%         data_file = [data_file(1:k(i)-1),'\',data_file(k(i):end)];
%     end
end

k = strfind(results_file_bash,'.mat');
if isempty(k) 
    k=length(results_file_bash); 
else
    k = k(end);
end

LOG = ['-logfile ',results_file_bash(1:k-1),'_log.txt'];

if Merged
    SCR = [' -batch "addpath(genpath(''../src''));  ssit.parest.bgFitModel2DataMerged([],[],''',...
        data_file,''',''',results_file,''',',num2str(Merged),')" &'];
else
    SCR = [' -batch "addpath(genpath(''../src''));  ssit.parest.bgFitModel2Data([],''',...
        data_file,''',''',results_file,''',',num2str(Merged),')" &'];
end

command2 = [ML,ARGS,LOG,SCR];

[command1, command2]

system([command1, command2])
