clear
close all
addpath(genpath('../../../SSIT/src'));

%% Create Simple Gene Regulatory network
% onGene <-> offGene
% onGene -> mRNA + onGene
% mRNA -> 0
% mRNA of gene 1 can either influence the kon and koff

%% No Influence genes (ni) - Build Model
niModel = SSIT;    
niModel.species = {'G1offGene'; 'G1onGene'; 'G1mRNA'; 'G2offGene'; 'G2onGene'; ...
                   'G2mRNA'}; % TODO: Change the order of these from smallest to largest
niModel.stoichiometry = [-1,1,0,0,  0,0,0,0; ...
                         1,-1,0,0,  0,0,0,0; ...
                         0,0,1,-1,  0,0,0,0; ...
                         0,0,0,0,  -1,1,0,0; ...
                         0,0,0,0,   1,-1,0,0;...
                         0,0,0,0,   0,0,1,-1;
                        ]; 
% niModel.propensityFunctions = {'G1_kon * G1_offGene';'G1_koff * G1_onGene';...
%                                  'G1_kr * G1_onGene';'G1_gr * G1_mRNA'; ...
%                                  'G2_kon * G2_offGene';'G2_koff * G2_onGene';...
%                                  'G2_kr * G2_onGene';'G2_gr * G2_mRNA'; ...
%                                  }; 
% niModel.parameters = ({'G1_kon',30; 'G1_koff',30; 'G1_kr',100; 'G1_gr',0.005; ...
%                          'G2_kon',30; 'G2_koff',30; 'G2_kr',100; 'G2_gr',0.005;});

niModel.propensityFunctions = {'kon * G1offGene';'koff * G1onGene';...
                                 'kr * G1onGene';'gr * G1mRNA'; ...
                                 'kon * G2offGene';'koff * G2onGene';...
                                 'kr * G2onGene';'gr * G2mRNA'; ...
                                 }; 
niModel.parameters = ({'kon',30; 'koff',30; 'kr',1; 'gr',0.005;});

niModel.initialCondition = [1;0;0;1;0;0]; 
% niModel.summarizeModel
niModel.tSpan = linspace(0,1000,6);

niModel = niModel.formPropensitiesGeneral('No_influence_2Genes');  % This line is nes

%% No Influence genes (ni) - Solve FSP
niModel.fspOptions.verbose = true;
niModel.solutionScheme = 'FSP';    % Set solutions scheme to FSP.
niModel.fspOptions.fspTol = 1e-5;  % Set FSP error tolerance.
tic
[FSPsoln,niModel.fspOptions.bounds,niModel] = niModel.solve;  % Solve the FSP analysis
toc

niModel.makePlot(FSPsoln,'meansAndDevs',[],[],1,{'linewidth',3,'color',[0,1,1]}) % Make plot of mean vs. time.
return
%% No Influence genes (ni) - Solve SSA
niModel.solutionScheme = 'SSA';  % Set solution scheme to SSA.
niModel.ssaOptions.Nexp = 1;   % Number of independent data sets to generate.
niModel.ssaOptions.nSimsPerExpt = 1; % Number of cells to include at each time point for each data set.
% niModel.ssaOptions.applyPDO = true; % Include the distortion in the SSA data.
SSASoln = niModel.solve;


%% 

niModel.makePlot(SSASoln,'trajectories',[],[],4) % Make some plots.
% F1.makePlot(FSPsoln,'meansAndDevs',[],[],4,...
%     {'linewidth',4,'color',[0,1,1],'Marker','s','MarkerSize',20}) % Add FSP Solution to plot.