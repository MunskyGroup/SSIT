clear
close all 
addpath(genpath('../../../../src'));
addpath(genpath('../../Model'));

% c = parcluster('local');
% c.NumWorkers = 4;
% parpool(c,4)



nChains = 2;
nSamples = 100;

NCells = 200;
makePlot = false;

kon_true = 1e-3; % define these in log space 
koff_true = 7.5e-5;
w_true = 0.0025;
kex_true = 750;
kr_true = 7.5e4;
D_true = [0.01,5,4]; % bound, full, part
gam_true =[0.035;0.0025;0.001];

if false
    trueResults = Simulate(kon_true, koff_true, w_true, kex_true, kr_true, D_true, gam_true, NCells, makePlot);
    save('trueResults', 'trueResults')
end

load('trueResults.mat')

delta = [kon_true, koff_true, w_true, kex_true, kr_true, D_true, gam_true'] * 0.25;
Proposal = @(x) ...
    [x(1) + (rand()-0.5)*2*delta(1), ... % use randn in log *0.1 to move an order of mag
     x(2) + (rand()-0.5)*2*delta(2), ...
     x(3) + (rand()-0.5)*2*delta(3), ...
     x(4) + (rand()-0.5)*2*delta(4), ...
     x(5) + (rand()-0.5)*2*delta(5), ...
     x(6) + (rand()-0.5)*2*delta(6), ...
     x(7) + (rand()-0.5)*2*delta(7), ...
     x(8) + (rand()-0.5)*2*delta(8), ...
     x(9) + (rand()-0.5)*2*delta(9), ...
     x(10) + (rand()-0.5)*2*delta(10), ...
     x(11) + (rand()-0.5)*2*delta(11)];

x0 = [kon_true, koff_true, w_true, kex_true, kr_true, D_true, gam_true'];
distanceMetric = @(x,y) DistanceFunction(x,y, 20);
Simulator = @(x) Simulate(x(1), x(2), x(3), x(4), x(5), squeeze(x(6:8)), squeeze(x(9:11))', NCells, makePlot);

[Distances, Accepted, Parameters] = MHABC(nChains, nSamples, trueResults, Simulator, Proposal, x0, distanceMetric);

save('ABC_Results.mat', "Parameters","Accepted","Distances")


