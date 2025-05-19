%% Calc p(x, t) for simulations
clear
close all 
addpath(genpath('../../../../src'));
addpath(genpath('../../Model'))

s = settings;
s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light";


%% Simulate 
makePlot = false;
NCells = 2000;
kon = 1e-3;
koff = 7.5e-5;
w = 0.0025;
kex = 750;
kr = 7.5e4;
D = [0.01,5,4];
gam =[0.035;0.0025;0.001];

results = Simulate(kon, koff, w, kex, kr, D, gam, NCells, makePlot);

RNACountBins = results.binCounts; 
TSbins = results.binTS;

%% Display the P(x,t)
for ibinloc = 1:5
    iCellsInTSBin = (TSbins == ibinloc);
    figure()
    hold on
    for iBin = 2
        h = histogram(RNACountBins(iBin,iCellsInTSBin), ...
                                        'BinMethod', 'auto', ...
                                        'Normalization', 'probability');
        pdf = h.Values;
        edges = h.BinEdges;
    end
end


%% test function

P = Prob(results);



























%% Functions
function P = Prob(results, edges)

    if nargin < 2 || isempty(edges)
        edges = [0,10,20,30,40,50,60,70,80,90,100];
    end


    RNACountBins = results.binCounts; 
    TSbins = results.binTS;

    P = cell([5, size(RNACountBins, 1)]);

    for ibinloc = 1:5
        iCellsInTSBin = (TSbins == ibinloc);
        for iBin = 1:size(RNACountBins, 1)
            data = RNACountBins(iBin, iCellsInTSBin);
            [pdf, ~] = histcounts(data, 'BinEdges', edges, 'Normalization', 'probability');
            P{ibinloc, iBin} = pdf;
        end
    end

end






















