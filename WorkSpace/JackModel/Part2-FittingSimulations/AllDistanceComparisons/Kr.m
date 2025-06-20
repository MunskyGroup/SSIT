clear
close all 
addpath(genpath('../../../../src'));
addpath(genpath('../../Model'))

%%
num_edges = 100; % to divide and compare pdfs with this number of edges for spatial distributions 
num_param_values = 30; % number of different parameter values per parameter i.e. linspace(0, 1, num_param_values)
NCells = 2000; % Num cells per simulation
num_sims_per_parameter = 10;


kr_start = 7.5e2;
kr_end = 7.5e6;
sim_kr = true;
calc_kr = true;

kon_default = 1e-3;
koff_default = 7.5e-5;
w_default = 0.0025;
kex_default = 750;
kr_default = 7.5e4;
D_default = [0.01,5,4]; % bound, full, part
gam_default =[0.035;0.0025;0.001];
makePlot = false;


%% Simulate: Theta = kr
kr = logspace(log10(7.5e-6), log10(7.5e-2), num_param_values);

if sim_kr
    disp('Simulating kr \n')
    results = cell([length(kr), num_sims_per_parameter]);
    for i = 1:num_param_values
        for n = 1:num_sims_per_parameter
           fprintf('Simulating n = %i, kr = %.2e \n', n, kr(i))
            results{i,n} = Simulate(kon_default, koff_default, w_default, kex_default, kr(i), D_default, gam_default, NCells, makePlot);
        end
    end
    save("kr", "results")
end
if calc_kr
    fprintf('Calculations for kr \n')
results = load("kr.mat", "results").results;
[distances, mean_distances, var_distances] = main(results, 'kr', kr, ceil(num_param_values/2), false, num_edges);
end


