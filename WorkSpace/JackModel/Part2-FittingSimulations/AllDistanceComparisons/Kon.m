clear
close all 
addpath(genpath('../../../../src'));
addpath(genpath('../../Model'))

%%
num_edges = 100; % to divide and compare pdfs with this number of edges for spatial distributions 
num_param_values = 30; % number of different parameter values per parameter i.e. linspace(0, 1, num_param_values)
NCells = 2000; % Num cells per simulation
num_sims_per_parameter = 10;


kon_start = -7;
kon_end = -2;
sim_kon = true;
calc_kon = true;

kon_default = 1e-3;
koff_default = 7.5e-5;
w_default = 0.0025;
kex_default = 750;
kr_default = 7.5e4;
D_default = [0.01,5,4]; % bound, full, part
gam_default =[0.035;0.0025;0.001];
makePlot = false;


%% Simulate: Theta = kon
kon = logspace(kon_start, kon_end, num_param_values);

if sim_kon
    fprintf('Simulating on \n')
    results = cell([num_param_values, num_sims_per_parameter]);
    for i = 1:num_param_values
        for n = 1:num_sims_per_parameter
            fprintf('Simulating n = %i, kon = %.2e \n', n, kon(i))
            results{i,n} = Simulate(kon(i), koff_default, w_default, kex_default, kr_default, D_default, gam_default, NCells, makePlot);
        end
    end
    save("Kon", "results")
end
if calc_kon
    fprintf('Calculations for kon \n')
    results = load("Kon.mat", "results").results;
    [distances, mean_distances, var_distances] = main(results, 'kon', kon, ceil(num_param_values/2), false, num_edges);
end


