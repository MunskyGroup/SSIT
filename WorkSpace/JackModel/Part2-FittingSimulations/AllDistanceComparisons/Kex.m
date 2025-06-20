clear
close all 
addpath(genpath('../../../../src'));
addpath(genpath('../../Model'))

%%
num_edges = 100; % to divide and compare pdfs with this number of edges for spatial distributions 
num_param_values = 30; % number of different parameter values per parameter i.e. linspace(0, 1, num_param_values)
NCells = 2000; % Num cells per simulation
num_sims_per_parameter = 10;

kex_start = 750*10^-3;
kex_end = 750*10^2;
sim_kex = true;
calc_kex = true;

kon_default = 1e-3;
koff_default = 7.5e-5;
w_default = 0.0025;
kex_default = 750;
kr_default = 7.5e4;
D_default = [0.01,5,4]; % bound, full, part
gam_default =[0.035;0.0025;0.001];
makePlot = false;


%% Simulate: Theta = kex
kex = logspace(log10(kex_start), log10(kex_end), num_param_values);

if sim_kex
    disp('Simulating kex \n')
    results = cell([length(kex), num_sims_per_parameter]);
    for i = 1:num_param_values
        for n = 1:num_sims_per_parameter
           fprintf('Simulating n = %i, kex = %.2e \n', n, kex(i))
            results{i,n} = Simulate(kon_default, koff_default, w_default, kex_default, kex(i), D_default, gam_default, NCells, makePlot);
        end
    end
    save("kex", "results")
end
if calc_kex
    fprintf('Calculations for kex \n')
results = load("kex.mat", "results").results;
[distances, mean_distances, var_distances] = main(results, 'kex', kex, ceil(num_param_values/2), false, num_edges);
end


