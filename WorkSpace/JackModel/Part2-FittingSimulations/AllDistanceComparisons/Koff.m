clear
close all 
addpath(genpath('../../../../src'));
addpath(genpath('../../Model'))

%%
num_edges = 100; % to divide and compare pdfs with this number of edges for spatial distributions 
num_param_values = 30; % number of different parameter values per parameter i.e. linspace(0, 1, num_param_values)
NCells = 2000; % Num cells per simulation
num_sims_per_parameter = 10;

koff_start = 7.5e-6;
koff_end = 7.5e-2;
sim_koff = true;
calc_koff = true;

kon_default = 1e-3;
koff_default = 7.5e-5;
w_default = 0.0025;
kex_default = 750;
kr_default = 7.5e4;
D_default = [0.01,5,4]; % bound, full, part
gam_default =[0.035;0.0025;0.001];
makePlot = false;


%% Simulate: Theta = koff
koff = logspace(log10(koff_start), log10(koff_end), num_param_values);

if sim_koff
    disp('Simulating koff \n')
    results = cell([length(koff), num_sims_per_parameter]);
    for i = 1:num_param_values
        for n = 1:num_sims_per_parameter
           fprintf('Simulating n = %i, koff = %.2e \n', n, koff(i))
            results{i,n} = Simulate(kon_default, koff(i), w_default, kex_default, kr_default, D_default, gam_default, NCells, makePlot);
        end
    end
    save("koff", "results")
end
if calc_koff
    fprintf('Calculations for koff \n')
results = load("koff.mat", "results").results;
[distances, mean_distances, var_distances] = main(results, 'koff', koff, ceil(num_param_values/2), false, num_edges);
end

