clear
close all 
addpath(genpath('../../../../src'));
addpath(genpath('../../Model'))

%%
num_edges = 100; % to divide and compare pdfs with this number of edges for spatial distributions 
num_param_values = 30; % number of different parameter values per parameter i.e. linspace(0, 1, num_param_values)
NCells = 2000; % Num cells per simulation
num_sims_per_parameter = 10;


w_start = 0.0025*10^-3;
w_end = 0.0025*10^1;
sim_w = true;
calc_w = true;

kon_default = 1e-3;
koff_default = 7.5e-5;
w_default = 0.0025;
kex_default = 750;
kr_default = 7.5e4;
D_default = [0.01,5,4]; % bound, full, part
gam_default =[0.035;0.0025;0.001];
makePlot = false;


%% Simulate: Theta = w
w = logspace(log10(w_start), log10(w_end), num_param_values);

if sim_w
    disp('Simulating w \n')
    results = cell([length(w), num_sims_per_parameter]);
    for i = 1:num_param_values
        for n = 1:num_sims_per_parameter
           fprintf('Simulating n = %i, w = %.2e \n', n, w(i))
            results{i,n} = Simulate(kon_default, koff_default, w(i), kex_default, kr_default, D_default, gam_default, NCells, makePlot);
        end
    end
    save("w", "results")
end
if calc_w
    fprintf('Calculations for w \n')
results = load("w.mat", "results").results;
[distances, mean_distances, var_distances] = main(results, 'w', w, ceil(num_param_values/2), false, num_edges);
end


