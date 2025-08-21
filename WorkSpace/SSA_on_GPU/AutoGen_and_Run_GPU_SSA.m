clear all
fun_name = 'qbSS_AutoGen_SSA_Brian';    
% Name of the automatic matlab function to be generated.

k = [10,1,10,1,2,5,2,2];  
% Parameters for the model.

w = {{'k1'},{'k2*x1'},{'k3*x1'},{'k4*x2'},{'k5*x2'},{'k6*x3'},{'k7*x3'},{'k8*x4'}};
% Propensity functions for the model.

S = [1 -1 0  0 0 0 0 0;...   % Stoichiometry matrix.
    0  0 1 -1 -1 0 0 0;...
    0  0  0  0 1 -1 -1 0;...
    0  0  0  0 0 0  1 -1]; 

x0 = [5;0;0;0]; 
% initial condition.

tstop = 500; 
%final time.

N_t = 3;     
% Number of intermediate time points at which to print results.

tprint = linspace(0,tstop,N_t);  
% times of the actual results.

% Call code to write a GPU friendly SSA code.
Write_GPU_SSA(k,w,S,tprint,fun_name);
return
%%

N_run = 1000000; 
% Number of runs to perform.

N_split=gpuDeviceCount;
% Number of graphics cards over which to split the SSA runs.

fun = str2func(fun_name);
% Convert the function name string to a function handle.

tic
[X]=fun(x0,N_run,N_split);
% Call the GPU SSA code previously generated.
toc

%% Plot some results.
figure(1); clf
for i=1:length(x0)
    for j=1:N_t        
        subplot(length(x0),N_t,(i-1)*N_t+j);        
        hist(squeeze(X(i,j,:)))
        set(gca,'fontsize',12);
        if i==1
            title(['Time = ',num2str(tprint(j))]);
        end
        xlabel(['Pop. of Species ',num2str(i)]);
        if j==1
            ylabel(['Probability']);
        end
    end
end