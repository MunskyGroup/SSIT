%%
clear all
fun_name = 'TEST';    % Name of the automatic matlab function to be generated.
% k = [10,1,10,1,2,5,2,2];  % Parameters
% w = {{'k1'},{'k2*x1'},{'k3*x1'},{'k4*x2'},{'k5*x2'},{'k6*x3'},{'k7*x3'},{'k8*x4'}};
% Propensity functions.

sbmlobj = sbmlimport('00146-sbml-l1v2.xml');
[S,objSpecies,objReactions]= getstoichmatrix(sbmlobj);

for i=1:length(sbmlobj.Parameters)
    k(i) = sbmlobj.Parameters(i).Value;
end

for i=1:length(sbmlobj.Species)
    x0(i) = sbmlobj.Species(i).InitialAmount;
end
x0 = x0*10/min(x0(x0>0));

for i =1:size(S,2)
    txt = ['k',num2str(i)];
    for j=1:size(sbmlobj.Reactions(i).Reactants)
        txt = [txt,'*x',sbmlobj.Reactions(i).Reactants(j).Name(2:end)];
    end   
    w{i} = txt;
end

tstop = 5; %final time.
N_t = 10;     % Number of intermediate time points for results.
tprint = linspace(0,tstop,N_t);  % times of the actual results.
N_run = 1000; % Number of runs to perform.
N_split=4;

% Call code to create the GPU friendly SSA.
Write_GPU_SSA(k,w,S,tprint,fun_name);
%%
tic
[X]=TEST(x0,N_run,N_split);
toc

%%
figure(1); clf
for i=1:length(x0)
    for j=1:N_t
        subplot(length(x0),N_t,(i-1)*N_t+j);
        hist(squeeze(X(i,j,:)))
    end
end
