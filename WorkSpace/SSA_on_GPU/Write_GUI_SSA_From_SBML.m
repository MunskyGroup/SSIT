%%
clear all
reset(gpuDevice([]))
clear all
fun_name = 'qbSS_AutoGen_from_SBML_SSA';    % Name of the automatic matlab function to be generated.

sbmlobj = sbmlimport('SBML_test_cases/00146/00146-sbml-l1v2.xml');
[S,objSpecies,objReactions]= getstoichmatrix(sbmlobj);

for i=1:length(sbmlobj.Parameters)
    k(i) = sbmlobj.Parameters(i).Value;
end

for i=1:length(sbmlobj.Species)
    x0(i) = sbmlobj.Species(i).InitialAmount;
end
Vol = min(x0(x0>0))/10000;
x0 = x0/Vol;

for i =1:size(S,2)
    txt = ['k',num2str(i)];
    Vol_Power = length(sbmlobj.Reactions(i).Reactants)-1;
    for j=1:size(sbmlobj.Reactions(i).Reactants)
        txt = [txt,'*x',sbmlobj.Reactions(i).Reactants(j).Name(2:end)];
    end   
    w{i} = [txt,'*',num2str(Vol^(Vol_Power))];
end
tstop = 1000; %final time.
N_t = 5;     % Number of intermediate time points for results.
tprint = linspace(0,tstop,N_t);  % times of the actual results.
N_run = 1000; % Number of runs to perform.
N_split=2;

% Call code to create the GPU friendly SSA.
Write_GPU_SSA(k,w,S,tprint,fun_name);

%%
fun = str2func(fun_name);
% Convert the function name string to a function handle.
tic
[X]=fun(x0,N_run,N_split);
toc

%%
figure(1); clf
for i=1:length(x0)
    for j=1:N_t
        subplot(length(x0),N_t,(i-1)*N_t+j);
        hist(squeeze(X(i,j,:)))
    end
end
