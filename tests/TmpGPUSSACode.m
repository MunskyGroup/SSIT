function [X]=TmpGPUSSACode(x0,N_run,parametersIn,useGPU)
%This is an automatically generated MATLAB SSA Program.
%The tools used to generate this file were covered at the.
%2016 q-bio Summer School at the Colorado State University.
arguments
   x0
   N_run
   parametersIn
   useGPU=CPU
end

Nspec = 1; % Number of species.
Nt = 21; % Number of time points.
X = zeros(Nspec,Nt,N_run); % Initialize matrix of results. 
if strcmp(useGPU,'GPU')
   g = parallel.gpu.RandStream('Philox4x32-10','Seed',0); % Set seed for RNG on GPU.
   parallel.gpu.RandStream.setGlobalStream(g); % Apply RNG seed to GPU.
   x1_0_GPU = x0(1)*gpuArray.ones(1,N_run); % Specific Initial Conditions.
   TmpGPUSSACode_SSA_GPU = @(x1_0_GPU)TmpGPUSSACode_SSA(x1_0_GPU,parametersIn);
   [x1_1,x1_2,x1_3,x1_4,x1_5,x1_6,x1_7,x1_8,x1_9,x1_10,x1_11,x1_12,x1_13,x1_14,x1_15,x1_16,x1_17,x1_18,x1_19,x1_20,x1_21] = arrayfun(@TmpGPUSSACode_SSA_GPU,x1_0_GPU);
   X(1,1,:) =gather(x1_1);
   X(1,2,:) =gather(x1_2);
   X(1,3,:) =gather(x1_3);
   X(1,4,:) =gather(x1_4);
   X(1,5,:) =gather(x1_5);
   X(1,6,:) =gather(x1_6);
   X(1,7,:) =gather(x1_7);
   X(1,8,:) =gather(x1_8);
   X(1,9,:) =gather(x1_9);
   X(1,10,:) =gather(x1_10);
   X(1,11,:) =gather(x1_11);
   X(1,12,:) =gather(x1_12);
   X(1,13,:) =gather(x1_13);
   X(1,14,:) =gather(x1_14);
   X(1,15,:) =gather(x1_15);
   X(1,16,:) =gather(x1_16);
   X(1,17,:) =gather(x1_17);
   X(1,18,:) =gather(x1_18);
   X(1,19,:) =gather(x1_19);
   X(1,20,:) =gather(x1_20);
   X(1,21,:) =gather(x1_21);
elseif strcmp(useGPU,'Parallel')
   x1_0 = x0(1); % Specific Initial Conditions.
  parfor i = 1:N_run
    [x] = collectFun(x1_0,parametersIn);
    X(:,:,i) = reshape(x,[Nt,Nspec])';
  end
elseif strcmp(useGPU,'Series')
  for i = 1:N_run
   x1_0 = x0(1); % Specific Initial Conditions.
    [x1_1,x1_2,x1_3,x1_4,x1_5,x1_6,x1_7,x1_8,x1_9,x1_10,x1_11,x1_12,x1_13,x1_14,x1_15,x1_16,x1_17,x1_18,x1_19,x1_20,x1_21] = TmpGPUSSACode_SSA(x1_0,parametersIn);
    X(:,:,i) = reshape(  [x1_1,x1_2,x1_3,x1_4,x1_5,x1_6,x1_7,x1_8,x1_9,x1_10,x1_11,x1_12,x1_13,x1_14,x1_15,x1_16,x1_17,x1_18,x1_19,x1_20,x1_21],[Nt,Nspec])';
  end
end
end


function [x1_1,x1_2,x1_3,x1_4,x1_5,x1_6,x1_7,x1_8,x1_9,x1_10,x1_11,x1_12,x1_13,x1_14,x1_15,x1_16,x1_17,x1_18,x1_19,x1_20,x1_21] = TmpGPUSSACode_SSA(x1,parametersIn)
%First we define the parameters.
k1=parametersIn(1);
k2=parametersIn(2);

%Initialize the time.
t=0;
%Initialize the state.
x1new=x1;
%Start the SSA.
tstop = 0;   %Next time to print results.
while t<tstop   %Next time to print results.
  %Update the state to reflect the last reaction.
  x1=x1new;
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  r2w0=rand*w0;
  if r2w0<w1
    x1new=x1+(1);
  elseif r2w0<w1+w2
    x1new=x1+(-1);
  end
end
x1_1=x1;

tstop = 0.1;   %Next time to print results.
while t<tstop   %Next time to print results.
  %Update the state to reflect the last reaction.
  x1=x1new;
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  r2w0=rand*w0;
  if r2w0<w1
    x1new=x1+(1);
  elseif r2w0<w1+w2
    x1new=x1+(-1);
  end
end
x1_2=x1;

tstop = 0.2;   %Next time to print results.
while t<tstop   %Next time to print results.
  %Update the state to reflect the last reaction.
  x1=x1new;
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  r2w0=rand*w0;
  if r2w0<w1
    x1new=x1+(1);
  elseif r2w0<w1+w2
    x1new=x1+(-1);
  end
end
x1_3=x1;

tstop = 0.3;   %Next time to print results.
while t<tstop   %Next time to print results.
  %Update the state to reflect the last reaction.
  x1=x1new;
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  r2w0=rand*w0;
  if r2w0<w1
    x1new=x1+(1);
  elseif r2w0<w1+w2
    x1new=x1+(-1);
  end
end
x1_4=x1;

tstop = 0.4;   %Next time to print results.
while t<tstop   %Next time to print results.
  %Update the state to reflect the last reaction.
  x1=x1new;
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  r2w0=rand*w0;
  if r2w0<w1
    x1new=x1+(1);
  elseif r2w0<w1+w2
    x1new=x1+(-1);
  end
end
x1_5=x1;

tstop = 0.5;   %Next time to print results.
while t<tstop   %Next time to print results.
  %Update the state to reflect the last reaction.
  x1=x1new;
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  r2w0=rand*w0;
  if r2w0<w1
    x1new=x1+(1);
  elseif r2w0<w1+w2
    x1new=x1+(-1);
  end
end
x1_6=x1;

tstop = 0.6;   %Next time to print results.
while t<tstop   %Next time to print results.
  %Update the state to reflect the last reaction.
  x1=x1new;
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  r2w0=rand*w0;
  if r2w0<w1
    x1new=x1+(1);
  elseif r2w0<w1+w2
    x1new=x1+(-1);
  end
end
x1_7=x1;

tstop = 0.7;   %Next time to print results.
while t<tstop   %Next time to print results.
  %Update the state to reflect the last reaction.
  x1=x1new;
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  r2w0=rand*w0;
  if r2w0<w1
    x1new=x1+(1);
  elseif r2w0<w1+w2
    x1new=x1+(-1);
  end
end
x1_8=x1;

tstop = 0.8;   %Next time to print results.
while t<tstop   %Next time to print results.
  %Update the state to reflect the last reaction.
  x1=x1new;
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  r2w0=rand*w0;
  if r2w0<w1
    x1new=x1+(1);
  elseif r2w0<w1+w2
    x1new=x1+(-1);
  end
end
x1_9=x1;

tstop = 0.9;   %Next time to print results.
while t<tstop   %Next time to print results.
  %Update the state to reflect the last reaction.
  x1=x1new;
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  r2w0=rand*w0;
  if r2w0<w1
    x1new=x1+(1);
  elseif r2w0<w1+w2
    x1new=x1+(-1);
  end
end
x1_10=x1;

tstop = 1;   %Next time to print results.
while t<tstop   %Next time to print results.
  %Update the state to reflect the last reaction.
  x1=x1new;
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  r2w0=rand*w0;
  if r2w0<w1
    x1new=x1+(1);
  elseif r2w0<w1+w2
    x1new=x1+(-1);
  end
end
x1_11=x1;

tstop = 1.1;   %Next time to print results.
while t<tstop   %Next time to print results.
  %Update the state to reflect the last reaction.
  x1=x1new;
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  r2w0=rand*w0;
  if r2w0<w1
    x1new=x1+(1);
  elseif r2w0<w1+w2
    x1new=x1+(-1);
  end
end
x1_12=x1;

tstop = 1.2;   %Next time to print results.
while t<tstop   %Next time to print results.
  %Update the state to reflect the last reaction.
  x1=x1new;
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  r2w0=rand*w0;
  if r2w0<w1
    x1new=x1+(1);
  elseif r2w0<w1+w2
    x1new=x1+(-1);
  end
end
x1_13=x1;

tstop = 1.3;   %Next time to print results.
while t<tstop   %Next time to print results.
  %Update the state to reflect the last reaction.
  x1=x1new;
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  r2w0=rand*w0;
  if r2w0<w1
    x1new=x1+(1);
  elseif r2w0<w1+w2
    x1new=x1+(-1);
  end
end
x1_14=x1;

tstop = 1.4;   %Next time to print results.
while t<tstop   %Next time to print results.
  %Update the state to reflect the last reaction.
  x1=x1new;
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  r2w0=rand*w0;
  if r2w0<w1
    x1new=x1+(1);
  elseif r2w0<w1+w2
    x1new=x1+(-1);
  end
end
x1_15=x1;

tstop = 1.5;   %Next time to print results.
while t<tstop   %Next time to print results.
  %Update the state to reflect the last reaction.
  x1=x1new;
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  r2w0=rand*w0;
  if r2w0<w1
    x1new=x1+(1);
  elseif r2w0<w1+w2
    x1new=x1+(-1);
  end
end
x1_16=x1;

tstop = 1.6;   %Next time to print results.
while t<tstop   %Next time to print results.
  %Update the state to reflect the last reaction.
  x1=x1new;
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  r2w0=rand*w0;
  if r2w0<w1
    x1new=x1+(1);
  elseif r2w0<w1+w2
    x1new=x1+(-1);
  end
end
x1_17=x1;

tstop = 1.7;   %Next time to print results.
while t<tstop   %Next time to print results.
  %Update the state to reflect the last reaction.
  x1=x1new;
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  r2w0=rand*w0;
  if r2w0<w1
    x1new=x1+(1);
  elseif r2w0<w1+w2
    x1new=x1+(-1);
  end
end
x1_18=x1;

tstop = 1.8;   %Next time to print results.
while t<tstop   %Next time to print results.
  %Update the state to reflect the last reaction.
  x1=x1new;
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  r2w0=rand*w0;
  if r2w0<w1
    x1new=x1+(1);
  elseif r2w0<w1+w2
    x1new=x1+(-1);
  end
end
x1_19=x1;

tstop = 1.9;   %Next time to print results.
while t<tstop   %Next time to print results.
  %Update the state to reflect the last reaction.
  x1=x1new;
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  r2w0=rand*w0;
  if r2w0<w1
    x1new=x1+(1);
  elseif r2w0<w1+w2
    x1new=x1+(-1);
  end
end
x1_20=x1;

tstop = 2;   %Next time to print results.
while t<tstop   %Next time to print results.
  %Update the state to reflect the last reaction.
  x1=x1new;
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  r2w0=rand*w0;
  if r2w0<w1
    x1new=x1+(1);
  elseif r2w0<w1+w2
    x1new=x1+(-1);
  end
end
x1_21=x1;

end
function [x] = collectFun(x1,parametersIn)
% This function runs the SSA and gathers the results into a single matrix.
[x1_1,x1_2,x1_3,x1_4,x1_5,x1_6,x1_7,x1_8,x1_9,x1_10,x1_11,x1_12,x1_13,x1_14,x1_15,x1_16,x1_17,x1_18,x1_19,x1_20,x1_21] = TmpGPUSSACode_SSA(x1,parametersIn);
x=[x1_1,x1_2,x1_3,x1_4,x1_5,x1_6,x1_7,x1_8,x1_9,x1_10,x1_11,x1_12,x1_13,x1_14,x1_15,x1_16,x1_17,x1_18,x1_19,x1_20,x1_21];
end
