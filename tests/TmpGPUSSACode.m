function [X]=TmpGPUSSACode(x0,N_run,useGPU)
%This is an automatically generated MATLAB SSA Program.
%The tools used to generate this file were covered at the.
%2016 q-bio Summer School at the Colorado State University.
arguments
   x0
   N_run
   useGPU=CPU
end

Nspec = 1; % Number of species.
Nt = 11; % Number of time points.
X = zeros(Nspec,Nt,N_run); % Initialize matrix of results. 
if strcmp(useGPU,'GPU')
   g = parallel.gpu.RandStream('Philox4x32-10','Seed',0); % Set seed for RNG on GPU.
   parallel.gpu.RandStream.setGlobalStream(g); % Apply RNG seed to GPU.
   x1_0_GPU = x0(1)*gpuArray.ones(1,N_run); % Specific Initial Conditions.
   [x1_1,x1_2,x1_3,x1_4,x1_5,x1_6,x1_7,x1_8,x1_9,x1_10,x1_11] = arrayfun(@TmpGPUSSACode_SSA,x1_0_GPU);
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
elseif strcmp(useGPU,'Parallel')
   x1_0 = x0(1); % Specific Initial Conditions.
  parfor i = 1:N_run
    [x1_1,x1_2,x1_3,x1_4,x1_5,x1_6,x1_7,x1_8,x1_9,x1_10,x1_11] = TmpGPUSSACode_SSA(x1_0);
    X(:,:,i) = reshape(  [x1_1,x1_2,x1_3,x1_4,x1_5,x1_6,x1_7,x1_8,x1_9,x1_10,x1_11],[Nspec,Nt]);
  end
elseif strcmp(useGPU,'Series')
  for i = 1:N_run
   x1_0 = x0(1); % Specific Initial Conditions.
    [x1_1,x1_2,x1_3,x1_4,x1_5,x1_6,x1_7,x1_8,x1_9,x1_10,x1_11] = TmpGPUSSACode_SSA(x1_0);
    X(:,:,i) = reshape(  [x1_1,x1_2,x1_3,x1_4,x1_5,x1_6,x1_7,x1_8,x1_9,x1_10,x1_11],[Nspec,Nt]);
  end
end


function [x1_1,x1_2,x1_3,x1_4,x1_5,x1_6,x1_7,x1_8,x1_9,x1_10,x1_11] = TmpGPUSSACode_SSA(x1)
%First we define the parameters.
k1=10;
k2=1;

%Initialize the time.
t=0;
%Start the SSA.
tstop = 0;   %Next time to print results.
while t<tstop   %Next time to print results.
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  if t<=tstop
    r2w0=rand*w0;
    if r2w0<w1
      x1=x1+(1);
    elseif r2w0<w1+w2
      x1=x1+(-1);
    end
  end
end
x1_1=x1;

tstop = 1;   %Next time to print results.
while t<tstop   %Next time to print results.
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  if t<=tstop
    r2w0=rand*w0;
    if r2w0<w1
      x1=x1+(1);
    elseif r2w0<w1+w2
      x1=x1+(-1);
    end
  end
end
x1_2=x1;

tstop = 2;   %Next time to print results.
while t<tstop   %Next time to print results.
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  if t<=tstop
    r2w0=rand*w0;
    if r2w0<w1
      x1=x1+(1);
    elseif r2w0<w1+w2
      x1=x1+(-1);
    end
  end
end
x1_3=x1;

tstop = 3;   %Next time to print results.
while t<tstop   %Next time to print results.
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  if t<=tstop
    r2w0=rand*w0;
    if r2w0<w1
      x1=x1+(1);
    elseif r2w0<w1+w2
      x1=x1+(-1);
    end
  end
end
x1_4=x1;

tstop = 4;   %Next time to print results.
while t<tstop   %Next time to print results.
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  if t<=tstop
    r2w0=rand*w0;
    if r2w0<w1
      x1=x1+(1);
    elseif r2w0<w1+w2
      x1=x1+(-1);
    end
  end
end
x1_5=x1;

tstop = 5;   %Next time to print results.
while t<tstop   %Next time to print results.
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  if t<=tstop
    r2w0=rand*w0;
    if r2w0<w1
      x1=x1+(1);
    elseif r2w0<w1+w2
      x1=x1+(-1);
    end
  end
end
x1_6=x1;

tstop = 6;   %Next time to print results.
while t<tstop   %Next time to print results.
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  if t<=tstop
    r2w0=rand*w0;
    if r2w0<w1
      x1=x1+(1);
    elseif r2w0<w1+w2
      x1=x1+(-1);
    end
  end
end
x1_7=x1;

tstop = 7;   %Next time to print results.
while t<tstop   %Next time to print results.
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  if t<=tstop
    r2w0=rand*w0;
    if r2w0<w1
      x1=x1+(1);
    elseif r2w0<w1+w2
      x1=x1+(-1);
    end
  end
end
x1_8=x1;

tstop = 8;   %Next time to print results.
while t<tstop   %Next time to print results.
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  if t<=tstop
    r2w0=rand*w0;
    if r2w0<w1
      x1=x1+(1);
    elseif r2w0<w1+w2
      x1=x1+(-1);
    end
  end
end
x1_9=x1;

tstop = 9;   %Next time to print results.
while t<tstop   %Next time to print results.
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  if t<=tstop
    r2w0=rand*w0;
    if r2w0<w1
      x1=x1+(1);
    elseif r2w0<w1+w2
      x1=x1+(-1);
    end
  end
end
x1_10=x1;

tstop = 10;   %Next time to print results.
while t<tstop   %Next time to print results.
  w1=k1;
  w2=k2*x1;
  w0=0+w1+w2;
  t = t-1/w0*log(rand);
  if t<=tstop
    r2w0=rand*w0;
    if r2w0<w1
      x1=x1+(1);
    elseif r2w0<w1+w2
      x1=x1+(-1);
    end
  end
end
x1_11=x1;

