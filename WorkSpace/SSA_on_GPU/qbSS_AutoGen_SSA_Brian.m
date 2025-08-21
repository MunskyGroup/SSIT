function [X]=qbSS_AutoGen_SSA_Brian(x0,N_run,N_split)
%This is an automatically generated MATLAB SSA Program.
%The tools used to generate this file were covered at the.
%2016 q-bio Summer School at the Colorado State University.

Nspec = 4; % Number of species.
Nt = 3; % Number of time points.
n_run = N_run/N_split; % Number of runs per CPU.
X = zeros(Nspec,Nt,N_run); % Initialize matrix of results. 
x1_0_GPU = x0(1)*gpuArray.ones(1,n_run); % Specific Initial Conditions.
x2_0_GPU = x0(2)*gpuArray.ones(1,n_run); % Specific Initial Conditions.
x3_0_GPU = x0(3)*gpuArray.ones(1,n_run); % Specific Initial Conditions.
x4_0_GPU = x0(4)*gpuArray.ones(1,n_run); % Specific Initial Conditions.
parfor i = 1:N_split  % Run loop over different GPU cards.
  [x1_1,x1_2,x1_3,x2_1,x2_2,x2_3,x3_1,x3_2,x3_3,x4_1,x4_2,x4_3] = arrayfun(@qbSS_AutoGen_SSA_Brian_SSA,x1_0_GPU,x2_0_GPU,x3_0_GPU,x4_0_GPU);
  xp1_1(:,i)=gather(x1_1);
  xp1_2(:,i)=gather(x1_2);
  xp1_3(:,i)=gather(x1_3);
  xp2_1(:,i)=gather(x2_1);
  xp2_2(:,i)=gather(x2_2);
  xp2_3(:,i)=gather(x2_3);
  xp3_1(:,i)=gather(x3_1);
  xp3_2(:,i)=gather(x3_2);
  xp3_3(:,i)=gather(x3_3);
  xp4_1(:,i)=gather(x4_1);
  xp4_2(:,i)=gather(x4_2);
  xp4_3(:,i)=gather(x4_3);
end
for k=1:N_split
  X(1,1,(k-1)*n_run+1:k*n_run) = xp1_1(:,k);
  X(1,2,(k-1)*n_run+1:k*n_run) = xp1_2(:,k);
  X(1,3,(k-1)*n_run+1:k*n_run) = xp1_3(:,k);
  X(2,1,(k-1)*n_run+1:k*n_run) = xp2_1(:,k);
  X(2,2,(k-1)*n_run+1:k*n_run) = xp2_2(:,k);
  X(2,3,(k-1)*n_run+1:k*n_run) = xp2_3(:,k);
  X(3,1,(k-1)*n_run+1:k*n_run) = xp3_1(:,k);
  X(3,2,(k-1)*n_run+1:k*n_run) = xp3_2(:,k);
  X(3,3,(k-1)*n_run+1:k*n_run) = xp3_3(:,k);
  X(4,1,(k-1)*n_run+1:k*n_run) = xp4_1(:,k);
  X(4,2,(k-1)*n_run+1:k*n_run) = xp4_2(:,k);
  X(4,3,(k-1)*n_run+1:k*n_run) = xp4_3(:,k);
end


function [x1_1,x1_2,x1_3,x2_1,x2_2,x2_3,x3_1,x3_2,x3_3,x4_1,x4_2,x4_3] = qbSS_AutoGen_SSA_Brian_SSA(x1,x2,x3,x4)
%First we define the parameters.
k1=10;
k2=1;
k3=10;
k4=1;
k5=2;
k6=5;
k7=2;
k8=2;

%Initialize the time.
t=0;
%Start the SSA.
tstop = 0;   %Next time to print results.
while t<tstop;   %Next time to print results.
  w1=k1;
  w2=k2*x1;
  w3=k3*x1;
  w4=k4*x2;
  w5=k5*x2;
  w6=k6*x3;
  w7=k7*x3;
  w8=k8*x4;
  w0=0+w1+w2+w3+w4+w5+w6+w7+w8;
  t = t-1/w0*log(rand);
  if t<=tstop
    r2w0=rand*w0;
    if r2w0<w1;
      x1=x1+(1);
    elseif r2w0<w1+w2;
      x1=x1+(-1);
    elseif r2w0<w1+w2+w3;
      x2=x2+(1);
    elseif r2w0<w1+w2+w3+w4;
      x2=x2+(-1);
    elseif r2w0<w1+w2+w3+w4+w5;
      x2=x2+(-1);
      x3=x3+(1);
    elseif r2w0<w1+w2+w3+w4+w5+w6;
      x3=x3+(-1);
    elseif r2w0<w1+w2+w3+w4+w5+w6+w7;
      x3=x3+(-1);
      x4=x4+(1);
    elseif r2w0<w1+w2+w3+w4+w5+w6+w7+w8;
      x4=x4+(-1);
    end
  end
end
x1_1=x1;
x2_1=x2;
x3_1=x3;
x4_1=x4;

tstop = 250;   %Next time to print results.
while t<tstop;   %Next time to print results.
  w1=k1;
  w2=k2*x1;
  w3=k3*x1;
  w4=k4*x2;
  w5=k5*x2;
  w6=k6*x3;
  w7=k7*x3;
  w8=k8*x4;
  w0=0+w1+w2+w3+w4+w5+w6+w7+w8;
  t = t-1/w0*log(rand);
  if t<=tstop
    r2w0=rand*w0;
    if r2w0<w1;
      x1=x1+(1);
    elseif r2w0<w1+w2;
      x1=x1+(-1);
    elseif r2w0<w1+w2+w3;
      x2=x2+(1);
    elseif r2w0<w1+w2+w3+w4;
      x2=x2+(-1);
    elseif r2w0<w1+w2+w3+w4+w5;
      x2=x2+(-1);
      x3=x3+(1);
    elseif r2w0<w1+w2+w3+w4+w5+w6;
      x3=x3+(-1);
    elseif r2w0<w1+w2+w3+w4+w5+w6+w7;
      x3=x3+(-1);
      x4=x4+(1);
    elseif r2w0<w1+w2+w3+w4+w5+w6+w7+w8;
      x4=x4+(-1);
    end
  end
end
x1_2=x1;
x2_2=x2;
x3_2=x3;
x4_2=x4;

tstop = 500;   %Next time to print results.
while t<tstop;   %Next time to print results.
  w1=k1;
  w2=k2*x1;
  w3=k3*x1;
  w4=k4*x2;
  w5=k5*x2;
  w6=k6*x3;
  w7=k7*x3;
  w8=k8*x4;
  w0=0+w1+w2+w3+w4+w5+w6+w7+w8;
  t = t-1/w0*log(rand);
  if t<=tstop
    r2w0=rand*w0;
    if r2w0<w1;
      x1=x1+(1);
    elseif r2w0<w1+w2;
      x1=x1+(-1);
    elseif r2w0<w1+w2+w3;
      x2=x2+(1);
    elseif r2w0<w1+w2+w3+w4;
      x2=x2+(-1);
    elseif r2w0<w1+w2+w3+w4+w5;
      x2=x2+(-1);
      x3=x3+(1);
    elseif r2w0<w1+w2+w3+w4+w5+w6;
      x3=x3+(-1);
    elseif r2w0<w1+w2+w3+w4+w5+w6+w7;
      x3=x3+(-1);
      x4=x4+(1);
    elseif r2w0<w1+w2+w3+w4+w5+w6+w7+w8;
      x4=x4+(-1);
    end
  end
end
x1_3=x1;
x2_3=x2;
x3_3=x3;
x4_3=x4;

