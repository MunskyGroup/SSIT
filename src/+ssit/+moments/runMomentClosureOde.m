clear all
close all

% example problem for linear propensities
S = [1 0 0; -1 0 0;0 1 0;0 -1 0;0 0 1;0 0 -1]'; % stoich matrix
syms x1 x2 x3 real % create symbolic variables for species
w = [50; x1; 12; 0.5*x2; 3; 0.2*x3]; % propensity functions
x = [x1;x2;x3]; % species
n = size(x,1); % number of species

% Call function to generate moment equations
ssit.moments.writeFunForMomentsODE(S,w,x,'tmpMomentOdeRHS')

% run ode for moments equations
x0 = [20;40;60];
v0 = [x0(1) x0(2) x0(3),...
    x0(1)+x0(1)^2 x0(1)*x0(2) x0(1)*x0(3),...
    x0(2)+x0(2)^2 x0(2)*x0(3),...
    x0(3)+x0(3)^2]'; % create initial conditions
momOde = @(t,v)tmpMomentOdeRHS(v); % call on function written for the ode
[t,vvt] = ode45(momOde,linspace(0,10,100),v0); % solve the ode from time 0 to 10

% plot results
figure(1); clf;
plot(t,vvt(:,1:n)); hold on % plot mean of each species
plot(t,vvt(:,n+1)-vvt(:,1).^2,'r--'); % plot variance of species 1
plot(t,vvt(:,n+n+1)-vvt(:,2).^2,'g--'); % plot variance of species 2
plot(t,vvt(:,n+n+(n-1)+1)-vvt(:,3).^2,'c--'); % plot variance of species 3
plot(t,vvt(:,n+2)-vvt(:,1).*vvt(:,2),'m--'); % plot the covariance of species 1 and 2
legend('\mu_1','\mu_2','\mu_3','\sigma_{11}','\sigma_{22}','\sigma_{33}','\sigma_{12}')

%%
clear all
% example problem for non-linear propensities
S = [1; -1]'; % stoich matrix
syms x1 real
w = [50; x1*(x1-1)]; % propensity functions
x = [x1]; % species
n = size(x,1); % number of species

% Call function to generate moment equations
ssit.moments.writeFunForMomentsODE(S,w,x,'tmpMomentOdeRHS')

% run ode for moments equations
x0 = [20];
v0 = [x0(1),x0(1)+x0(1)^2]'; % create initial conditions

momOde = @(t,v)tmpMomentOdeRHS(v); % call on function written for the ode
[t,vvt] = ode45(momOde,linspace(0,10,100),v0); % solve the ode from time 0 to 10

% plot results
figure(1); clf;
plot(t,vvt(:,1:n)); hold on % plot mean of each species
% plot(t,vvt(:,n+1)-vvt(:,1).^2,'r--'); % plot variance of species 1
legend('\mu_1','\sigma_{11}')

%% SIR model check

clear all

S = [-1 1 0; 0 -1 1]'; % stoich matrix
syms x1 x2 x3 real
w = [0.0057*x1*x2; 0.330*x2]; % propensity functions
x = [x1;x2;x3]; % species
n = size(x,1); % number of species

% Call function to generate moment equations
ssit.moments.writeFunForMomentsODE(S,w,x,'tmpMomentOdeRHS')

% run ode for moments equations
x0 = [99;5;0];
v0 = [x0(1) x0(2) x0(3),...
    x0(1)^2 x0(1)*x0(2) x0(1)*x0(3),...
    x0(2)^2 x0(2)*x0(3),...
    x0(3)^2]'; % create initial conditions
momOde = @(t,v)tmpMomentOdeRHS(v); % call on function written for the ode
[t,vvt] = ode45(momOde,linspace(0,60,100),v0); % solve the ode from time 0 to 10

% plot results
figure(1); clf;
plot(t,vvt(:,1:n)); hold on % plot mean of each species
% plot(t,vvt(:,n+1)-vvt(:,1).^2,'r--'); % plot variance of species 1
% plot(t,vvt(:,n+n+1)-vvt(:,2).^2,'g--'); % plot variance of species 2
% plot(t,vvt(:,n+n+(n-1)+1)-vvt(:,3).^2,'c--'); % plot variance of species 3
% plot(t,vvt(:,n+2)-vvt(:,1).*vvt(:,2),'m--'); % plot the covariance of species 1 and 2
legend('\mu_1','\mu_2','\mu_3','\sigma_{11}','\sigma_{22}','\sigma_{33}','\sigma_{12}')
%% poisson model check

clear all

S = [1;-1]'; % stoich matrix
syms x1  real
w = [40; x1]; % propensity functions
x = [x1]; % species
n = size(x,1); % number of species

% Call function to generate moment equations
ssit.moments.writeFunForMomentsODE(S,w,x,'tmpMomentOdeRHS')

% run ode for moments equations
x0 = [0;0;0];
v0 = [x0(1),x0(1)+x0(1)^2]'; % create initial conditions
momOde = @(t,v)tmpMomentOdeRHS(v); % call on function written for the ode
[t,vvt] = ode45(momOde,linspace(0,10,100),v0); % solve the ode from time 0 to 10

% plot results
figure(1); clf;
plot(t,vvt(:,1:n)); hold on % plot mean of each species
plot(t,vvt(:,n+1)-vvt(:,1).^2,'r--'); % plot variance of species 1
legend('\mu_1','\sigma_{11}')

%% autoregulation model check

clear all

S = [1 0; -1 0; 0 1; 0 -1]'; % stoich matrix
syms x1 x2 real
w = [5/(1+0*x2); 1*x1; 5*x1; 1*x2]; % propensity functions
x = [x1;x2]; % species
n = size(x,1); % number of species

% Call function to generate moment equations
ssit.moments.writeFunForMomentsODE(S,w,x,'tmpMomentOdeRHS')

% run ode for moments equations
x0 = [10;10];
v0 = [x0(1) x0(2),...
    x0(1)+x0(1)^2 x0(1)*x0(2),...
    x0(2)+x0(2)^2]'; % create initial conditions
momOde = @(t,v)tmpMomentOdeRHS(v); % call on function written for the ode
[t,vvt] = ode45(momOde,linspace(0,10,100),v0); % solve the ode from time 0 to 10

% plot results
figure(1); clf;
plot(t,vvt(:,1:n)); hold on % plot mean of each species
plot(t,vvt(:,n+1)-vvt(:,1).^2,'r--'); % plot variance of species 1
plot(t,vvt(:,n+n+1)-vvt(:,2).^2,'g--'); % plot variance of species 2
% plot(t,vvt(:,n+2)-vvt(:,1).*vvt(:,2),'m--'); % plot the covariance of species 1 and 2
legend('\mu_1','\mu_2','\sigma_{11}','\sigma_{22}','\sigma_{12}')

%% gene expression negative feedback loop model check

clear all

S = [1 0;-1 0; 0 1; 0 -1]'; % stoich matrix
syms x1 x2 real
w = [0.05*(1-x1); 0.001*x1*x2^2; 2288*x1*70; 1*x2]; % propensity functions
x = [x1;x2]; % species
n = size(x,1); % number of species

% Call function to generate moment equations
ssit.moments.writeFunForMomentsODE(S,w,x,'tmpMomentOdeRHS')

% run ode for moments equations
x0 = [1;20];
v0 = [x0(1) x0(2),...
    x0(1)+x0(1)^2 x0(1)*x0(2),...
    x0(2)+x0(2)^2]'; % create initial conditions
momOde = @(t,v)tmpMomentOdeRHS(v); % call on function written for the ode
[t,vvt] = ode45(momOde,linspace(0,10,100),v0); % solve the ode from time 0 to 10

% plot results
figure(1); clf;
plot(t,vvt(:,1:n)); hold on % plot mean of each species
% plot(t,vvt(:,n+1)-vvt(:,1).^2,'r--'); % plot variance of species 1
% plot(t,vvt(:,n+n+1)-vvt(:,2).^2,'g--'); % plot variance of species 2
% plot(t,vvt(:,n+2)-vvt(:,1).*vvt(:,2),'m--'); % plot the covariance of species 1 and 2
legend('\mu_1','\mu_2','\sigma_{11}','\sigma_{22}','\sigma_{12}')

% source: https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=7226883&casa_token=0CZze2V_0ckAAAAA:7ic9GZMhqCidYLsDnq3SNteEjTCA279zvwlXFg-QpWGIWZhgKVGiwz-p-j68Ea3A_IwMZRycNxk&tag=1

%% two-gene circuit with acitivator and repressor model check

clear all

S = [1 0 0 0;-1 0 0 0; 0 1 0 0; 0 -1 0 0; 0 0 1 0; 0 0 -1 0; 0 0 0 1; 0 0 0 -1]'; % stoich matrix
syms x1 x2 x3 x4 real
w = [0.05*(1-x1)*x4; 4*x1; 0.3*(1-x2); 0.05*x2*x3; 50*x1*5; 1*x3; 350*x2*5; 1*x4]; % propensity functions
x = [x1;x2;x3;x4]; % species
n = size(x,1); % number of species

% Call function to generate moment equations
ssit.moments.writeFunForMomentsODE(S,w,x,'tmpMomentOdeRHS')

% run ode for moments equations
x0 = [5;5;5;5];
v0 = [x0(1) x0(2) x0(3) x0(4),...
    x0(1)+x0(1)^2 x0(1)*x0(2) x0(1)*x0(3) x0(1)*x0(4),...
    x0(2)+x0(2)^2 x0(2)*x0(3) x0(2)*x0(4),...
    x0(3)+x0(3)^2 x0(3)*x0(4),...
    x0(4)+x0(4)^2]'; % create initial conditions
momOde = @(t,v)tmpMomentOdeRHS(v); % call on function written for the ode
[t,vvt] = ode45(momOde,linspace(0,10,100),v0); % solve the ode from time 0 to 10

% plot results
figure(1); clf;
plot(t,vvt(:,1:n)); hold on % plot mean of each species
% plot(t,vvt(:,n+1)-vvt(:,1).^2,'r--'); % plot variance of species 1
% plot(t,vvt(:,n+n+1)-vvt(:,2).^2,'g--'); % plot variance of species 2
% plot(t,vvt(:,n+n+(n-1)+1)-vvt(:,3).^2,'c--'); % plot variance of species 3
% plot(t,vvt(:,n+n+(n-1)+3)-vvt(:,4).^2,'c--'); % plot variance of species 4
% plot(t,vvt(:,n+2)-vvt(:,1).*vvt(:,2),'m--'); % plot the covariance of species 1 and 2
legend('\mu_1','\mu_2','\mu_3','\mu_4','\sigma_{11}','\sigma_{22}','\sigma_{33}','\sigma_{12}')


% source: https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=7226883&casa_token=0CZze2V_0ckAAAAA:7ic9GZMhqCidYLsDnq3SNteEjTCA279zvwlXFg-QpWGIWZhgKVGiwz-p-j68Ea3A_IwMZRycNxk&tag=1

%% enzymatic reaction model check

clear all

S = [1 0; -1 0; 0 1; 0 -1]'; % stoich matrix
syms x1 x2 real
w = [3; x1; 2+1*x1; x2]; % propensity functions
x = [x1;x2]; % species
n = size(x,1); % number of species

% Call function to generate moment equations
ssit.moments.writeFunForMomentsODE(S,w,x,'tmpMomentOdeRHS')

% run ode for moments equations
x0 = [0;0];
v0 = [x0(1) x0(2),...
    x0(1)+x0(1)^2 x0(1)*x0(2),...
    x0(2)+x0(2)^2]'; % create initial conditions
momOde = @(t,v)tmpMomentOdeRHS(v); % call on function written for the ode
[t,vvt] = ode45(momOde,linspace(0,10,100),v0); % solve the ode from time 0 to 10

% plot results
figure(1); clf;
plot(t,vvt(:,1:n)); hold on % plot mean of each species
plot(t,vvt(:,n+1)-vvt(:,1).^2,'r--'); % plot variance of species 1
plot(t,vvt(:,n+n+1)-vvt(:,2).^2,'g--'); % plot variance of species 2
% plot(t,vvt(:,n+2)-vvt(:,1).*vvt(:,2),'m--'); % plot the covariance of species 1 and 2
legend('\mu_1','\mu_2','\sigma_{11}','\sigma_{22}','\sigma_{12}')

%%  chemical kinetics model check

clear all

S = [40 0; -2 0; 0 15; 0 -1]'; % stoich matrix
syms x1 x2 real
w = [100; 0.8/1*x1^2; 0.02/1*x1^2; 15*x2]; % propensity functions
x = [x1;x2]; % species
n = size(x,1); % number of species

% Call function to generate moment equations
ssit.moments.writeFunForMomentsODE(S,w,x,'tmpMomentOdeRHS')

% run ode for moments equations
x0 = [40;40];
v0 = [x0(1) x0(2),...
    x0(1)+x0(1)^2 x0(1)*x0(2),...
    x0(2)+x0(2)^2]'; % create initial conditions
momOde = @(t,v)tmpMomentOdeRHS(v); % call on function written for the ode
[t,vvt] = ode45(momOde,linspace(0,0.15,100),v0); % solve the ode from time 0 to 10

% plot results
figure(1); clf;
plot(t,vvt(:,1:n)); hold on % plot mean of each species
% plot(t,vvt(:,n+1)-vvt(:,1).^2,'r--'); % plot variance of species 1
% plot(t,vvt(:,n+n+1)-vvt(:,2).^2,'g--'); % plot variance of species 2
% plot(t,vvt(:,n+2)-vvt(:,1).*vvt(:,2),'m--'); % plot the covariance of species 1 and 2
legend('\mu_1','\mu_2','\sigma_{11}','\sigma_{22}','\sigma_{12}')

% source: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.159.296&rep=rep1&type=pdf

%%  dimerisation model check

clear all

S = [-2 1; 2 -1]'; % stoich matrix
syms x1 x2 real
w = [(0.00166*x1*(x1-1))/2; (0.2*(301-x1))/2]; % propensity functions
x = [x1;x2]; % species
n = size(x,1); % number of species

% Call function to generate moment equations
ssit.moments.writeFunForMomentsODE(S,w,x,'tmpMomentOdeRHS')

% run ode for moments equations
x0 = [301;0];
v0 = [x0(1) x0(2),...
    x0(1)+x0(1)^2 x0(1)*x0(2),...
    x0(2)+x0(2)^2]'; % create initial conditions
momOde = @(t,v)tmpMomentOdeRHS(v); % call on function written for the ode
[t,vvt] = ode45(momOde,linspace(0,10,100),v0); % solve the ode from time 0 to 10

% plot results
figure(1); clf;
plot(t,vvt(:,1:n)); hold on % plot mean of each species
% plot(t,vvt(:,n+1)-vvt(:,1).^2,'r--'); % plot variance of species 1
% plot(t,vvt(:,n+n+1)-vvt(:,2).^2,'g--'); % plot variance of species 2
% plot(t,vvt(:,n+2)-vvt(:,1).*vvt(:,2),'m--'); % plot the covariance of species 1 and 2
legend('\mu_1','\mu_2','\sigma_{11}','\sigma_{22}','\sigma_{12}')

% source: https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4762250&casa_token=avIXhRCUziYAAAAA:s5dRwNBU_dn95phRgrWjWWq9Vi4fL0--nC51-HqC8d_yI9YU8Qk8KRPZOSv4MEpU7rxlreT6WSw