% Testing the ODE solver
clear all,clc

addpath(genpath('src'));

tspan = linspace(0, 10, 10);
x0 = [10 0 0]; % initial state
stoichMatrix = [1 0 0;
    0 1 0;
    0 0 1;
    -1 0 0;
    0 -1 0;
    0 0 -1]';

pars = { 'a', 5;
    'kn0', 0;
    'kn1', 25;
    'n', 6;
    'g', 1};
inputs = {};

propens = { ...
    'kn0+kn1*(1/(1+a*(x2^n)))',...
    'kn0+kn1*(1/(1+a*(x3^n)))',...
    'kn0+kn1*(1/(1+a*(x1^n)))',...
    'g*x1',...
    'g*x2',...
    'g*x3'
    };
[t_ode,ode_solutions] = ssit.moments.solve_ode(x0, tspan, stoichMatrix, propens, pars, inputs)

plot(t_ode,ode_solutions)


