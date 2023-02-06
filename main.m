clear
clc

global problem

rng(123);

% Choose the problem to be solved

problem = 'BK1';

% Initial parameters

[n,m,l,u,x0] = inip;

% Define the uncertainty parameter

delta = 2/100 * norm(x0);

% Set problem data

[dimA,A,b] = datas(n,m,delta);

% Call the solvers

[x,info] = CondG(n,m,l,u,x0,dimA,A,b);

[x,info] = ProxGrad(n,m,l,u,x0,dimA,A,b);