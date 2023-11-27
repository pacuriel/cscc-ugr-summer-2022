clear all, close all, clc

% Attempting local sensitivity analysis on  5-D system (bee-only model)

% Initial Conditions
x0 = [12500; 0; 0; 0; 0;];

% Time span (days)
tspan = 0:1:500; % 10 years

% Zero column vector for varroacide treatment
zero_vec = zeros(5,1);
zero_vec2 = zeros(3,1);

% Increasing tolerances for ode45
options = odeset('RelTol',1e-3,'AbsTol',1e-5);

% Solving system with default (Table 1) parameter values 
[t,x] = ode45(@(t,x) sys_eqs(t,x,0,0,0,0.02272, ...
0,0,0,0,0,0,0,0,0,zero_vec2,zero_vec),tspan,x0,options);

plot(t,x(:,1))
