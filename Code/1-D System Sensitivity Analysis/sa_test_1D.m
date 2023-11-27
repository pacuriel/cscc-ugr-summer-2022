clear all, close all, clc

% Attempting local sensitivity analysis on default model on 1D system

% Initial Conditions
x0 = [12500;];

% Time span (days)
tspan = 0:1:500; % 10 years

% Increasing tolerances for ode45
% options = odeset('RelTol',1e-3,'AbsTol',1e-5);

d_1 = 0.02272; % Spring value

% d_1 +/- 10 percent 
d_1_LB = d_1 - d_1*0.1;
d_1_UB = d_1 + d_1*0.1;

delta_d_1 = 10^(-6);
 
new_d_1 = d_1_LB; % Initializing to d_1_LB


avg_values = zeros(1,numel(tspan));

% Solving system with default (Table 1) parameter values 




% Solving system with Spring, LB, and UB d_1 values
[t,x] = ode45(@(t,x) sys_eqs_1D(t,x,d_1),tspan,x0);
[t_LB,x_LB] = ode45(@(t,x) sys_eqs_1D(t,x,d_1_LB),tspan,x0);
[t_UB,x_UB] = ode45(@(t,x) sys_eqs_1D(t,x,d_1_UB),tspan,x0);

% Calculating avgerages
avg_spring = sum(x(:,1))/length(tspan)
avg_LB = sum(x_LB(:,1))/length(tspan)
avg_UB = sum(x_UB(:,1))/length(tspan)

% Calculating percent change in avgs. 
percent_change_LB = 100*((avg_LB/avg_spring) - 1)
percent_change_UB = 100*((avg_UB/avg_spring) - 1)

i = 1; % Iterator value

% while loop to calculate avgs of sims
% while (new_d_1 <= d_1_UB)
% 
%     % Solving w/ new_d_1 val
%     [t1,x1] = ode45(@(t,x) sys_eqs_1D(t,x,new_d_1),tspan,x0);
%    
%     % Getting avg values
%     avg_values(i) = sum(x1(:,1))/length(tspan);
% 
%     % Incrementing new_d_1 and i
%     new_d_1 = new_d_1 + delta_d_1;
%     i = i + 1;
% end

 plot(t,x(:,1),t_LB,x_LB(:,1),t_UB,x_UB(:,1),LineWidth=1)
% legend('Spring','LB','UB')