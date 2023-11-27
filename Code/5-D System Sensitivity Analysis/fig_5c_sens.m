clear all, close all, clc

% Sensitivity analysis on bee-only model with 10% change

% Initial Conditions
x0 = [12500; 5000; 0; 10; 50;];

% Time span (days)
tspan = 0:1:3650; % 10 years

% Zero column vector to not apply varroacide treatment
zero_vec = zeros(5,1);

% gamma_i values (Rate mites kills bees)
gamma_1 = 10^(-7);
gamma_2 = gamma_1;
gamma_3 = 0.0000002;

% Initialzing p-value
p = 0.1; 

% Initializing zero matrices to store LB/UB avgs. and percent changes
% Column 1 = LB avgs./percentages, Column 2 = UB avgs./percentages, Rows = parameters
% LB = lower bound, UB = upper bound
avg_populations = zeros(22,2);
percent_changes = zeros(22,2);

i = 1; % Initializing iterator value

% Increasing tolerances for ode45
options = odeset('RelTol',1e-3,'AbsTol',1e-5);

% Solving system with default (Table 1) parameter values 
[t,x] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Calculating avgerage total bee population
avg_pop_table_1 = sum(x(:,1)+x(:,2)+x(:,3))/length(tspan)

%% Running simulation with +/- 10% for each parameter

% beta_1

% Solving system with lower bound beta_1 value
[t,x_LB] = ode45(@(t,x) system_eqs(t,x,(beta_1_fcn(t) - 0.1*beta_1_fcn(t)),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound beta_1 value
[t,x_UB] = ode45(@(t,x) system_eqs(t,x,(beta_1_fcn(t) + 0.1*beta_1_fcn(t)),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% beta_2

% Solving system with lower bound beta_2 value
[t,x_LB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),(beta_2_fcn(t) - 0.1*beta_2_fcn(t)),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound beta_2 value
[t,x_UB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),(beta_2_fcn(t) + 0.1*beta_2_fcn(t)),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% beta_3

% Solving system with lower bound beta_3 value
[t,x_LB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),(beta_3_fcn(t) - 0.1*beta_3_fcn(t)),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound beta_3 value
[t,x_UB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),(beta_3_fcn(t) + 0.1*beta_3_fcn(t)),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% d_1

% Solving system with lower bound d_1 value
[t,x_LB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),(d_1_fcn(t) - 0.1*d_1_fcn(t)), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound d_1 value
[t,x_UB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),(d_1_fcn(t) + 0.1*d_1_fcn(t)), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% d_2

% Solving system with lower bound d_2 value
[t,x_LB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
(d_2_fcn(t) - 0.1*d_2_fcn(t)),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound d_2 value
[t,x_UB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
(d_2_fcn(t) + 0.1*d_2_fcn(t)),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% d_3

% Solving system with lower bound d_3 value
[t,x_LB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),(d_3_fcn(t) - 0.1*d_3_fcn(t)),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound d_3 value
[t,x_UB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),(d_3_fcn(t) + 0.1*d_3_fcn(t)),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% mu

% Solving system with lower bound mu value
[t,x_LB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),(mu_fcn(t) - 0.1*mu_fcn(t)),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound mu value
[t,x_UB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),(mu_fcn(t) + 0.1*mu_fcn(t)),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% little_k

% Solving system with lower bound little_k value
[t,x_LB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),(little_k_fcn(t) - 0.1*little_k_fcn(t)),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound little_k value
[t,x_UB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),(little_k_fcn(t) + 0.1*little_k_fcn(t)),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% r

% Solving system with lower bound r value
[t,x_LB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),(r_fcn(t) - 0.1*r_fcn(t)),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound r value
[t,x_UB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),(r_fcn(t) + 0.1*r_fcn(t)),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% alpha

% Solving system with lower bound alpha value
[t,x_LB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),(alpha_fcn(t) - 0.1*alpha_fcn(t)),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound alpha value
[t,x_UB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),(alpha_fcn(t) + 0.1*alpha_fcn(t)),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% K

% Solving system with lower bound K value
[t,x_LB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),(K_fcn(t) - 0.1*K_fcn(t)),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound K value
[t,x_UB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),(K_fcn(t) + 0.1*K_fcn(t)),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options); %#ok<*ASGLU> 

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% sigma_1

% Solving system with lower bound sigma_1 value
[t,x_LB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),(sigma_1_fcn(t) - 0.1*sigma_1_fcn(t)), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound sigma_1 value
[t,x_UB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),(sigma_1_fcn(t) + 0.1*sigma_1_fcn(t)), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% sigma_2

% Solving system with lower bound sigma_2 value
[t,x_LB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
(sigma_2_fcn(t) - 0.1*sigma_2_fcn(t)),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound sigma_2 value
[t,x_UB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
(sigma_2_fcn(t) + 0.1*sigma_2_fcn(t)),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% gamma_1

% Solving system with lower bound gamma_1 value
[t,x_LB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),(gamma_1 - 0.1*gamma_1),gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound gamma_1 value
[t,x_UB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),(gamma_1 + 0.1*gamma_1),gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% gamma_2

% Solving system with lower bound gamma_2 value
[t,x_LB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,(gamma_2 - 0.1*gamma_2),gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound gamma_2 value
[t,x_UB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,(gamma_2 + 0.1*gamma_2),gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% gamma_3

% Solving system with lower bound gamma_3 value
[t,x_LB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,(gamma_3 - 0.1*gamma_3),delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound gamma_3 value
[t,x_UB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,(gamma_3 + 0.1*gamma_3),delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% delta_1

% Solving system with lower bound delta_1 value
[t,x_LB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,(delta_1_fcn(t) - 0.1*delta_1_fcn(t)),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound delta_1 value
[t,x_UB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,(delta_1_fcn(t) + 0.1*delta_1_fcn(t)),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% delta_2

% Solving system with lower bound delta_2 value
[t,x_LB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),(delta_2_fcn(t) - 0.1*delta_2_fcn(t)),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound delta_2 value
[t,x_UB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),(delta_2_fcn(t) + 0.1*delta_2_fcn(t)),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% delta_3

% Solving system with lower bound delta_3 value
[t,x_LB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),(delta_3_fcn(t) - 0.1*delta_3_fcn(t)),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound delta_3 value
[t,x_UB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),(delta_3_fcn(t) + 0.1*delta_3_fcn(t)),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% delta_4

% Solving system with lower bound delta_4 value
[t,x_LB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),(delta_4_fcn(t) - 0.1*delta_4_fcn(t)), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound delta_4 value
[t,x_UB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),(delta_4_fcn(t) + 0.1*delta_4_fcn(t)), ...
delta_5_fcn(t),p),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% delta_5

% Solving system with lower bound delta_5 value
[t,x_LB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t),(delta_5_fcn(t) - 0.1*delta_5_fcn(t)), ...
p),tspan,x0,options);

% Solving system with upper bound delta_5 value
[t,x_UB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t),(delta_5_fcn(t) + 0.1*delta_5_fcn(t)), ...
p),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% p

% Solving system with lower bound p value
[t,x_LB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t),delta_5_fcn(t), ...
(p - 0.1*p)),tspan,x0,options);

% Solving system with upper bound p value
[t,x_UB] = ode45(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t),delta_5_fcn(t), ...
(p + 0.1*p)),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% Save percent_changes matrix as an Excel file
xlswrite('fig_5c_sens_nonnegative.xlsx',percent_changes)