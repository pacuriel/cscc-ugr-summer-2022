clear all, close all, clc

% Sensitivity analysis on bee-only model with 10% change

% Initial Conditions
x0 = [12500; 5000; 0; 0; 0;];

% Time span (days)
tspan = 0:1:3650; % 10 years

% Zero column vector for varroacide treatment
zero_vec = zeros(5,1);

% Initializing zero matrices to store LB/UB avgs. and percent changes
% Column 1 = LB avgs./percentages, Column 2 = UB avgs./percentages, Rows = parameters
% LB = lower bound, UB = upper bound
avg_populations = zeros(6,2);
percent_changes = zeros(6,2);

i = 1; % Initializing iterator value

% Increasing tolerances for ode45
options = odeset('RelTol',1e-3,'AbsTol',1e-5);

% Solving system with default (Table 1) parameter values 
[t,x] = ode45(@(t,x) sys_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_i_fcn(t),zero_vec),tspan,x0,options);

% Calculating avgerage total bee population
avg_pop_table_1 = sum(x(:,1)+x(:,2))/length(tspan)

% d_1

% Solving system with lower bound d_1 value
[t,x_LB] = ode45(@(t,x) sys_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),(d_1_fcn(t) - 0.1*d_1_fcn(t)), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_i_fcn(t),zero_vec),tspan,x0,options);

% Solving system with upper bound d_1 value
[t,x_UB] = ode45(@(t,x) sys_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),(d_1_fcn(t) + 0.1*d_1_fcn(t)), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_i_fcn(t),zero_vec),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% d_2

% Solving system with lower bound d_2 value
[t,x_LB] = ode45(@(t,x) sys_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
(d_2_fcn(t) - 0.1*d_2_fcn(t)),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_i_fcn(t),zero_vec),tspan,x0,options);

% Solving system with upper bound d_2 value
[t,x_UB] = ode45(@(t,x) sys_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
(d_2_fcn(t) + 0.1*d_2_fcn(t)),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_i_fcn(t),zero_vec),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% mu

% Solving system with lower bound mu value
[t,x_LB] = ode45(@(t,x) sys_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),(mu_fcn(t) - 0.1*mu_fcn(t)),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_i_fcn(t),zero_vec),tspan,x0,options);

% Solving system with upper bound mu value
[t,x_UB] = ode45(@(t,x) sys_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),(mu_fcn(t) + 0.1*mu_fcn(t)),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_i_fcn(t),zero_vec),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% K

% Solving system with lower bound mu value
[t,x_LB] = ode45(@(t,x) sys_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),(K_fcn(t) - 0.1*K_fcn(t)),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_i_fcn(t),zero_vec),tspan,x0,options);

% Solving system with upper bound mu value
[t,x_UB] = ode45(@(t,x) sys_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),(K_fcn(t) + 0.1*K_fcn(t)),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_i_fcn(t),zero_vec),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% sigma_1

% Solving system with lower bound mu value
[t,x_LB] = ode45(@(t,x) sys_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),(sigma_1_fcn(t) - 0.1*sigma_1_fcn(t)), ...
sigma_2_fcn(t),gamma_i_fcn(t),zero_vec),tspan,x0,options);

% Solving system with upper bound mu value
[t,x_UB] = ode45(@(t,x) sys_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),(sigma_1_fcn(t) + 0.1*sigma_1_fcn(t)), ...
sigma_2_fcn(t),gamma_i_fcn(t),zero_vec),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2))/length(tspan))/avg_pop_table_1) - 1);

i = i + 1; % Incrementing iterator

% sigma_2

% Solving system with lower bound mu value
[t,x_LB] = ode45(@(t,x) sys_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
(sigma_2_fcn(t) - 0.1*sigma_2_fcn(t)),gamma_i_fcn(t),zero_vec),tspan,x0,options);

% Solving system with upper bound mu value
[t,x_UB] = ode45(@(t,x) sys_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
(sigma_2_fcn(t) + 0.1*sigma_2_fcn(t)),gamma_i_fcn(t),zero_vec),tspan,x0,options);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2))/length(tspan)

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*(((sum(x_LB(:,1)+x_LB(:,2))/length(tspan))/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*(((sum(x_UB(:,1)+x_UB(:,2))/length(tspan))/avg_pop_table_1) - 1)

i = i + 1; % Incrementing iterator

% Save percent_changes matrix as an Excel file
%xlswrite('bee_only_model_SA.xlsx',percent_changes)