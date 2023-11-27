% Function sensitivity analysis on complete model w/ periodic coefficients
function [avg_populations,percent_changes,days_survived] = sens_periodic_5D(x0,tspan,percent_change)

% Zero column vector to not apply varroacide treatment
zero_vec = zeros(5,1);

% gamma_i values (Rate mites kills bees)
gamma_1 = 10^(-7);
gamma_2 = gamma_1;
gamma_3 = 0.0000002;

% Initialzing p-value (homing failure)
p = 0.2; 

% Initializing zero matrices to store LB/UB avgs. and percent changes
% Col1 = LB avgs./percentages, Col2 = UB avgs./percentages, Col3 = Total change, Rows = parameters
% LB = lower bound, UB = upper bound
avg_populations = zeros(22,2);
percent_changes = zeros(22,3);
days_survived = zeros(22,3);

i = 1; % Initializing iterator value
epsilon = 10^(-5); % Epsilon to track colony failure

% Increasing tolerances for ode15s
options = odeset('RelTol',1e-9,'AbsTol',1e-8,'NonNegative',1:5);

% Solving system with default (Table 1) parameter values 
[t,x] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Storing days survived by colony
if (find(x(:,1) <= epsilon,1,'first'))
    days_survived_table_1(i,1) = find(x(:,1) <= epsilon,1,'first');
else
    days_survived_table_1(i,1) = size(x,1);
end

% Calculating avgerage total bee population
avg_pop_table_1 = sum(x(:,1)+x(:,2)+x(:,3))/length(tspan)


%Plotting system solutions
% figure(1);
% hold on
% plot(t,x(:,1),'b',LineWidth=1);
% plot(t,x(:,2),'k--',LineWidth=1);
% plot(t,x(:,3),'-og','MarkerSize',2,'MarkerFaceColor','g',LineWidth=1)
% plot(t,x(:,4),'-.','Color',[0.4940 0.1840 0.556],LineWidth=1);
% plot(t,x(:,5),'r:',LineWidth=1);
% xlabel('Time (days)')
% ylabel('Population')
% legend('x_h','x_f','y','m','n');
% grid on
%% Running simulation with +/- 10% for each parameter

% beta_1

% Solving system with lower bound beta_1 value
[t,x_LB] = ode15s(@(t,x) system_eqs(t,x,(beta_1_fcn(t) - percent_change*beta_1_fcn(t)),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound beta_1 value
[t,x_UB] = ode15s(@(t,x) system_eqs(t,x,(beta_1_fcn(t) + percent_change*beta_1_fcn(t)),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_table_1);
i = i + 1; % Incrementing iterator

% beta_2

% Solving system with lower bound beta_2 value
[t,x_LB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),(beta_2_fcn(t) - percent_change*beta_2_fcn(t)),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound beta_2 value
[t,x_UB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),(beta_2_fcn(t) + percent_change*beta_2_fcn(t)),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_table_1);
i = i + 1; % Incrementing iterator

% beta_3

% Solving system with lower bound beta_3 value
[t,x_LB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),(beta_3_fcn(t) - percent_change*beta_3_fcn(t)),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound beta_3 value
[t,x_UB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),(beta_3_fcn(t) + percent_change*beta_3_fcn(t)),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_table_1);
i = i + 1; % Incrementing iterator

% d_1

% Solving system with lower bound d_1 value
[t,x_LB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),(d_1_fcn(t) - percent_change*d_1_fcn(t)), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound d_1 value
[t,x_UB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),(d_1_fcn(t) + percent_change*d_1_fcn(t)), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_table_1);
i = i + 1; % Incrementing iterator

% d_2

% Solving system with lower bound d_2 value
[t,x_LB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
(d_2_fcn(t) - percent_change*d_2_fcn(t)),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound d_2 value
[t,x_UB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
(d_2_fcn(t) + percent_change*d_2_fcn(t)),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_table_1);
i = i + 1; % Incrementing iterator

% d_3

% Solving system with lower bound d_3 value
[t,x_LB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),(d_3_fcn(t) - percent_change*d_3_fcn(t)),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound d_3 value
[t,x_UB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),(d_3_fcn(t) + percent_change*d_3_fcn(t)),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_table_1);
i = i + 1; % Incrementing iterator

% mu

% Solving system with lower bound mu value
[t,x_LB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),(mu_fcn(t) - percent_change*mu_fcn(t)),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound mu value
[t,x_UB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),(mu_fcn(t) + percent_change*mu_fcn(t)),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_table_1);
i = i + 1; % Incrementing iterator

% little_k

% Solving system with lower bound little_k value
[t,x_LB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),(little_k_fcn(t) - percent_change*little_k_fcn(t)),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound little_k value
[t,x_UB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),(little_k_fcn(t) + percent_change*little_k_fcn(t)),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_table_1);
i = i + 1; % Incrementing iterator

% r

% Solving system with lower bound r value
[t,x_LB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),(r_fcn(t) - percent_change*r_fcn(t)),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound r value
[t,x_UB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),(r_fcn(t) + percent_change*r_fcn(t)),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_table_1);
i = i + 1; % Incrementing iterator

% alpha

% Solving system with lower bound alpha value
[t,x_LB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),(alpha_fcn(t) - percent_change*alpha_fcn(t)),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound alpha value
[t,x_UB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),(alpha_fcn(t) + percent_change*alpha_fcn(t)),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_table_1);
i = i + 1; % Incrementing iterator

% K

% Solving system with lower bound K value
[t,x_LB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),(K_fcn(t) - percent_change*K_fcn(t)),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound K value
[t,x_UB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),(K_fcn(t) + percent_change*K_fcn(t)),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options); %#ok<*ASGLU> 

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_table_1);
i = i + 1; % Incrementing iterator

% sigma_1

% Solving system with lower bound sigma_1 value
[t,x_LB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),(sigma_1_fcn(t) - percent_change*sigma_1_fcn(t)), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound sigma_1 value
[t,x_UB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),(sigma_1_fcn(t) + percent_change*sigma_1_fcn(t)), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_table_1);
i = i + 1; % Incrementing iterator

% sigma_2

% Solving system with lower bound sigma_2 value
[t,x_LB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
(sigma_2_fcn(t) - percent_change*sigma_2_fcn(t)),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound sigma_2 value
[t,x_UB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
(sigma_2_fcn(t) + percent_change*sigma_2_fcn(t)),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_table_1);
i = i + 1; % Incrementing iterator

% gamma_1

% Solving system with lower bound gamma_1 value
[t,x_LB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),(gamma_1 - percent_change*gamma_1),gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound gamma_1 value
[t,x_UB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),(gamma_1 + percent_change*gamma_1),gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_table_1);
i = i + 1; % Incrementing iterator

% gamma_2

% Solving system with lower bound gamma_2 value
[t,x_LB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,(gamma_2 - percent_change*gamma_2),gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound gamma_2 value
[t,x_UB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,(gamma_2 + percent_change*gamma_2),gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_table_1);
i = i + 1; % Incrementing iterator

% gamma_3

% Solving system with lower bound gamma_3 value
[t,x_LB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,(gamma_3 - percent_change*gamma_3),delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound gamma_3 value
[t,x_UB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,(gamma_3 + percent_change*gamma_3),delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_table_1);
i = i + 1; % Incrementing iterator

% delta_1

% Solving system with lower bound delta_1 value
[t,x_LB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,(delta_1_fcn(t) - percent_change*delta_1_fcn(t)),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound delta_1 value
[t,x_UB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,(delta_1_fcn(t) + percent_change*delta_1_fcn(t)),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_table_1);
i = i + 1; % Incrementing iterator

% delta_2

% Solving system with lower bound delta_2 value
[t,x_LB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),(delta_2_fcn(t) - percent_change*delta_2_fcn(t)),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound delta_2 value
[t,x_UB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),(delta_2_fcn(t) + percent_change*delta_2_fcn(t)),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_table_1);
i = i + 1; % Incrementing iterator

% delta_3

% Solving system with lower bound delta_3 value
[t,x_LB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),(delta_3_fcn(t) - percent_change*delta_3_fcn(t)),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound delta_3 value
[t,x_UB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),(delta_3_fcn(t) + percent_change*delta_3_fcn(t)),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_table_1);
i = i + 1; % Incrementing iterator

% delta_4

% Solving system with lower bound delta_4 value
[t,x_LB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),(delta_4_fcn(t) - percent_change*delta_4_fcn(t)), ...
delta_5_fcn(t),p),tspan,x0,options);

% Solving system with upper bound delta_4 value
[t,x_UB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),(delta_4_fcn(t) + percent_change*delta_4_fcn(t)), ...
delta_5_fcn(t),p),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_table_1);
i = i + 1; % Incrementing iterator

% delta_5

% Solving system with lower bound delta_5 value
[t,x_LB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t),(delta_5_fcn(t) - percent_change*delta_5_fcn(t)), ...
p),tspan,x0,options);

% Solving system with upper bound delta_5 value
[t,x_UB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t),(delta_5_fcn(t) + percent_change*delta_5_fcn(t)), ...
p),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_table_1);
i = i + 1; % Incrementing iterator

% p

% Solving system with lower bound p value
[t,x_LB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t),delta_5_fcn(t), ...
(p - percent_change*p)),tspan,x0,options);

% Solving system with upper bound p value
[t,x_UB] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t),delta_5_fcn(t), ...
(p + percent_change*p)),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_table_1);

% Save percent_changes matrix as an Excel file
%xlswrite('sens_avg_pop_crit_p.xlsx',percent_changes)
 
% Function to update avg. population/days survived data
function [avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_table_1)

epsilon = 10^(-5);

% Storing average bee populations
avg_populations(i,1) = sum(x_LB(:,1)+x_LB(:,2)+x_LB(:,3))/length(tspan);
avg_populations(i,2) = sum(x_UB(:,1)+x_UB(:,2)+x_UB(:,3))/length(tspan);

% Storing percent changes in avg. populations
percent_changes(i,1) = 100*((avg_populations(i,1)/avg_pop_table_1) - 1);
percent_changes(i,2) = 100*((avg_populations(i,2)/avg_pop_table_1) - 1);
percent_changes(i,3) = abs(percent_changes(i,1)) + abs(percent_changes(i,2));

% Storing days survived by colony
if (find(x_LB(:,1) <= epsilon,1,'first'))
    days_survived(i,1) = find(x_LB(:,1) <= epsilon,1,'first');
else
    days_survived(i,1) = size(x_LB,1);
end

if (find(x_UB(:,1) <= epsilon,1,'first'))
    days_survived(i,2) = find(x_UB(:,1) <= epsilon,1,'first');
else
    days_survived(i,2) = size(x_UB,1);
end
% Total days survived per parameter
days_survived(i,3) = days_survived(i,1) + days_survived(i,2);

