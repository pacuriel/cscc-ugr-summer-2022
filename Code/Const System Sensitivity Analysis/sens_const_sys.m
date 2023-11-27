% Function sensitivity analysis on model w/ constant coefficients
function [avg_populations,percent_changes,days_survived] = sens_const_sys(x0,tspan,percent_change)

% Parameter values (Spring)
% beta_1 = 0.1984;
% beta_2 = 0.1593;
% beta_3 = 0.04959;
% d_1 = 0.02272;
% d_2 = 0.02272;
% d_3 = 0.2;
% mu = 500;
% k = 0.000075;
% r = 0.0165;
% alpha = 0.4784;
% K = 8000;
% sigma_1 = 0.25;
% sigma_2 = 0.75;

% Parameter values (Fall)
beta_1 = 0.1900;
beta_2 = 0.1489;
beta_3 = 0.04750;
d_1 = 0.02272;
d_2 = 0.02272;
d_3 = 0.2;
mu = 500;
little_k = 0.000075;
r = 0.0045;
alpha = 0.5;
K = 8000;
sigma_1 = 0.25;
sigma_2 = 0.75;

% Rate mites kill bees
gamma_1 = 10^(-7);
gamma_2 = gamma_1;
gamma_3 = 0.0000002; % gamma_3 > gamma_1,2

% Initializing zero matrices to store LB/UB avgs. and percent changes
% Col1 = LB avgs./percentages, Col2 = UB avgs./percentages, Col3 = Total change, Rows = parameters
% LB = lower bound, UB = upper bound
avg_populations = zeros(16,3);
percent_changes = zeros(16,3);
days_survived = zeros(16,3);

% Names of parameters in order
param_names = ["beta_1"; "beta_2"; "beta_3"; "d_1"; "d_2"; "d_3"; "mu"; 
    "little_k"; "r"; "alpha"; "K"; "sigma_1"; "sigma_2"; "gamma_1"; "gamma_2";
    "gamma_3"; "delta_1"; "delta_2"; "delta_3"; "delta_4"; ""; "p";];

i = 1; % Initializing iterator value
epsilon = 10^(-5); % Epsilon to track colony failure

% Increasing tolerances for ode15s
options = odeset('RelTol',1e-7,'AbsTol',1e-7,'NonNegative',1:5);

% Solving system with default (Table 1) parameter values 
[t,x] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,d_1, ...
d_2,d_3,mu,little_k,r,alpha,K,sigma_1, ...
sigma_2,gamma_1,gamma_2,gamma_3 ...
),tspan,x0,options);

% Storing days survived by colony
if (find(x(:,1) <= epsilon,1,'first'))
    days_survived_table_1(i,1) = find(x(:,1) <= epsilon,1,'first');
else
    days_survived_table_1(i,1) = size(x,1);
end

% Calculating avgerage total bee population
avg_pop_Fall = sum(x(:,1)+x(:,2)+x(:,3))/length(tspan)

%% Running simulation with +/- 10% for each parameter

% beta_1

% Solving system with lower bound beta_1 value
[t,x_LB] = ode15s(@(t,x) const_sys_eqs(t,x,(beta_1 - percent_change*beta_1),beta_2,beta_3,d_1, ...
d_2,d_3,mu,little_k,r,alpha,K,sigma_1, ...
sigma_2,gamma_1,gamma_2,gamma_3 ...
),tspan,x0,options);

% Solving system with upper bound beta_1 value
[t,x_UB] = ode15s(@(t,x) const_sys_eqs(t,x,(beta_1 + percent_change*beta_1),beta_2,beta_3,d_1, ...
d_2,d_3,mu,little_k,r,alpha,K,sigma_1, ...
sigma_2,gamma_1,gamma_2,gamma_3 ...
),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_Fall);
i = i + 1; % Incrementing iterator

% beta_2

% Solving system with lower bound beta_2 value
[t,x_LB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,(beta_2 - percent_change*beta_2),beta_3,d_1, ...
d_2,d_3,mu,little_k,r,alpha,K,sigma_1, ...
sigma_2,gamma_1,gamma_2,gamma_3 ...
),tspan,x0,options);

% Solving system with upper bound beta_2 value
[t,x_UB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,(beta_2 + percent_change*beta_2),beta_3,d_1, ...
d_2,d_3,mu,little_k,r,alpha,K,sigma_1, ...
sigma_2,gamma_1,gamma_2,gamma_3 ...
),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_Fall);
i = i + 1; % Incrementing iterator

% beta_3

% Solving system with lower bound beta_3 value
[t,x_LB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,(beta_3 - percent_change*beta_3),d_1, ...
d_2,d_3,mu,little_k,r,alpha,K,sigma_1, ...
sigma_2,gamma_1,gamma_2,gamma_3 ...
),tspan,x0,options);

% Solving system with upper bound beta_3 value
[t,x_UB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,(beta_3 + percent_change*beta_3),d_1, ...
d_2,d_3,mu,little_k,r,alpha,K,sigma_1, ...
sigma_2,gamma_1,gamma_2,gamma_3 ...
),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_Fall);
i = i + 1; % Incrementing iterator

% d_1

% Solving system with lower bound d_1 value
[t,x_LB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,(d_1 - percent_change*d_1), ...
d_2,d_3,mu,little_k,r,alpha,K,sigma_1, ...
sigma_2,gamma_1,gamma_2,gamma_3 ...
),tspan,x0,options);

% Solving system with upper bound d_1 value
[t,x_UB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,(d_1 + percent_change*d_1), ...
d_2,d_3,mu,little_k,r,alpha,K,sigma_1, ...
sigma_2,gamma_1,gamma_2,gamma_3 ...
),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_Fall);
i = i + 1; % Incrementing iterator

% d_2

% Solving system with lower bound d_2 value
[t,x_LB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,d_1, ...
(d_2 - percent_change*d_2),d_3,mu,little_k,r,alpha,K,sigma_1, ...
sigma_2,gamma_1,gamma_2,gamma_3 ...
),tspan,x0,options);

% Solving system with upper bound d_2 value
[t,x_UB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,d_1, ...
(d_2 + percent_change*d_2),d_3,mu,little_k,r,alpha,K,sigma_1, ...
sigma_2,gamma_1,gamma_2,gamma_3 ...
),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_Fall);
i = i + 1; % Incrementing iterator

% d_3

% Solving system with lower bound d_3 value
[t,x_LB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,d_1, ...
d_2,(d_3 - percent_change*d_3),mu,little_k,r,alpha,K,sigma_1, ...
sigma_2,gamma_1,gamma_2,gamma_3 ...
),tspan,x0,options);

% Solving system with upper bound d_3 value
[t,x_UB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,d_1, ...
d_2,(d_3 + percent_change*d_3),mu,little_k,r,alpha,K,sigma_1, ...
sigma_2,gamma_1,gamma_2,gamma_3 ...
),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_Fall);
i = i + 1; % Incrementing iterator

% mu

% Solving system with lower bound mu value
[t,x_LB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,d_1, ...
d_2,d_3,(mu - percent_change*mu),little_k,r,alpha,K,sigma_1, ...
sigma_2,gamma_1,gamma_2,gamma_3 ...
),tspan,x0,options);

% Solving system with upper bound mu value
[t,x_UB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,d_1, ...
d_2,d_3,(mu + percent_change*mu),little_k,r,alpha,K,sigma_1, ...
sigma_2,gamma_1,gamma_2,gamma_3 ...
),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_Fall);
i = i + 1; % Incrementing iterator

% little_k

% Solving system with lower bound little_k value
[t,x_LB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,d_1, ...
d_2,d_3,mu,(little_k - percent_change*little_k),r,alpha,K,sigma_1, ...
sigma_2,gamma_1,gamma_2,gamma_3 ...
),tspan,x0,options);

% Solving system with upper bound little_k value
[t,x_UB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,d_1, ...
d_2,d_3,mu,(little_k + percent_change*little_k),r,alpha,K,sigma_1, ...
sigma_2,gamma_1,gamma_2,gamma_3 ...
),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_Fall);
i = i + 1; % Incrementing iterator

% r

% Solving system with lower bound r value
[t,x_LB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,d_1, ...
d_2,d_3,mu,little_k,(r - percent_change*r),alpha,K,sigma_1, ...
sigma_2,gamma_1,gamma_2,gamma_3 ...
),tspan,x0,options);

% Solving system with upper bound r value
[t,x_UB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,d_1, ...
d_2,d_3,mu,little_k,(r + percent_change*r),alpha,K,sigma_1, ...
sigma_2,gamma_1,gamma_2,gamma_3 ...
),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_Fall);
i = i + 1; % Incrementing iterator

% alpha

% Solving system with lower bound alpha value
[t,x_LB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,d_1, ...
d_2,d_3,mu,little_k,r,(alpha - percent_change*alpha),K,sigma_1, ...
sigma_2,gamma_1,gamma_2,gamma_3 ...
),tspan,x0,options);

% Solving system with upper bound alpha value
[t,x_UB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,d_1, ...
d_2,d_3,mu,little_k,r,(alpha + percent_change*alpha),K,sigma_1, ...
sigma_2,gamma_1,gamma_2,gamma_3 ...
),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_Fall);
i = i + 1; % Incrementing iterator

% K

% Solving system with lower bound K value
[t,x_LB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,d_1, ...
d_2,d_3,mu,little_k,r,alpha,(K - percent_change*K),sigma_1, ...
sigma_2,gamma_1,gamma_2,gamma_3 ...
),tspan,x0,options);

% Solving system with upper bound K value
[t,x_UB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,d_1, ...
d_2,d_3,mu,little_k,r,alpha,(K + percent_change*K),sigma_1, ...
sigma_2,gamma_1,gamma_2,gamma_3 ...
),tspan,x0,options); %#ok<*ASGLU> 

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_Fall);
i = i + 1; % Incrementing iterator

% sigma_1

% Solving system with lower bound sigma_1 value
[t,x_LB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,d_1, ...
d_2,d_3,mu,little_k,r,alpha,K,(sigma_1 - percent_change*sigma_1), ...
sigma_2,gamma_1,gamma_2,gamma_3 ...
),tspan,x0,options);

% Solving system with upper bound sigma_1 value
[t,x_UB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,d_1, ...
d_2,d_3,mu,little_k,r,alpha,K,(sigma_1 + percent_change*sigma_1), ...
sigma_2,gamma_1,gamma_2,gamma_3 ...
),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_Fall);
i = i + 1; % Incrementing iterator

% sigma_2

% Solving system with lower bound sigma_2 value
[t,x_LB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,d_1, ...
d_2,d_3,mu,little_k,r,alpha,K,sigma_1, ...
(sigma_2 - percent_change*sigma_2),gamma_1,gamma_2,gamma_3 ...
),tspan,x0,options);

% Solving system with upper bound sigma_2 value
[t,x_UB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,d_1, ...
d_2,d_3,mu,little_k,r,alpha,K,sigma_1, ...
(sigma_2 + percent_change*sigma_2),gamma_1,gamma_2,gamma_3),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_Fall);
i = i + 1; % Incrementing iterator

% gamma_1

% Solving system with lower bound gamma_1 value
[t,x_LB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,d_1, ...
d_2,d_3,mu,little_k,r,alpha,K,sigma_1, ...
sigma_2,(gamma_1 - percent_change*gamma_1),gamma_2,gamma_3 ...
),tspan,x0,options);

% Solving system with upper bound gamma_1 value
[t,x_UB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,d_1,d_2,d_3,mu,little_k,r,alpha,K,sigma_1, ...
sigma_2,(gamma_1 + percent_change*gamma_1),gamma_2,gamma_3),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_Fall);
i = i + 1; % Incrementing iterator

% gamma_2

% Solving system with lower bound gamma_2 value
[t,x_LB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,d_1, ...
d_2,d_3,mu,little_k,r,alpha,K,sigma_1, ...
sigma_2,gamma_1,(gamma_2 - percent_change*gamma_2),gamma_3 ...
),tspan,x0,options);

% Solving system with upper bound gamma_2 value
[t,x_UB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,d_1, ...
d_2,d_3,mu,little_k,r,alpha,K,sigma_1, ...
sigma_2,gamma_1,(gamma_2 + percent_change*gamma_2),gamma_3 ...
),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_Fall);
i = i + 1; % Incrementing iterator

% gamma_3

% Solving system with lower bound gamma_3 value
[t,x_LB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,d_1, ...
d_2,d_3,mu,little_k,r,alpha,K,sigma_1, ...
sigma_2,gamma_1,gamma_2,(gamma_3 - percent_change*gamma_3) ...
),tspan,x0,options);

% Solving system with upper bound gamma_3 value
[t,x_UB] = ode15s(@(t,x) const_sys_eqs(t,x,beta_1,beta_2,beta_3,d_1, ...
d_2,d_3,mu,little_k,r,alpha,K,sigma_1, ...
sigma_2,gamma_1,gamma_2,(gamma_3 + percent_change*gamma_3) ...
),tspan,x0,options);

[avg_populations,percent_changes,days_survived] = calc_sens(x_LB,x_UB,avg_populations,percent_changes,days_survived,i,tspan,avg_pop_Fall);

% Save percent_changes matrix as an Excel file
xlswrite('const_Fall_sens.xlsx',percent_changes)

% Function to update avg. population data
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
