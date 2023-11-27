clear all, close all, clc

% Attempting local sensitivity analysis on default model for Figure 3a
% (bee-only model) with slightly modified d_1 values

% Initial Conditions
x0 = [12500; 5000; 0; 0; 0;];

% Time span (days)
tspan = 0:1:3650; % 10 years

% Zero column vector for varroacide treatment
zero_vec = zeros(5,1);

% Increasing tolerances for ode45
options = odeset('RelTol',1e-3,'AbsTol',1e-5);

% Solving system with default (Table 1) parameter values 
[t,x] = ode45(@(t,x) systemEqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_i_fcn(t),0,zero_vec),tspan,x0,options);

% Solving system with modified d_1 values
[t1,x1] = ode45(@(t,x) systemEqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),(d_1_fcn(t) - 0.01), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_i_fcn(t),0,zero_vec),tspan,x0,options);

[t2,x2] = ode45(@(t,x) systemEqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),(d_1_fcn(t) + 0.01), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_i_fcn(t),0,zero_vec),tspan,x0,options);

% Calculating total avg population for 3 sims
total_avg_pop = (x(:,1) + x(:,2) + x(:,3))/length(tspan);
total_avg_pop_decreased = (x1(:,1) + x1(:,2) + x1(:,3))/length(tspan);
total_avg_pop_increased = (x2(:,1) + x2(:,2) + x2(:,3))/length(tspan);

% Attempting finite difference to calculate sensitivity: (default - modified)/(change in parameter)
sens_decreased = (total_avg_pop_decreased - total_avg_pop)/0.1;
%sens_increased = (total_avg_pop_increased - total_avg_pop)/0.1;

% % Plotting system solutions
% figure(1);
% hold on
% plot(t,x1(:,1),'b',LineWidth=1);
% plot(t,x1(:,2),'k--',LineWidth=1);
% plot(t,x1(:,3),'-og','MarkerSize',2,'MarkerFaceColor','g',LineWidth=1);
% plot(t,x1(:,4),'-.','Color',[0.4940 0.1840 0.556],LineWidth=1);
% plot(t,x1(:,5),'r:',LineWidth=1);
% xlabel('Time (days)')
% ylabel('Population')
% ylim([0 30000])
% xlim([0 3650])
% legend('x_h','x_f','y','m','n');
% grid on

% Function to input model equations
function dx = systemEqs(t,x,beta_1,beta_2,beta_3,d_1,d_2,d_3,mu,k,r,alpha, ...
    K,sigma_1,sigma_2,gamma_i,p,delta_i)

% Epsilon (avoid div. by zero)
epsilon = 10^(-5);

healthy_bees = x(1) + x(2);

i = 2; % i-value greater than 1
% Function g (brood maintenance term w/ above i-value)
g = (healthy_bees^i)/(K^i + healthy_bees^i + epsilon);

% Function h (affects eclosion rate based on virus-carrying mites x(4))
h = exp(-x(4)*k);

% Function R
R = sigma_1 - sigma_2*(x(2)/(x(1) + x(2) + epsilon));

% Sum of all bees/mites
bees_sum = x(1) + x(2) + x(3);
mites_sum = x(4) + x(5);

% Proportions of bees
prop_hive = x(1)/(bees_sum + epsilon);
prop_forager = x(2)/(bees_sum + epsilon);
prop_healthy = healthy_bees/(bees_sum + epsilon);
prop_infected = x(3)/(bees_sum + epsilon);

% Mite logistic growth term
mite_logistic = 1 - (mites_sum/(alpha*bees_sum + epsilon));

% Setting period of exposure 
temp_t = t;
while (temp_t > 365)
    temp_t = temp_t - 365;
end
if ((temp_t < 5) || (temp_t > 35)) 
    p = 0; % Overriding set p-value
end

% 5-D complete system (eqs. 1-5)
dx = [
    mu*g*h - beta_1*x(4)*prop_hive - (d_1 + delta_i(1))*x(1) - gamma_i(1)*mites_sum*x(1) - x(1)*R;
    x(1)*R - beta_1*x(4)*prop_forager - (p + d_2 + delta_i(2))*x(2) - gamma_i(2)*mites_sum*x(2);
    beta_1*x(4)*prop_healthy - (d_3 + delta_i(3))*x(3) - gamma_i(3)*mites_sum*x(3);
    r*x(4)*mite_logistic + beta_2*x(5)*prop_infected - beta_3*x(4)*prop_healthy - delta_i(4)*x(4);
    r*x(5)*mite_logistic - beta_2*x(5)*prop_infected + beta_3*x(4)*prop_healthy - delta_i(5)*x(5);
    ];
end