clear all, close all, clc

% Sensitivity analysis on full model

% Initial Conditions
x0 = [12500; 5000; 0; 0; 0;];

% Time span (days)
tspan = 0:1:3650; % 10 years

% Zero column vector to not apply varroacide treatment
zero_vec = zeros(5,1);

% gamma_i values (Rate mites kills bees)
gamma_1 = 10^(-7);
gamma_2 = gamma_1;
gamma_3 = 0.0000002;

% Initialzing p-value (homing failure)
p = 0.1; 

% Increasing tolerances for ode15s
options = odeset('RelTol',1e-3,'AbsTol',1e-5,'NonNegative',1:5);

% Solving system with default (Table 1) parameter values 
[t,x] = ode15s(@(t,x) system_eqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,delta_1_fcn(t),delta_2_fcn(t),delta_3_fcn(t),delta_4_fcn(t), ...
delta_5_fcn(t),p),tspan,x0,options);

% Plotting system solutions
figure(1);
hold on
plot(t,x(:,1),'b',LineWidth=1);
plot(t,x(:,2),'k--',LineWidth=1);
plot(t,x(:,3),'-og','MarkerSize',2,'MarkerFaceColor','g',LineWidth=1)
plot(t,x(:,4),'-.','Color',[0.4940 0.1840 0.556],LineWidth=1);
plot(t,x(:,5),'r:',LineWidth=1);
xlabel('Time (days)')
ylabel('Population')
ylim([0 30000])
xlim([0 3000])
legend('x_h','x_f','y','m','n');
grid on