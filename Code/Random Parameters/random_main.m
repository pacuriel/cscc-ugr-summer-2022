clear all, close all, clc

% Initial Conditions
x0 = [12500; 5000; 0; 10; 50;];

% Time span (days)
tspan = 0:1:3650; % 10 years

% Obtaining random parameter values
parameter_vals = parameter_avgs(tspan);

% Assigning parameters their own column vector
beta_1 = parameter_vals(:,1);
beta_2 = parameter_vals(:,2);
beta_3 = parameter_vals(:,3);
d_1 = parameter_vals(:,4);
d_2 = parameter_vals(:,5);
d_3 = parameter_vals(:,6);
mu = parameter_vals(:,7);
little_k = parameter_vals(:,8);
r = parameter_vals(:,9);
alpha = parameter_vals(:,10);
K = parameter_vals(:,11);
sigma_1 = parameter_vals(:,12);
sigma_2 = parameter_vals(:,13);
gamma_1 = parameter_vals(:,14);
gamma_2 = parameter_vals(:,15);
gamma_3 = parameter_vals(:,16);
delta_1 = parameter_vals(:,17);
delta_2 = parameter_vals(:,18);
delta_3 = parameter_vals(:,19);
delta_4 = parameter_vals(:,20);
delta_5 = parameter_vals(:,21);
p = parameter_vals(:,22);

% Increasing tolerances for ode15s
options = odeset('RelTol',1e-5,'AbsTol',1e-7,'NonNegative',1:5);

% Solving system with default (Table 1) parameter values 
[t,x] = ode15s(@(t,x) system_eqs(t,x,beta_1,beta_2,beta_3,d_1,d_2,d_3,mu, ...
    little_k,r,alpha,K,sigma_1,sigma_2,gamma_1,gamma_2,gamma_3,delta_1, ...
    delta_2,delta_3,delta_4,delta_5,p),tspan,x0,options);

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
xlim([0 length(tspan)])
legend('x_h','x_f','y','m','n');
grid on