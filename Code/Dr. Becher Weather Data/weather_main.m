clear all, close all, clc

% Initial Conditions
x0 = [12500; 5000; 0; 0; 0;];

% Time span (days)
tspan = 0:1:365; % 10 years

% Function to assign values to parameter vectors
[beta_1,beta_2,beta_3,d_1,d_2,d_3,mu,little_k,r,alpha,K,sigma_1,sigma_2,...
    gamma_1,gamma_2,gamma_3,delta_1,delta_2,delta_3,delta_4,delta_5,p] = parameter_avgs(tspan); 

% Columns = (tmax,tmin,sun,relh,t_avg,julian)
weather_data = readmatrix('Rothamstead_data_cleaned.csv');
%%

plot(tspan,weather_data(1:366,5),LineWidth=2);
xlim([0 length(tspan)]);
xlabel('Time (days)')
ylabel('Avg. Temp. (C)')
title('Rothamsted, UK (1995)');
grid on




%%
% Modifying p values based on weather data
p = mod_p(tspan,p,weather_data);

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