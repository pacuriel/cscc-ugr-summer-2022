clear all, close all, clc

% Time span
tspan = 0:1:365;

% Obtaining parameter values
d_1 = d_1_fcn(tspan);
mu = mu_fcn(tspan);
K = K_fcn(tspan);
%%
% Plotting d_1
figure(1); 
hold on
plot(tspan,d_1,'LineWidth',2);
legend('d_1');
ylim([0 0.25])
xlim([0 730])
ylabel('d_1(t)')
xlabel('Time t (days)')
grid on
%%
% Plotting mu
figure(2); 
plot(tspan,mu,'LineWidth',3);
ylim([0 2000])
xlim([0 length(tspan)])
ylabel('\mu(t)')
xlabel('Time t (days)')
title('Max. eclosion rate \mu over 1 year')
grid on
%%
figure(3); % K
plot(tspan,K,'LineWidth',2);
ylim([0 13000])
xlim([0 730])
ylabel('K(t)')
xlabel('Time t (days)')
grid on
