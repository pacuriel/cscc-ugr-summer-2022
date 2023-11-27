clear all, close all, clc
tspan = 0:1:365; % Time span

%%
figure(1); % beta_i
hold on
plot(tspan,beta_1(tspan));
plot(tspan,beta_2(tspan));
plot(tspan,beta_3(tspan));
legend('\beta_1','\beta_2','\beta_3');
ylim([0 0.2])
xlim([0 730])
ylabel('\beta_i(t)')
xlabel('Time t (days)')
grid on

%%
figure(2); % d_i
hold on
plot(tspan,d_1_fcn(tspan),'LineWidth',4);
% plot(tspan,d_2_(tspan));
% plot(tspan,d_3_(tspan));
legend('d_1','d_2','d_3');
ylim([0 0.2])
xlim([0 365])
ylabel('d_i(t)')
xlabel('Time t (days)')
grid on

%%
figure(3); % mu
plot(tspan,mu(tspan),'LineWidth',4);
ylim([0 2000])
xlim([0 730])
ylabel('\mu(t)')
xlabel('Time t (days)')
grid on

%%
figure(4); % k
plot(tspan,little_k(tspan),'LineWidth',4);
xlim([0 730])
ylabel('k(t)')
xlabel('Time t (days)')
grid on

%%
figure(5); % r
plot(tspan,r(tspan),'LineWidth',4);
ylim([0 0.02])
xlim([0 730])
ylabel('r(t)')
xlabel('Time t (days)')
grid on

%%
figure(6); % alpha
plot(tspan,alpha(tspan),'LineWidth',4);
ylim([0 0.6])
xlim([0 730])
ylabel('\alpha(t)')
xlabel('Time t (days)')
grid on

%%
figure(7); % K
plot(tspan,K(tspan),'LineWidth',4);
ylim([0 12000])
xlim([0 730])
ylabel('K(t)')
xlabel('Time t (days)')
grid on

%% 
figure(8); % sigma_i
hold on
plot(tspan,sigma_1_fcn(tspan));
plot(tspan,sigma_2_fcn(tspan));
legend('\sigma_1','\sigma_2');
ylim([0 1])
xlim([0 730])
ylabel('\sigma_i(t)')
xlabel('Time t (days)')
grid on

%%
figure(9); % Plotting all/most parameters 
hold on
plot(tspan,beta_1(tspan));
plot(tspan,beta_2(tspan));
plot(tspan,beta_3(tspan));
plot(tspan,d_1(tspan));
plot(tspan,d_2(tspan));
plot(tspan,d_3(tspan));
plot(tspan,little_k(tspan));
plot(tspan,r(tspan));
plot(tspan,alpha(tspan));
plot(tspan,sigma_1(tspan));
plot(tspan,sigma_2(tspan));
legend('\beta_1','\beta_2','\beta_3','d_1','d_2','d_3','k','r', ...
    '\alpha','\sigma_1','\sigma_2');