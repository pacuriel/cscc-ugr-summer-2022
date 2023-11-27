clear all, close all, clc

% Parameter values (Spring)
beta_1 = 0.1984;
beta_2 = 0.1593;
beta_3 = 0.04959;
d_1 = 0.02272;
d_2 = 0.02272;
d_3 = 0.2;
mu = 500;
k = 0.000075;
r = 0.0165;
alpha = 0.4784;
K = 8000;
sigma_1 = 0.25;
sigma_2 = 0.75;

% Initial Conditions
x0 = [12000; 5000;];

% Time span
tspan = 0:0.1:100;

% Solving system with ode45
[t,x] = ode45(@(t,x) systemEqs(t,x,d_1,d_2,mu,K,sigma_1,sigma_2),tspan,x0);

% Plotting system solutions
figure(1);
plot(t,x(:,1),'b',t,x(:,2),'r');
xlabel('Time t')
ylabel('Bee populations')
legend('x_h','x_f')
grid on

% % Solution curve (x_h-x_f)
% figure(2);
% plot(x(:,1),x(:,2),'k');
% xlabel('x_h')
% ylabel('x_f')
% legend('Solution curve')
% grid on

% Function to input model equations
function dx = systemEqs(t,x,d_1,d_2,mu,K,sigma_1,sigma_2)
i = 2; % i-value greater than 1
% Function g (w/ above i-value)
g = ((x(1) + x(2))^i)/(K^i + (x(1) + x(2))^i);

% Function R
R = sigma_1 - sigma_2*(x(2)/(x(1) + x(2)));

% Setting period of exposure 
p = 0; % Initialzing p-value
if ((t < 5) || (t > 35)) 
    p = 0;
end

% 2-D disease-free model (eqs. 8 and 9)
dx = [
    mu*g - d_1*x(1) - x(1)*R;
    x(1)*R - (p + d_2)*x(2);
    ];
end