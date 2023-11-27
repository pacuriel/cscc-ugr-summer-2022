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
p = 0.13;

% Rate mites kill bees
gamma_1 = 10^(-7);
gamma_2 = gamma_1;
gamma_3 = 0.9; % gamma_3 > gamma_1,2

% Initial Conditions
x0 = [13000; 5558; 0; 0; 0;];

x0_2 = [12500; 5000; 0; 10; 500;];

% Time span
tspan = 0:1:100;

% Solving disease-free system
[t1,x1] = ode45(@(t,x) systemEqs(t,x,beta_1,beta_2,beta_3,d_1,d_2,d_3,mu,k,r, ...
    alpha, K,sigma_1,sigma_2,gamma_1,gamma_2,gamma_3,p),tspan,x0);

% Solving disease-infested system
[t2,x2] = ode45(@(t,x) systemEqs(t,x,beta_1,beta_2,beta_3,d_1,d_2,d_3,mu,k,r, ...
    alpha, K,sigma_1,sigma_2,gamma_1,gamma_2,gamma_3,p),tspan,x0_2);

% Trying to plot Fig. 2
figure(2);
hold on
% Shaded rectangle denotes period of exposure
rectangle('Position',[5 0 30 25000],'FaceColor',[.9 .9 .9],'EdgeColor',[.9 .9 .9])
plot(t1,(x1(:,1) + x1(:,2) + x1(:,3)),'r');
plot(t2,(x2(:,1) + x2(:,2) + x2(:,3)),'b--');
xlabel('Time (in days)')
ylabel('Total bee population size')
ylim([0 25000])
legend('Disease-free', 'Disease-infested')
grid on

% Function to input model equations
function dx = systemEqs(t,x,beta_1,beta_2,beta_3,d_1,d_2,d_3,mu,k,r,alpha, ...
    K,sigma_1,sigma_2,gamma_1,gamma_2,gamma_3,p)
i = 2; % i-value greater than 1
% Function g (w/ above i-value)
g = ((x(1) + x(2))^i)/(K^i + (x(1) + x(2))^i);

% Function h (affects eclosion rate based on virus-carrying mites x(4))
h = exp(-x(4)*k);

% Function R
R = sigma_1 - sigma_2*(x(2)/(x(1) + x(2)));

% Delta values (Not applying varroacide treatment)
delta_1 = 0;
delta_2 = delta_1;
delta_3 = delta_1;
delta_4 = delta_1;
delta_5 = delta_1;

% Setting period of exposure 
if ((t < 5) || (t > 35)) 
    p = 0;
end

% 5-D complete system (eqs. 1-5)
dx = [
    mu*g*h - beta_1*x(4)*(x(1)/(x(1) + x(2) + x(3))) - (d_1 + delta_1)*x(1) - gamma_1*(x(4) + x(5))*x(1) - x(1)*R;
    x(1)*R - beta_1*x(4)*(x(2)/(x(1) + x(2) + x(3))) - (p + d_2 + delta_2)*x(2) - gamma_2*(x(4) + x(5))*x(2);
    beta_1*x(4)*(x(1) + x(2))/(x(1) + x(2) + x(3)) - (d_3 + delta_3)*x(3) - gamma_3*(x(4) + x(5))*x(3);
    r*x(4)*(1 - ((x(4) + x(5))/(alpha*(x(1) + x(2) + x(3))))) + beta_2*x(5)*(x(3)/(x(1) + x(2) + x(3))) - beta_3*x(4)*(x(1) + x(2))/(x(1) + x(2) + x(3)) - delta_4*x(4);
    r*x(5)*(1 - ((x(4) + x(5))/(alpha*(x(1) + x(2) + x(3))))) - beta_2*x(5)*(x(3)/(x(1) + x(2) + x(3))) + beta_3*x(4)*(x(1) + x(2))/(x(1) + x(2) + x(3)) - delta_5*x(5);
    ];
end