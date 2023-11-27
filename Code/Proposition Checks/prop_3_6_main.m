clear all, close all, clc

% Parameter values (Spring)
beta_1 = 0.1984;
beta_2 = 0.1593;
beta_3 = 0.04959;
d_1 = 0.02272;
d_2 = 0.02272;
d_3 = 0.2;
mu = 500;
little_k = 0.000075;
r = 0.0165;
alpha = 0.4784;
K = 8000;
sigma_1 = 0.25;
sigma_2 = 0.75;

% % Parameter values (Fall)
% beta_1 = 0.1900;
% beta_2 = 0.1489;
% beta_3 = 0.04750;
% d_1 = 0.02272;
% d_2 = 0.02272;
% d_3 = 0.2;
% mu = 500;
% little_k = 0.000075;
% r = 0.0045;
% alpha = 0.5;
% K = 8000;
% sigma_1 = 0.25;
% sigma_2 = 0.75;

% Rate mites kill bees
gamma_1 = 10^(-7);
gamma_2 = gamma_1;
gamma_3 = 0.9; % gamma_3 > gamma_1,2

delta_i_zero = zeros(5,1); % Zero vector for no varroacide treatment
delta_i = [0.005; 0.005; 0.005; 0.5; 1;]; % Delta values

p = 0; % Initialzing p-value
i = 2; % i-value greater than 1

% Initial Conditions
x0 = [12000; 5581; 0; 0; 0;];

% Time span
tspan = 0:1:1000;

% Solving system with ode15s
[t,x] = ode15s(@(t,x) systemEqs(t,x,beta_1,beta_2,beta_3,d_1,d_2,d_3,mu,little_k,r, ...
    alpha,K,sigma_1,sigma_2,gamma_1,gamma_2,gamma_3,delta_i,p,i),tspan,x0);

% Plotting system solutions
plot_solutions(t,x)

% Checking stability based on propositions
check_prop(t,x,beta_3,d_1,d_2,mu,r,K,sigma_1,sigma_2,i,delta_i,p)


bool_vec = eqn_13(t,x,d_2,sigma_1,sigma_2,p);
%% No varroacide treatment (should be unstable)
% Solving system with ode15s
[t,x] = ode15s(@(t,x) systemEqs(t,x,beta_1,beta_2,beta_3,d_1,d_2,d_3,mu,little_k,r, ...
    alpha,K,sigma_1,sigma_2,gamma_1,gamma_2,gamma_3,delta_i_zero,p,i),tspan,x0);
% Plotting system solutions
plot_solutions(t,x)
% Checking stability based on propositions
check_prop(t,x,beta_3,d_1,d_2,mu,r,K,sigma_1,sigma_2,i,delta_i_zero,p)


%% Increased r-value (should be unstable)

r = 0.6;

% Solving system with ode15s
[t,x] = ode15s(@(t,x) systemEqs(t,x,beta_1,beta_2,beta_3,d_1,d_2,d_3,mu,little_k,r, ...
    alpha,K,sigma_1,sigma_2,gamma_1,gamma_2,gamma_3,delta_i,p,i),tspan,x0);
% Plotting system solutions
plot_solutions(t,x)
% Checking stability based on propositions
check_prop(t,x,beta_3,d_1,d_2,mu,r,K,sigma_1,sigma_2,i,delta_i,p)

%% Increased p-value

p = 0.6;

% Solving system with ode15s
[t,x] = ode15s(@(t,x) systemEqs(t,x,beta_1,beta_2,beta_3,d_1,d_2,d_3,mu,little_k,r, ...
    alpha,K,sigma_1,sigma_2,gamma_1,gamma_2,gamma_3,delta_i,p,i),tspan,x0);
% Plotting system solutions
plot_solutions(t,x)
% Checking stability based on propositions
check_prop(t,x,beta_3,d_1,d_2,mu,r,K,sigma_1,sigma_2,i,delta_i,p)

%% Increased mu-value

mu = 1000;

% Solving system with ode15s
[t,x] = ode15s(@(t,x) systemEqs(t,x,beta_1,beta_2,beta_3,d_1,d_2,d_3,mu,little_k,r, ...
    alpha,K,sigma_1,sigma_2,gamma_1,gamma_2,gamma_3,delta_i,p,i),tspan,x0);
% Plotting system solutions
plot_solutions(t,x)
% Checking stability based on propositions
check_prop(t,x,beta_3,d_1,d_2,mu,r,K,sigma_1,sigma_2,i,delta_i,p)

%% Decreased mu-value

mu = 200;

% Solving system with ode15s
[t,x] = ode15s(@(t,x) systemEqs(t,x,beta_1,beta_2,beta_3,d_1,d_2,d_3,mu,little_k,r, ...
    alpha,K,sigma_1,sigma_2,gamma_1,gamma_2,gamma_3,delta_i,p,i),tspan,x0);
% Plotting system solutions
plot_solutions(t,x)
% Checking stability based on propositions
check_prop(t,x,beta_3,d_1,d_2,mu,r,K,sigma_1,sigma_2,i,delta_i,p)

%% Increased K-value

K = 12000;

% Solving system with ode15s
[t,x] = ode15s(@(t,x) systemEqs(t,x,beta_1,beta_2,beta_3,d_1,d_2,d_3,mu,little_k,r, ...
    alpha,K,sigma_1,sigma_2,gamma_1,gamma_2,gamma_3,delta_i,p,i),tspan,x0);
% Plotting system solutions
plot_solutions(t,x)
% Checking stability based on propositions
check_prop(t,x,beta_3,d_1,d_2,mu,r,K,sigma_1,sigma_2,i,delta_i,p)

%%
function plot_pops = plot_solutions(t,x)
    hold on
    plot(t,x(:,1),'b',LineWidth=1);
    plot(t,x(:,2),'k--',LineWidth=1);
    plot(t,x(:,3),'-og','MarkerSize',2,'MarkerFaceColor','g',LineWidth=1)
    plot(t,x(:,4),'-.','Color',[0.4940 0.1840 0.556],LineWidth=1);
    plot(t,x(:,5),'r:',LineWidth=1);
    xlabel('Time (days)')
    ylabel('Population')
    legend('x_h','x_f','y','m','n');
    grid on
end

function bool_vec = eqn_13(t,x,d_2,sigma_1,sigma_2,p)

epsilon = 10^(-3);

% Equation 10 in 2017 paper
F_term = (sigma_1 - sigma_2 - p - d_2)/(p + d_2);
F = (F_term + sqrt(F_term^2 + 4*sigma_1/(p + d_2))) / 2;

% Checking equation 13 in for loop
for i = 1:length(t)
    if (x(i,2) == F*x(i,1) + epsilon) || (x(i,2) == F*x(i,1) - epsilon)
        bool_vec(i,1) = true;
    else
        bool_vec(i,1) = false;
    end
end
end


