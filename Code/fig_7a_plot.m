clear all, close all, clc

% Rate mites eat bees
gamma_1 = 10^(-7);
gamma_2 = gamma_1;
gamma_3 = 0.0000002; % gamma_3 > gamma_1,2

% Initial Conditions
x0 = [12500; 5000; 0; 0; 0;];

% Time span (days)
tspan = 0:1:7300; % 20 years

% Zero column vector for varroacide treatment
zero_vec = zeros(5,1);

% Increasing tolerances for ode45
options = odeset('RelTol',1e-3,'AbsTol',1e-5);

% p-values to plot
p_vals = 0:0.01:0.2;

% Initializing matrix to store avg bee populations
avg_bee_pop = zeros(1,numel(p_vals));

i = 1; % Initialzing iterator
% while loop to get average bee populations for varying p-values
while (i <= numel(p_vals))
    
    % Solving system with ode45
    [t,x] = ode45(@(t,x) systemEqs(t,x,beta_1_fcn(t),beta_2_fcn(t),beta_3_fcn(t),d_1_fcn(t), ...
    d_2_fcn(t),d_3_fcn(t),mu_fcn(t),little_k_fcn(t),r_fcn(t),alpha_fcn(t),K_fcn(t),sigma_1_fcn(t), ...
    sigma_2_fcn(t),gamma_1,gamma_2,gamma_3,p_vals(i),zero_vec),tspan,x0,options);

    j = 1; % Initialzing iterator
    % while loop to get average based on last year
    while (j < numel(tspan))
        % Condition to terminate while loop
        if (numel(x(:,1)) ~= numel(tspan))
            break
        end

        if (tspan(j) == 6935) % End of 19th year
     
            % Getting sum/avg of all bees during last year of sim
            sum_all_bees = sum(x(j:numel(tspan),1)) + sum(x(j:numel(tspan),2)) + sum(x(j:numel(tspan),3));
            
            % Storing avg bees for set p-value
            avg_bee_pop(1,i) = sum_all_bees/365;
        end
        j = j + 1; % Incrementing iterator
    end

%     % if average is undefined, set to zero
%     if (isnan(avg_bee_pop(1,i)))
%         avg_bee_pop(1,i) = 0;
%     end

    i = i + 1 ; % Incrementing iterator
end

% Plotting Fig_7a
figure(1); 
plot(p_vals,avg_bee_pop,LineWidth=2);
ylim([0 30000])
xlim([0 0.2])
xlabel('Forager loss due to homing failure (p)')
ylabel('Average total bee population')
legend('Disease-free')
grid on 

% Plotting system solutions
% figure(1);
% hold on
% plot(t,x(:,1),'b',LineWidth=1);
% plot(t,x(:,2),'k--',LineWidth=1);
% plot(t,x(:,3),'-og','MarkerSize',2,'MarkerFaceColor','g',LineWidth=1);
% plot(t,x(:,4),'-.','Color',[0.4940 0.1840 0.556],LineWidth=1);
% plot(t,x(:,5),'r:',LineWidth=1);
% xlabel('Time (days)')
% ylabel('Population')
% ylim([0 30000])
% xlim([0 7300])
% legend('x_h','x_f','y','m','n');
% grid on

% Function to input model equations
function dx = systemEqs(t,x,beta_1,beta_2,beta_3,d_1,d_2,d_3,mu,k,r,alpha, ...
    K,sigma_1,sigma_2,gamma_1,gamma_2,gamma_3,p,delta_i)

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
    mu*g*h - beta_1*x(4)*prop_hive - (d_1 + delta_i(1))*x(1) - gamma_1*mites_sum*x(1) - x(1)*R;
    x(1)*R - beta_1*x(4)*prop_forager - (p + d_2 + delta_i(2))*x(2) - gamma_2*mites_sum*x(2);
    beta_1*x(4)*prop_healthy - (d_3 + delta_i(3))*x(3) - gamma_3*mites_sum*x(3);
    r*x(4)*mite_logistic + beta_2*x(5)*prop_infected - beta_3*x(4)*prop_healthy - delta_i(4)*x(4);
    r*x(5)*mite_logistic - beta_2*x(5)*prop_infected + beta_3*x(4)*prop_healthy - delta_i(5)*x(5);
    ];
end