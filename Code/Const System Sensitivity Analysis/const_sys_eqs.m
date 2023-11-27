% Function to input model equations
function dx = const_sys_eqs(t,x,beta_1,beta_2,beta_3,d_1,d_2,d_3,mu,k,r,alpha, ...
    K,sigma_1,sigma_2,gamma_1,gamma_2,gamma_3)
i = 2; % i-value greater than 1
% Function g (w/ above i-value)
g = ((x(1) + x(2))^i)/(K^i + (x(1) + x(2))^i);

% Function h (affects eclosion rate based on virus-carrying mites x(4))
h = exp(-x(4)*k);

% Epsilon (avoid div. by zero)
epsilon = 10^(-5);

% Function R
R = sigma_1 - sigma_2*(x(2)/(x(1) + x(2) + epsilon));

% Delta values (Not applying varroacide treatment)
delta_1 = 0;
delta_2 = delta_1;
delta_3 = delta_1;
delta_4 = delta_1;
delta_5 = delta_1;

healthy_bees = x(1) + x(2);

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
p = 0;
temp_t = t;
while (temp_t > 365)
    temp_t = temp_t - 365;
end
if ((temp_t < 5) || (temp_t > 35)) 
    p = 0;
end

% 5-D complete system (eqs. 1-5)
dx = [
    mu*g*h - beta_1*x(4)*prop_hive - (d_1 + delta_1)*x(1) - gamma_1*mites_sum*x(1) - x(1)*R;
    x(1)*R - beta_1*x(4)*prop_forager - (p + d_2 + delta_2)*x(2) - gamma_2*mites_sum*x(2);
    beta_1*x(4)*prop_healthy - (d_3 + delta_3)*x(3) - gamma_3*mites_sum*x(3);
    r*x(4)*mite_logistic + beta_2*x(5)*prop_infected - beta_3*x(4)*prop_healthy - delta_4*x(4);
    r*x(5)*mite_logistic - beta_2*x(5)*prop_infected + beta_3*x(4)*prop_healthy - delta_5*x(5);
    ];