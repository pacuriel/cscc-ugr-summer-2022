% Function to input model equations
function dx = system_eqs(t,x,beta_1,beta_2,beta_3,d_1,d_2,d_3,mu,little_k,r,alpha, ...
    K,sigma_1,sigma_2,gamma_1,gamma_2,gamma_3,delta_1,delta_2,delta_3, ...
    delta_4,delta_5,p)

% Nearest day-value
day = round(t) + 1;

% Epsilon (avoid div. by zero)
epsilon = 10^(-5);

% Storing total healthy bees
healthy_bees = x(1) + x(2);

i = 2; % i-value greater than 1
% Function g (brood maintenance term w/ above i-value)
g = (healthy_bees^i)/(K(day)^i + healthy_bees^i);

% Function h (affects eclosion rate based on virus-carrying mites x(4))
h = exp(-x(4)*little_k(day));

% Function R
R = sigma_1(day) - sigma_2(day)*(x(2)/(healthy_bees + epsilon));

% Sum of all bees/mites
bees_sum = x(1) + x(2) + x(3);
mites_sum = x(4) + x(5);

% Proportions of bees
prop_hive = x(1)/(bees_sum + epsilon);
prop_forager = x(2)/(bees_sum + epsilon);
prop_healthy = healthy_bees/(bees_sum + epsilon);
prop_infected = x(3)/(bees_sum + epsilon);

% Mite logistic growth term
mite_logistic = 1 - (mites_sum/(alpha(day)*bees_sum + epsilon));

% Setting period of exposure 
temp_t = t;
while (temp_t > 365)
    temp_t = temp_t - 365;
end
if ((temp_t < 5) || (temp_t > 35)) 
    p(day) = 0;
end

% Setting varroacide treatment
if ~((30 <= t && t < 31) || (60 <= t && t < 61) || (90 <= t && t < 91))
    delta_1(day) = 0;
    delta_2(day) = 0;
    delta_3(day) = 0;
    delta_4(day) = 0;
    delta_5(day) = 0;
end
    
% 5-D complete system (eqs. 1-5)
dx = [
    mu(day)*g*h - beta_1(day)*x(4)*prop_hive - (d_1(day) + delta_1(day))*x(1) - gamma_1(day)*mites_sum*x(1) - x(1)*R;
    x(1)*R - beta_1(day)*x(4)*prop_forager - (p(day) + d_2(day) + delta_2(day))*x(2) - gamma_2(day)*mites_sum*x(2);
    beta_1(day)*x(4)*prop_healthy - (d_3(day) + delta_3(day))*x(3) - gamma_3(day)*mites_sum*x(3);
    r(day)*x(4)*mite_logistic + beta_2(day)*x(5)*prop_infected - beta_3(day)*x(4)*prop_healthy - delta_4(day)*x(4);
    r(day)*x(5)*mite_logistic - beta_2(day)*x(5)*prop_infected + beta_3(day)*x(4)*prop_healthy - delta_5(day)*x(5);
    ];
end