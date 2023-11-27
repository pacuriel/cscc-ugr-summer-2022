% Parameter values (Spring)
beta_1 = 0.1984; % Rate uninfected (hive and forager) bees become infected 
beta_2 = 0.1593; % Rate uninefected mites acquire virus
beta_3 = 0.04959; % Rate infected mites lose their virus to uninfected host 
d_1 = 0.02272; % Natural death rates for hive, forager, infected (resp.)
d_2 = 0.02272;
d_3 = 0.2;
mu = 500; % Max eclosion rate (worker bees born per day)
k = 0.000075; % Non-negative real number for function h
r = 0.0165; % Max mite birth rate
alpha = 0.4784; % Avg number of mites that can be sustained per bee
K = 8000; % Size of bee colony when eclosion rate is half of max
sigma_1 = 0.25; % Max rate of recruitment (no forgagers present)
sigma_2 = 0.75; % Rate of social inhibition
p = 0.13; % Forager mortality (homing failure)
gamma_1 = 10^(-7); % Rates at which mite kill bees (hive,forager,infected)
gamma_2 = gamma_1;
gamma_3 = 0.0000002; % gamma_3 > gamma_1,2

% Function parameters
% g: represents a large number of workers needed to care for brood
% h: eclosion rate is affected by presence of virus-carrying mites
% R: effect of social inhibition on the recruitment rate