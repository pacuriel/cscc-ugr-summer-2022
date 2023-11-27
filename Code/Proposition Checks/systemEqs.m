% Function to input model equations
function dx = systemEqs(t,x,beta_1,beta_2,beta_3,d_1,d_2,d_3,mu,k,r,alpha, ...
    K,sigma_1,sigma_2,gamma_1,gamma_2,gamma_3,delta_i,p,i)
% Function g (w/ above i-value)
g = ((x(1) + x(2))^i)/(K^i + (x(1) + x(2))^i);

% Function h (affects eclosion rate based on virus-carrying mites x(4))
h = exp(-x(4)*k);

% Function R
R = sigma_1 - sigma_2*(x(2)/(x(1) + x(2)));

% Applying varroacide treatment
temp_t = t; % Temporary t-value 
while (temp_t > 365)
    temp_t = temp_t - 365;    
end
% Updating delta values to zero if not in specific 3 days
if ~(((30 <= temp_t) && (temp_t < 31)) ||  ((60 <= temp_t) && (temp_t < 61)) || ((90 <= temp_t) && (temp_t < 91)))
    delta_i = zeros(5,1);
end 

% Setting period of exposure 
if ((t < 5) || (t > 35)) 
    p = 0;
end

% 5-D complete system (eqs. 1-5)
dx = [
    mu*g*h - beta_1*x(4)*(x(1)/(x(1) + x(2) + x(3))) - (d_1 + delta_i(1))*x(1) - gamma_1*(x(4) + x(5))*x(1) - x(1)*R;
    x(1)*R - beta_1*x(4)*(x(2)/(x(1) + x(2) + x(3))) - (p + d_2 + delta_i(2))*x(2) - gamma_2*(x(4) + x(5))*x(2);
    beta_1*x(4)*(x(1) + x(2))/(x(1) + x(2) + x(3)) - (d_3 + delta_i(3))*x(3) - gamma_3*(x(4) + x(5))*x(3);
    r*x(4)*(1 - ((x(4) + x(5))/(alpha*(x(1) + x(2) + x(3))))) + beta_2*x(5)*(x(3)/(x(1) + x(2) + x(3))) - beta_3*x(4)*(x(1) + x(2))/(x(1) + x(2) + x(3)) - delta_i(4)*x(4);
    r*x(5)*(1 - ((x(4) + x(5))/(alpha*(x(1) + x(2) + x(3))))) - beta_2*x(5)*(x(3)/(x(1) + x(2) + x(3))) + beta_3*x(4)*(x(1) + x(2))/(x(1) + x(2) + x(3)) - delta_i(5)*x(5);
    ];
end