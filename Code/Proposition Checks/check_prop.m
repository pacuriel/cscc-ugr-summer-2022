% Returns true if stable and false if unstable.
function test_stability = check_prop(t,x,beta_3,d_1,d_2,mu,r,K,sigma_1,sigma_2,i,delta_i,p)

% Checking conditions (i) and (ii) in Prop. 3.6
if (r > beta_3 + delta_i(4)) || (r > delta_i(5))
    time = t(1);
    stable = false;
    output(time,stable);
    return
end



% If above did not return, then system inherits stability of 2-D system.
% (Prop. 3.1-3.3)

% Equation 10 in 2017 paper
F_term = (sigma_1 - sigma_2 - p - d_2)/(p + d_2);
F = (F_term + sqrt(F_term^2 + 4*sigma_1/(p + d_2))) / 2;
a = -mu / ((sigma_2*F)/(1+F) - d_1 - sigma_1);

% Inequality 11
if ((sigma_2*F)/(1+F) - d_1 - sigma_1) > 0
    disp("Inequality 11 holds. (0,0) is the only equilibrium.")
else
    disp("Inequality 11 does not hold.")
    if a > (K*i/(1+F))*(i-1)^(-1 + 1/i)
        disp("There exists two positive bee-only equilibria.")
    end
end

for idx = 1:length(t)
    % Inequality 21
    left_ineq_21 = (i*mu*K^i*x(idx,1)^(i-1)*(1 + F)^(i-1)) / (K^i + x(idx,1)^i*(1 + F)^i)^2;
    right_ineq_21 = (d_1 + (p + d_2)*F) / (1 + F)^2; 

    % Checking if inequality does not hold (unstable)
    if (left_ineq_21 > right_ineq_21)
        time = t(idx);
        stable = false;
        output(time,stable)
        return
    end
end

% Conditions were met (stable)
time = t(length(t));
stable = true;
output(time,stable)

function disp_output = output(time,stable)
    if (stable)
        A = ['System is stable.'];
        disp(A);
    else 
        A = ['System is unstable at t = ',num2str(time),'.'];
        disp(A);
end