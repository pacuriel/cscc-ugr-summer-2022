function mu_par = mu_fcn(t)

    mu_vals = [500; 1500; 500; 0;]; % Seasonal avgs.

    % Calling param_fcn to find periodic values
    mu_par = param_fcn(t,mu_vals);