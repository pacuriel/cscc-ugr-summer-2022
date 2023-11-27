function beta_2_par = beta_2_fcn(t)

    beta_2_vals = [0.1593; 0.1460; 0.1489; 0.04226;]; % Seasonal avgs.

    % Calling param_fcn to find periodic values
    beta_2_par = param_fcn(t,beta_2_vals);