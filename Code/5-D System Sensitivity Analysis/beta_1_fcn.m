function beta_1_par = beta_1_fcn(t)

    beta_1_vals = [0.1984; 0.1460; 0.1900; 0.03384;]; % Seasonal avgs.

    % Calling param_fcn to find periodic values
    beta_1_par = param_fcn(t,beta_1_vals);