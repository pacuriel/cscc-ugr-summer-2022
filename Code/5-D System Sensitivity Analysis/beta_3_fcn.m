function beta_3_par = beta_3_fcn(t)

    beta_3_vals = [0.04959; 0.03721; 0.04750; 0.008460;]; % Seasonal avgs.

    % Calling param_fcn to find periodic values
    beta_3_par = param_fcn(t,beta_3_vals);