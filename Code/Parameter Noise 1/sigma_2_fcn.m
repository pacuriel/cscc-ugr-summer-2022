function sigma_2_par = sigma_2_fcn(t)

    sigma_2_vals = [0.75; 0.75; 0.75; 0.75;]; % Seasonal avgs.

    % Calling param_fcn to find periodic values
    sigma_2_par = param_fcn(t,sigma_2_vals);