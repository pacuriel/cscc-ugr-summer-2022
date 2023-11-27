function sigma_1_par = sigma_1_fcn(t)

    sigma_1_vals = [0.25; 0.25; 0.25; 0;]; % Seasonal avgs.

        % Calling param_fcn to find periodic values
    sigma_1_par = param_fcn(t,sigma_1_vals);