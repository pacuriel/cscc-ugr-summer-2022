function alpha_par = alpha_fcn(t)

    alpha_vals = [0.4784; 0.5; 0.5; 0.4784;]; % Seasonal avgs.

        % Calling param_fcn to find periodic values
    alpha_par = param_fcn(t,alpha_vals);