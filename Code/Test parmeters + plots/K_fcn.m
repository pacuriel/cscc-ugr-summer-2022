function K_par = K_fcn(t)

    K_vals = [8000; 12000; 8000; 6000;]; % Seasonal avgs.

        % Calling param_fcn to find periodic values
    K_par = param_fcn(t,K_vals);