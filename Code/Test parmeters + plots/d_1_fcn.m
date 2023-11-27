function d_1_par = d_1_fcn(t)

    d_1_vals = [0.02272; 0.04; 0.2; 0.005263;]; % Seasonal avgs.

        % Calling param_fcn to find periodic values
    d_1_par = param_fcn(t,d_1_vals);