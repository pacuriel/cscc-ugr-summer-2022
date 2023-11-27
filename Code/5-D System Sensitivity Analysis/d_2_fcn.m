function d_2_par = d_2_fcn(t)

    d_2_vals = [0.02272; 0.04; 0.02272; 0.005263;]; % Seasonal avgs.

    % Calling param_fcn to find periodic values
    d_2_par = param_fcn(t,d_2_vals);