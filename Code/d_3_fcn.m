function d_3_par = d_3_fcn(t)

    d_3_vals = [0.2; 0.2; 0.2; 0.005300;]; % Seasonal avgs.

    % Calling param_fcn to find periodic values
    d_3_par = param_fcn(t,d_3_vals);