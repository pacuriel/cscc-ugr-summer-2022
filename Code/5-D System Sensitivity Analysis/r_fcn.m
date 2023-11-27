function r_par = r_fcn(t)

    r_vals = [0.0165; 0.0165; 0.0045; 0.0045;]; % Seasonal avgs.

    % Calling param_fcn to find periodic values
    r_par = param_fcn(t,r_vals);