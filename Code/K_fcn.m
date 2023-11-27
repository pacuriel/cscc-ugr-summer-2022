function K_par = K_fcn(t)

    K_vals = [8000; 12000; 8000; 6000;]; % Seasonal avgs.

    % Increased brood maintenance coefficient
    K_fig_3e = [14000; 16000; 14000; 9000;]; 

        % Calling param_fcn to find periodic values
    K_par = param_fcn(t,K_vals);