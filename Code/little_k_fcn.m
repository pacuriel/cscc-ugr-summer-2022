function k_par = little_k_fcn(t)

    k_vals = [0.000075; 0.00003125; 0.000075;]; % Seasonal avgs.

    k_par = zeros(numel(t),1); % Initalizing as zero vector

    % for loop to calculate parameter(t)
    for i = 1:numel(t)

        % Initialzing curr_time to current time (day)
        if (size(t,1) == 1)
            curr_time = t(1,i);
        elseif (size(t,2) == 1)
            curr_time = t(i,1);
        end

        % while loop to adjust temp_t if (t > 1 year)
        while (curr_time > 365) 
            curr_time = curr_time - 365;
        end

        % if-else statements to set parameter value
        if ((0 <= curr_time) & (curr_time < 91.25)) % Spring
            k_par(i,1) = k_vals(1,1);
        elseif ((91.25 <= curr_time) & (curr_time < 2*91.25)) % Summer
            k_par(i,1) = k_vals(2,1);
        elseif ((2*91.25 <= curr_time) & (curr_time < 3*91.25)) % Fall
            k_par(i,1) = k_vals(3,1);
        end
    end
