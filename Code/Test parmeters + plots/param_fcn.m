% Funtion to see seasonal paramerter values for a given time span
function parameter = param_fcn(t,parameter_vals)

parameter = zeros(numel(t),1); % Initalizing as zero vector

% for loop to calculate parameter(t)
    for i = 1:numel(t)

        % Initialzing curr_time to current time (day)
        if (size(t,1) == 1) % tspan has size (1xn)
            curr_time = t(1,i);
        elseif (size(t,2) == 1) % t in ode45 has size (nx1)
            curr_time = t(i,1);
        end
     
        % while loop to adjust temp_t if (t > 1 year)
        while (curr_time > 365) 
            curr_time = curr_time - 365;
        end

        % if-else statements to set parameter value
        if ((0 <= curr_time) && (curr_time < 91.25)) % Spring
            parameter(i,1) = parameter_vals(1,1);
        elseif ((91.25 <= curr_time) && (curr_time < 2*91.25)) % Summer
            parameter(i,1) = parameter_vals(2,1);
        elseif ((2*91.25 <= curr_time) && (curr_time < 3*91.25)) % Fall
            parameter(i,1) = parameter_vals(3,1);
        elseif ((3*91.25 <= curr_time) && (curr_time <= 4*91.25)) % Winter
            parameter(i,1) = parameter_vals(4,1);
        end
    end
