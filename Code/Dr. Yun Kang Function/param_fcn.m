% Funtion to see seasonal paramerter values for a given time span
function parameter_vals = param_fcn(t,parameter_avgs)

% Initialzing parameter_vals as zero matrix; rows = days/parameter values,
% columns = parameters
parameter_vals = zeros(length(t),size(parameter_avgs,1));

% for loop to calculate parameter(t)
    for i = 1:numel(t)

        curr_time = t(i); % Setting current time
     
        % while loop to adjust temp_t if (t > 1 year)
        while (curr_time > 365) 
            curr_time = curr_time - 365;
        end

        % for loop to set parameter values
        for j = 1:size(parameter_avgs,1)

            % if-else statements to set parameter value
            if ((0 <= curr_time) && (curr_time < 91.25)) % Spring
                parameter_vals(i,j) = parameter_avgs(j,1);
            elseif ((91.25 <= curr_time) && (curr_time < 2*91.25)) % Summer
                parameter_vals(i,j) = parameter_avgs(j,2);
            elseif ((2*91.25 <= curr_time) && (curr_time < 3*91.25)) % Fall
                parameter_vals(i,j) = parameter_avgs(j,3);
            elseif ((3*91.25 <= curr_time) && (curr_time <= 4*91.25)) % Winter
                parameter_vals(i,j) = parameter_avgs(j,4);
            end
        end


           
    end
