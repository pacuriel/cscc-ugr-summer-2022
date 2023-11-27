function parameter_vals = add_randomness(t,parameter_avgs)

% Initialzing parameter_vals as zero matrix; rows = days/parameter values,
% columns = parameters
parameter_vals = zeros(length(t),size(parameter_avgs,1));

% for loop to set values for each parameter on each day i
for i = 1:length(t)
    % Initializing year count to 0 (first year)
    year = 0;
    
    % Setting temp time as current time
    temp_t = t(i);
    
    % while loop to adjust temp_t if (t > 1 year) 
    while (temp_t > 365)
        year = year + 1; % Increment year count
        temp_t = temp_t - 365;
    end

    % Row value used to maintain seasonal averages
    row = 366*year + 1;

    % for loop to set values for each parameter j on day i
    for j=1:size(parameter_avgs,1)
        parameter_vals = get_rand_val(parameter_avgs,parameter_vals,i,j,temp_t,row);
   end
end

function param_vals = get_rand_val(parameter_avgs,parameter_vals,i,j,temp_t,row)

% Initialzing param_vals to parameter_vals
param_vals = parameter_vals;

% Stores seasonal avg.
avg = 0;

% Percent change for random parameter range
percent_change = 0.01;

% Percent of seasonal average for epsilon value
percent_epsilon = 0.005;

% Epsilon to keep avg parameter vals. within range of seasonal avg.
epsilon = 0;

% if-else statements to set parameter value
if ((0 <= temp_t) && (temp_t < 91.25)) % Spring
    avg = parameter_avgs(j,1);
    epsilon = avg*percent_epsilon;
    % Get a random number in range of +/-percent_change of seasonal avg.
    param_vals(i,j) = (avg - percent_change*avg) + ((avg...
        + percent_change*avg) - (avg - percent_change*avg))*rand;
    
    % While mean is not within +/-epsilon range of seasonal avg.,
    % obtain a new random value.
    while (mean(param_vals(row*1:i,j)) > avg + epsilon) || (mean(param_vals(row*1:i,j)) < avg - epsilon) 
        param_vals(i,j) = (avg - percent_change*avg) + ((avg...
            + percent_change*avg) - (avg - percent_change*avg))*rand;
    end
elseif ((91.25 <= temp_t) && (temp_t < 2*91.25)) % Summer
    avg = parameter_avgs(j,2);
    epsilon = avg*percent_epsilon;
    % Get a random number in range of +/-percent_change of seasonal avg.
    param_vals(i,j) = (avg - percent_change*avg) + ((avg...
        + percent_change*avg) - (avg - percent_change*avg))*rand;
    
    % While mean is not within +/-epsilon range of seasonal avg.,
    % obtain a new random value.
    while (mean(param_vals(row*93:i,j)) > avg + epsilon) || (mean(param_vals(row*93:i,j)) < avg - epsilon) 
        param_vals(i,j) = (avg - percent_change*avg) + ((avg...
            + percent_change*avg) - (avg - percent_change*avg))*rand;
    end
elseif ((2*91.25 <= temp_t) && (temp_t < 3*91.25)) % Fall
    avg = parameter_avgs(j,3);
    epsilon = avg*percent_epsilon;
    % Get a random number in range of +/-percent_change of seasonal avg.
    param_vals(i,j) = (avg - percent_change*avg) + ((avg...
        + percent_change*avg) - (avg - percent_change*avg))*rand;
    
    % While mean is not within +/-epsilon range of seasonal avg.,
    % obtain a new random value.
    while (mean(param_vals(row*184:i,j)) > avg + epsilon) || (mean(param_vals(row*184:i,j)) < avg - epsilon) 
        param_vals(i,j) = (avg - percent_change*avg) + ((avg...
            + percent_change*avg) - (avg - percent_change*avg))*rand;
    end
elseif ((3*91.25 <= temp_t) && (temp_t <= 4*91.25)) % Winter
    avg = parameter_avgs(j,4);
    epsilon = avg*percent_epsilon;
    % Get a random number in range of +/-percent_change of seasonal avg.
    param_vals(i,j) = (avg - percent_change*avg) + ((avg...
        + percent_change*avg) - (avg - percent_change*avg))*rand;
    
    % While mean is not within +/-epsilon range of seasonal avg.,
    % obtain a new random value.
    while (mean(param_vals(row*275:i,j)) > avg + epsilon) || (mean(param_vals(row*275:i,j)) < avg - epsilon) 
        param_vals(i,j) = (avg - percent_change*avg) + ((avg...
            + percent_change*avg) - (avg - percent_change*avg))*rand;
    end
end  