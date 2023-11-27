function parameter_vals = set_parameters(t,parameter_vec)

% Initialzing parameter_vals as zero matrix; rows = days/parameter values,
% columns = parameters
parameter_vals = zeros(length(t),size(parameter_vec,1));

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

    % for loop to set values for each parameter j on day i
    for j=1:size(parameter_vec,1)
        parameter_vals = get_rand_val(parameter_vec,parameter_vals,i,j,temp_t);
   end
end

function param_vals = get_rand_val(parameter_avgs,parameter_vals,i,j,temp_t)

% Initialzing param_vals to parameter_vals
param_vals = parameter_vals;

% if-else statements to set parameter value
if ((0 <= temp_t) && (temp_t < 91.25)) % Spring
    param_vals(i,j) = parameter_avgs(j,1);
elseif ((91.25 <= temp_t) && (temp_t < 2*91.25)) % Summer
    param_vals(i,j) = parameter_avgs(j,2);
elseif ((2*91.25 <= temp_t) && (temp_t < 3*91.25)) % Fall
    param_vals(i,j) = parameter_avgs(j,3);
elseif ((3*91.25 <= temp_t) && (temp_t <= 4*91.25)) % Winter
    param_vals(i,j) = parameter_avgs(j,4);
end  