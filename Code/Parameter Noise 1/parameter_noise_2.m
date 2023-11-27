function parameter_value = parameter_noise_2(table_1_val)

% Percent change for interval of random number
percent_change = 0.005;

% Standard deviation for normal distribution (5% of default value)
sd = percent_change*table_1_val;

% Setting parameter value as rand. number using normal distribution (mean = table_1_val)
parameter_value = normrnd(table_1_val,sd);