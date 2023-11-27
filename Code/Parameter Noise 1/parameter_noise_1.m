function parameter_value = parameter_noise_1(table_1_val)

% Percent change for interval of random number
percent_change = 0.1;

% Setting parameter value as rand. number within interval of +/-10% Table 1 value 
parameter_value = (table_1_val - percent_change*table_1_val) + ((table_1_val + percent_change*table_1_val) - (table_1_val - percent_change*table_1_val))*rand;