% Function to modify p-values based on weather data
function p_val = mod_p(t,p,weather_data)

% Year to use for daily weather data (1995-2021)
year = 2020; 

% year_ref to get desired year weather data
year_ref = year - 1995;

% Setting lower bound for foraging temperature
foraging_temp = 15;

% Getting desired year's weather data
row_min = 365*year_ref + 1; 
row_max = 365*(year_ref + 1);
year_data = weather_data(row_min:row_max,:);

% Initializing p_val to p
p_val = p;

% Winter = (0,46)U(365-45,365) in weather_data
% Starting weather_data during Spring based on above
data_iter = 46; 
n_year = 1; % Year counter

% while loop to update p-val
for i = 1:length(t)

    % If max. daily temp. < 15 degrees Celsius,
    % then bees do not forage (no homing failure). 
    if (year_data(data_iter,1) < foraging_temp)
        p_val(i) = 0;
    end
    
    % Increment data_iter
    data_iter = data_iter + 1; 
   
    % If reached end of weather data,
    % then go to start (2nd half of winter).
    if (data_iter == 366)
        data_iter = 1;
    end

    % Set data_iter back to start of Spring
    if (i == 365*n_year)
        data_iter = 46;
        n_year = n_year + 1;
    end 
end



