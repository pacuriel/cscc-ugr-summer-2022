function delta_3 = delta_3_fcn(t)

    temp_t = t; % Temporary t-value 
    while (temp_t > 365)
        temp_t = temp_t - 365;
    end
    % Updating delta values during specifc days of Spring
    if ((30 <= temp_t) && (temp_t < 31)) ||  ((60 <= temp_t) && (temp_t < 61)) || ((90 <= temp_t) && (temp_t < 91))
        delta_3 = 0.005;
    else 
        delta_3 = 0;
    end


 