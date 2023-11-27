function delta_i_par = delta_i_fcn(t)

    temp_t = t; % Temporary t-value 
    while (temp_t > 365)
        temp_t = temp_t - 365;
    end
    % Updating delta values during specifc days of Spring
    if ((30 <= temp_t) && (temp_t < 31)) ||  ((60 <= temp_t) && (temp_t < 61)) || ((90 <= temp_t) && (temp_t < 91))
        delta_i_par(1) = 0.005;
        delta_i_par(2) = delta_i_par(1);
        delta_i_par(3) = delta_i_par(1);
        delta_i_par(4) = 0.5;
        delta_i_par(5) = 1.5;
    else 
        delta_i_par(1) = 0;
        delta_i_par(2) = delta_i_par(1);
        delta_i_par(3) = delta_i_par(1);
        delta_i_par(4) = delta_i_par(1);
        delta_i_par(5) = delta_i_par(1);
    end


 