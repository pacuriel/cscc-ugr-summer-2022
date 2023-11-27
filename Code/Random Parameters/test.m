x = -3:3; 
y = [-1 -1 -1 0 1 1 1]; 
xq1 = -3:10^(-7):3;
p = pchip(x,y,xq1);

%%
A = zeros(500,1);

% Seasonal average of parameter
avg = 500;

% Percent change for random parameter range
percent_change = 0.01;

% Percent of seasonal average for epsilon value
percent = 0.005;

% Epsilon to keep avg parameter vals. within range of seasonal avg.
epsilon = avg*percent;

% for loop to obtain random parameter values
for i=1:size(A)
    % Get a random number in range of +/-percent_change of seasonal avg.
    A(i) = (avg - percent_change*avg) + ((avg + percent_change*avg) - (avg - percent_change*avg))*rand;
    
    mean(A(1:i))

    % While mean(A) is not within +/-epsilon range of seasonal avg.,
    % obtain a new random value.
    while (mean(A(1:i)) > avg + epsilon) || (mean(A(1:i)) < avg - epsilon)
        A(i) = (avg - percent_change*avg) + ((avg + percent_change*avg) - (avg - percent_change*avg))*rand;
    end
    
    

    % Below  chunk leads to infinite while loop (reason for epsilon range)
%     while (mean(A(1:i)) ~= avg)
%         A(i) = (avg - percent_change*avg) + ((avg + percent_change*avg) - (avg - percent_change*avg))*rand;
%     end

end

mean(A)