A = zeros(500,1);

avg = 500;

percent_change = 0.01;

percent = 0.005;

epsilon = avg*percent;

for i=1:size(A)
    A(i) = (avg - percent_change*avg) + ((avg + percent_change*avg) - (avg - percent_change*avg))*rand;
    while (mean(A(1:i)) > avg + epsilon) && (mean(A(1:i)) < avg - epsilon)
        A(i) = (avg - percent_change*avg) + ((avg + percent_change*avg) - (avg - percent_change*avg))*rand;
    end
    avg_check = mean(A(1:i));
end

mean(A)