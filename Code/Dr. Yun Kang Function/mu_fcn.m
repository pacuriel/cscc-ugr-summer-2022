% Based on cosine r(t) function in Dr. Yun Knag's papers

function mu_val = mu_fcn(t)

gamma = 365; % Length of season
mu_0 = 625; % Base egg-laying rate (avg. of seasonal avgs.)
psi = 75; % Time/day of max egg-laying rate 
epsilon = 1; % Strength of season

mu_val = zeros(length(t),1); % Initialzing mu value

% for loop to update mu values
for i = 1:length(t)

    % Cosine function from Dr. Yun Kang's papers
    mu_val(i,1) = mu_0*(1 + epsilon*cos(2*pi*(t(i) - psi) / gamma));
end