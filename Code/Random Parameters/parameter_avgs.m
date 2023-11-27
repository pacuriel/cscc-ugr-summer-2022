function parameter_vals = parameter_avgs(tspan)

% Seasonal parameter averages
alpha_vals = [0.4784 0.5 0.5 0.4784];
beta_1_vals = [0.1984 0.1460 0.1900 0.03384];
beta_2_vals = [0.1593 0.1460 0.1489 0.04226];
beta_3_vals = [0.04959 0.03721 0.04750 0.008460];
d_1_vals = [0.02272 0.04 0.02272 0.005263];
d_2_vals = [0.02272 0.04 0.02272 0.005263];
d_3_vals = [0.2 0.2 0.2 0.005300];
delta_1 = [0.005 0 0 0];
delta_2 = [0.005 0 0 0];
delta_3 = [0.005 0 0 0];
delta_4 = [0.5 0 0 0];
delta_5 = [1 0 0 0];
K_vals = [8000 12000 8000 6000];
little_k_vals = [0.000075 0.00003125 0.000075 0];
mu_vals = [500 1500 500 0];
r_vals = [0.0165 0.0165 0.0045 0.0045];
sigma_1_vals = [0.25 0.25 0.25 0];
sigma_2_vals = [0.75 0.75 0.75 0.75];
p_vals = [0.1 0 0 0];
gamma_1_vals = [10^(-7) 10^(-7) 10^(-7) 10^(-7)];
gamma_2_vals = [10^(-7) 10^(-7) 10^(-7) 10^(-7)];
gamma_3_vals = [0.0000002 0.0000002 0.0000002 0.0000002];

% Matrix where row = new parameter, column = seasonal avgs.
parameter_vec = [beta_1_vals; beta_2_vals; beta_3_vals; d_1_vals; d_2_vals;
    d_3_vals; mu_vals; little_k_vals; r_vals; alpha_vals; K_vals; 
    sigma_1_vals; sigma_2_vals; gamma_1_vals; gamma_2_vals; gamma_3_vals;
    delta_1; delta_2; delta_3; delta_4; delta_5; p_vals;];

% Calling add_randomness to obtain a matrix with random values while
% maintaining seasonal avgs. 
parameter_vals = add_randomness(tspan,parameter_vec);