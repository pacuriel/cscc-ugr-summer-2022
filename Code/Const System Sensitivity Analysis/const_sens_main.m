clear all, close all, clc

% Initial Conditions
x0 = [12000; 5000; 0; 10; 50;];

% Time span
tspan = 0:1:600;

% Percent change in parameter value
percent_change = 0.1;

parameter_names = ["\beta_1"; "\beta_2"; "\beta_3"; "d_1"; "d_2"; "d_3"; "\mu"; 
    "little_k"; "r"; "\alpha"; "K"; "\sigma_1"; "\sigma_2"; "\gamma_1"; "\gamma_2";
    "\gamma_3"; "\delta_1"; "\delta_2"; "\delta_3"; "\delta_4"; "\delta_5"; "p";];

[avg_populations,percent_changes,days_survived] = sens_const_sys(x0,tspan,percent_change);