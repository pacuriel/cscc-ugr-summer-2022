clear all, close all, clc

% Sensitivity analysis on full model

% Initial Conditions
x0 = [12500; 5000; 0; 10; 50;];

% Time span (days)
tspan = 0:1:3650; % 10 years

% Percent change in parameter value
percent_change = 0.1;

% Calling function to obtain sensitivity analysis data
[avg_populations,percent_changes,days_survived] = sens_periodic_5D(x0,tspan,percent_change);

%% Attempting to create tornado diagram

% Names of parameters in order analyzed
parameter_names = ["\beta_1"; "\beta_2"; "\beta_3"; "d_1"; "d_2"; "d_3"; "\mu"; 
    "k"; "r"; "\alpha"; "K"; "\sigma_1"; "\sigma_2"; "\gamma_1"; "\gamma_2";
    "\gamma_3"; "\delta_1"; "\delta_2"; "\delta_3"; "\delta_4"; "\delta_5"; "p";];
param_names = categorical(parameter_names);

%%

%%% Take data to excel and sort for tornado diagram
%%% Manually change data in percent_changes,param_names

% Plotting tornado diagram
% barh(percent_changes(1:22,1:2),'stacked')
% yticks(1:1:22)
% yticklabels(param_names)
% xlim([-85 240])
% grid on
% legend('-10%','+10%')
% title('%-Change in Avg. Total Population (p=0.166)')





% % Sorting days_survived by most days_survived
% days_survived_sorted = sortrows(days_survived,3,"descend");
% 
% Creating bar plot of days_survived
bar(param_names,days_survived(:,1:2))
title('Colony Lifespan (p=0.2)')
legend('-10%','+10%')
ylim([0 3651])