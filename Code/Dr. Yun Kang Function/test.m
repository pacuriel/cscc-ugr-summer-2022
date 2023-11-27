tspan = 1:1:365;
%%
plot(tspan,weather_avg_temps)
xlim([0 365])
ylim([min(weather_avg_temps) max(weather_avg_temps)])
grid on
xlabel(['Time (days)'])
ylabel(['Avg. Temp. (C)'])
title('Malahat, Canada (2006)')