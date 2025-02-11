clear all
close all 
clc 
syms h

%a)
%function
F = @(h) (exp(-h.^2) - 1)./(h.^2);

%symbolic func
F_sym = (exp(-h^2) - 1) / h^2;

%h values
h_values = linspace(0, 10, 1000);
results = F(h_values);

%limit as h -> 0
L = limit(F_sym, h, 0);
disp(['Limit as h -> 0: ', char(L)]);
%estimated limit from code = -1

%plotting
figure;
plot(h_values, results, 'LineStyle', '-');
xlabel('h');
ylabel('F(h)');
title('Plot of F(h) = (exp(-h^2) - 1) / h^2');
grid on;

%b)
T = taylor(F, h);
f_abs = abs(results - double(L));

%second plot
figure;
loglog(h_values, f_abs, 'LineStyle', '-', 'Color', 'b');
xlabel('h');
ylabel('|F(h) - L|');
title('Log-log plot of |F(h) - L| vs h');
grid on;

%linear region is where data behaviour replicates linear growth 
%i.e power-law relation
lin_region = (h_values > 1e-2 & h_values < 1e-1);

log_h = log(h_values(lin_region));
log_f_abs = log(f_abs(lin_region));

%using polyfit to esimate the line equation of the linear region;
%the slope of this line should be rougly = to p that im looking for
p = polyfit(log_h, log_f_abs, 1);
disp(['Estimated p from plot: ', num2str(p(1))]);