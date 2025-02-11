h = logspace(-10, 0, 100); % h values from 10^-10 to 1
F_h = (exp(-h.^2) - 1) ./ h.^2;

% Plot
figure;
loglog(h, F_h, 'b-', 'LineWidth', 1.5);
xlabel('h');
ylabel('F(h)');
title('Plot of F(h) vs. h');
grid on;

L = -1; % Analytical limit
abs_diff = abs(F_h - L);

% Plot
figure;
loglog(h, abs_diff, 'r-', 'LineWidth', 1.5);
xlabel('h');
ylabel('|F(h) - L|');
title('Plot of |F(h) - L| vs. h');
grid on;