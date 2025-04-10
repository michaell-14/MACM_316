close all;
clear all;
n_val = [2, 6, 10, 20, 50, 100, 200, 500]; % Various n values
m_val = [3, 5]; % Different m values

%Piecewise function definitions using cell arrays
f = {@(x) x.^2, @(x, m) x.^2 + ((x - 0.5).^m)};
a = [0, 0.5];
b = [0.5, pi];

syms x m real
f1_sym = x^2;
f2_sym = x^2 + (x - 0.5)^m;

m3_error = zeros(size(n_val));
m5_error = zeros(size(n_val));
h_3 = zeros(size(n_val));
h_5 = zeros(size(n_val));

% For each m value, calculate exact integral and error for each n
for m = 1:length(m_val)
    disp('-----------------------------------------')
    disp(['m = ', num2str(m_val(m))])
    I_exact = integral(f{1}, a(1), b(1)) + integral(@(x) f{2}(x, m_val(m)), a(2), b(2));
    Df1 = diff(f1_sym, 4)
    Df2 = diff(f2_sym, 4)

    %Convert symbolic to function
    Df1_func = matlabFunction(Df1);
    Df2_func = matlabFunction(Df2);

    for n = 1:length(n_val)
        I_total = 0;

        for t = 1:2
            h = (b(t)-a(t))/n_val(n);
            if m_val(m) == 3
                h_3(n) = h;
            else
                h_5(n) = h;
            end

            x = linspace(a(t), b(t), n_val(n) + 1);

            if t == 1
                f_x = f{1}(x);
                f4_max = 0; % Fourth derivative of x^2 is zero
            else
                f_x = f{2}(x, m_val(m));
                %Find maximum of fourth derivative on interval
                if m_val(m) == 3
                    f4_max = 0; % Fourth derivative is zero for m = 3
                elseif m_val(m) > 3
                    %Sample points to find approximate max of fourth derivative
                    test_points = linspace(a(t), b(t), 100);
                    f4_values = abs(Df2_func(test_points, m_val(m)));
                    f4_max = max(f4_values);
                else
                    %For m < 3, evaluate at the endpoint
                    f4_max = abs(Df2_func(b(t), m_val(m)));
                end
            end

            % Composite Simpson Rule
            sum_odd = sum(f_x(2:2:end-1));
            sum_even = sum(f_x(3:2:end-2));
            I_t = (h/3) * (f_x(1) + 4 * sum_odd + 2 * sum_even + f_x(end));
            I_total = I_total + I_t;

            % Error for this subinterval
            abs_simpson_error = (b(t)-a(t))^5/(180*n_val(n)^4) * f4_max; %for m = 3 should be 0
        end

        %error for plotting
        actual_error = abs(I_exact - I_total);
        if m_val(m) == 3
            m3_error(n) = actual_error;
        else
            m5_error(n) = actual_error;
        end

        disp(['For m = ', num2str(m_val(m)), '; n = ', num2str(n_val(n)), '; I = ', num2str(I_total)])
        disp(['I exact = ', num2str(I_exact)])
        disp(['Abs Simpsons Error: ', num2str(abs_simpson_error)])
        disp(['Abs Error: ', num2str(actual_error)])
        disp(' ')
    end
end

%-----------DISPLAYS AND ANALYSIS-----------------

disp(['N Values: ',num2str(n_val)])
disp(' ')
disp(['H for m = 3: ',num2str(h_3)])
disp(['Abs error for m = 3',num2str(m3_error)])
disp('------------------------------------------------------')
disp(['H for m = 5: ',num2str(h_5)])
disp(['Abs error for m = 5',num2str(m5_error)])

disp(' ')
disp('POLYFIT RESULTS:')

%plot for m=3
fig1 = figure;
loglog(h_3, m3_error, '-o', 'color', 'b')
xlabel('Subintervals Size (h)');
ylabel('Error');
title('Error vs Subinterval sizes for m = 3');
grid on;
print(fig1, "Error_vs_h3.png", '-dpng');

log_h3 = log(h_3);
log_error_m3 = log(m3_error);
p = polyfit(log_h3, log_error_m3, 1);
%disp(['Order of Convergence m = 3: h^',num2str(p(1))]); 
%For m = 3 the error is 0 as its a flat linear line, equating to O(1) convergence
disp('Order of Convergence m = 3: O(1)')

%plot for m=5
fig2 = figure;
loglog(h_5, m5_error, '-s','color', 'r')
xlabel('Subintervals Size (h)');
ylabel('Error');
title('Error vs Subinterval sizes for m = 5');
grid on;
print(fig2, "Error_vs_h5.png", '-dpng');

log_h5 = log(h_5);
log_error_m5 = log(m5_error);
q = polyfit(log_h5, log_error_m5, 1);
disp(['Order of Convergence m = 5: O(h^',num2str(q(1)), ')']);