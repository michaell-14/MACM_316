close all;
clear all;
n_val = [10, 20, 50, 100, 200, 500]; % Various n values
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

for m = 1:length(m_val)
    disp('-----------------------------------------')
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
            x = linspace(a(t), b(t), n_val(n) + 1); %Generating points 
            
            if t == 1
                f_x = f{1}(x); %f1(x) for first interval
            else
                f_x = f{2}(x, m_val(m)); %f2(x) for second interval
            end
            
            %Composite Simpson Rule formula
            sum_odd = sum(f_x(2:2:end-1));  %sum odd indices
            sum_even = sum(f_x(3:2:end-2)); %sum even indices
            
            I_t = (h/3) * (f_x(1) + 4 * sum_odd + 2 * sum_even + f_x(end));
            I_total = I_total + I_t;
            error = abs((I_exact-I_total)/I_exact);
            %error = max(abs(f4x))*((b(t)-a(t))/180)*h^4;
            if m_val(m) == 3
                m3_error(n) = error;
            elseif m_val(m) == 5
                m5_error(n) = error;
            end
        end
        disp(['For m = ', num2str(m_val(m)), '; n = ', num2str(n_val(n)), '; I = ', num2str(I_total)])
        disp(['I exact = ', num2str(I_exact)])
        disp(['Abs Error: ', num2str(error)])
        disp(' ') %just for new line clarity
    end
end
m3_error
m5_error