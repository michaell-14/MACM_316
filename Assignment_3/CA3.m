% Part A) Newton 
syms x

F_sym = (x.*sin(pi.*x) + cos(pi.*x).*sin(pi.*x));
f1_sym = diff(F_sym, x); %first derivativve
f2_sym = diff(f1_sym, x);%second derivative

F = matlabFunction(F_sym);
f1 = matlabFunction(f1_sym);
f2 = matlabFunction(f2_sym);

p0 = [0.77, 0.99]; %Intitial Points

%Stuff for double root%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%x_test = 0.7893;
x_test = 1; 
F_value = F(x_test);
f1_value = f1(x_test);
f2_value = f2(x_test);

if abs(F_value) < 1e-10 && abs(f1_value) < 1e-10 && abs(f2_value) > 1e-10
    fprintf('x = 1 is a double root\n\n'); %mulipicity == 2
elseif abs(F_value) < 1e-10 && abs(f1_value) > 1e-10
    fprintf('x = 1 is a simple root\n\n'); %mulipicity == 1
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tol = 1e-6; % Tolerance for accuracy 
max_iter = 100; % Max iterations/break case

root = [0.78983, 1]; %Found in Part A

% Newton method
for j = 1:length(p0)
    iter_history = []; %Reset for each root
    p_in = p0(j);
    error_ls = [];
    
    for i = 1:max_iter
        % Newton's approx.
        p1 = p_in - (F(p_in) / f1(p_in));

        error_pin = abs(p_in - root(j));  % Current error
        error_p1 = abs(p1 - root(j)); % Next error

        error_ls = [error_ls; error_pin]; %appends errors

        % Convergence Analysis: Log-Log regression
        if length(error_ls) > 2
               log_error = log(error_ls(1:end-1));
               log_error_next = log(error_ls(2:end));

               %line --> log(error_n+1) = alpha*(error_n) + C
               p = polyfit(log_error, log_error_next, 1);
               alpha = p(1); % Slope is the estimated order of convergence
               lambda = exp(p(2)); % Asymptotic error constant
        else
            alpha = NaN;
            lambda = NaN;
        end

        iter_history = [iter_history; i, p1, error_pin, error_p1, alpha, lambda]; 
        %iter_history = [iter_history; i, p1, error_pin, error_p1];

        if abs(p1- p_in) < tol 
            break; % Convergence met --> break loop
        end

        p_in = p1; % Overwrite for next iteration
    end 

    % Display results
    T = array2table(iter_history, 'VariableNames', {'Iteration', 'p_in', 'Error P_in', 'Error P1', 'Alpha', 'Lambda'});
    disp(T);
    disp(['For: ', num2str(p0(j)), ' Convergence reached at run: ', num2str(i)]);
    disp(['Convergence to: ', num2str(p1)])
    disp(' ')
end
