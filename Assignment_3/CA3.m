% Part A) Newton 
syms x

F_sym = (x.*sin(pi.*x) + cos(pi.*x).*sin(pi.*x));
f1_sym = diff(F_sym, x);

F = matlabFunction(F_sym);
f1 = matlabFunction(f1_sym);
p0 = [0.77, 0.99];
tol = 1e-10; % Tolerance for accuracy 
max_iter = 1000; % Max iterations/break case

root = [0.78983, 1]; % Found in Part A

% Newton method
for j = 1:length(p0)
    iter_history = []; % Reset history for each root
    p_in = p0(j);
    error_prev = NaN; % Stores error from i-1 (previous iteration)
    
    for i = 1:max_iter
        % Newton's approx.
        p1 = p_in - (F(p_in) / f1(p_in));

        error_pin = abs(p_in - root(j));  % Current error
        error_p1 = abs(p1 - root(j)); % Next error

        % Convergence Analysis: Compute Alpha & Lambda when enough data exists
        if i > 2
            alpha = log(error_p1 / error_pin) / log(error_pin / error_prev); % Order of convergence
            lambda = error_p1 / (error_pin^alpha); %Asymptotic error constant
        else
            alpha = NaN;  % Not defined for first two iterations
            lambda = NaN;
        end

        iter_history = [iter_history; i, p1, error_pin, error_p1, alpha, lambda]; 

        if abs(p1 - p_in) < tol
            break; % Convergence met --> break loop
        end
        
        % Update error
        error_prev = error_pin;
        p_in = p1; % Overwrite for next iteration
    end 

    % Display results
    T = array2table(iter_history, 'VariableNames', {'Iteration', 'p_in', 'Error P_in', 'Error P1', 'Alpha', 'Lambda'});
    disp(T);
    disp(['For: ', num2str(p0(j)), ' Convergence reached at run: ', num2str(i)]);
    disp(['Convergence to: ', num2str(p1)])
    disp(' ')
end
