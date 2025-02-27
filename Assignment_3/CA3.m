% Newton 
syms x

F_sym = (x.*sin(pi.*x) + cos(pi.*x).*sin(pi.*x));
f1_sym = diff(F_sym, x);

iter_history = [];

F = matlabFunction(F_sym);
f1 = matlabFunction(f1_sym);
p0 = [0.77, 0.99];
tol = 1e-7; %tolerance for accuracy 
max_iter = 1000; %max iterations/break case

%newton method
for j = 1:length(p0)
    p_in = p0(j);
    for i = 1:max_iter
        p1 = p_in - (F(p_in)./f1(p_in));
        error = abs(p1 - p_in); % Compute error
        iter_history = [iter_history; i, p1, F(p1), error]; 
        if abs(p1 - p_in) < tol
            break; %convergence met --> break loop
        end
        p_in = p1; %overwriting to continue next iteration
    end 
    disp(['For: ', num2str(p0(j)), ' Convergence reached at run: ', num2str(i)]);
    disp(['Convergence to: ', num2str(p1)])
    disp(' ')
    T = array2table(iter_history, 'VariableNames', {'Iteration', 'p_in', 'F(p_in)', 'Error'});
    disp(T);

    %figure;
    %semilogy(iter_history(:,1), iter_history(:,4), '-o', 'LineWidth', 1.5);
    %xlabel('Iteration');
    %ylabel('Error (log scale)');
    %title(['Convergence of Newton%s Method (p0 = ', num2str(p0(j)), ')']);
    %grid on;
end