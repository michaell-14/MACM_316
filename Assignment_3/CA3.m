% Newton 
syms x

F_sym = (x.*sin(pi.*x) + cos(pi.*x).*sin(pi.*x));
f1_sym = diff(F_sym, x);

F = matlabFunction(F_sym);
f1 = matlabFunction(f1_sym);

p0 = [0.77, 0.99];
tol = 1e-10; %tolerance for accuracy 
max_iter = 1000; %max iterations/break case

%newton method
for j = 1:length(p0)
    p_in = p0(j);
    for i = 1:max_iter
        p1 = p_in - (F(p_in)./f1(p_in));
        if abs(p1 - p_in) < tol
            break; %convergence met --> break loop
        end
        p_in = p1; %overwriting to continue next iteration
    end 
    disp(['For: ', num2str(p0(j)), ' Convergence reached at: ', num2str(i)]);
end
