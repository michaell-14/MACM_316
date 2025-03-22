close all;
clear all;

% Processing:
% x_n+1 = D_inv.*(L+U).*x + D'.*b_n
% Refactoring above: ^^^
% x_new = D_inv * (b_n - (L + U) * x);

tol = 1e-4;
c = [2, 4];  %[Part A, Part B]
%n = [5, 10, 25, 50, 100]; %low test
n = [10, 25, 50, 100, 500]; %base test
%n = [50, 100, 500, 1000, 5000]; %stress testing
max_iter = 10000;  % Safety limit for iterations

%empty ls for iterations, for plotting
iter_c2 = zeros(size(n));
iter_c4 = zeros(size(n));

for i = 1:length(c)
    disp('-----------------------------------')
    disp(['Jacobi for c=', num2str(c(i))])
    disp('-----------------------------------')
    for j = 1:length(n)
        % Construct the tridiagonal matrix and right-hand side
        h = 1 / (n(j) + 1);
        main_diag = (c(i)/h^2) * ones(n(j),1);
        off_diag = (-1/h^2) * ones(n(j),1);
        A_n = spdiags([off_diag main_diag off_diag], [-1 0 1], n(j), n(j));
        b_n = ones(n(j), 1);

        %Declaratives:
        D = diag(A_n);
        D_inv = diag(1./D);
        L = tril(A_n, -1); %lower tri 
        U = triu(A_n, 1); %upper tri
        x = zeros(n(j), 1);

        x_exact = A_n \ b_n;
        iter = 0;  %reset iteration

        %diag dominant check
        D_dd = abs(D);
        non_diag_sum = sum(abs(A_n), 2) - D_dd;
        while true
            if all(D_dd >= non_diag_sum)
                disp('diagonally dominant');
                break;
            else
                disp('NOT diagonally dominant');
                break;
            end
        end

        %Spectral Radius check
        spectral_radius = max(abs(eig(D_inv * (L + U))));
        disp(['Spectral Radius: ', num2str(spectral_radius)]);

        %Jacobi Iteration
        while iter < max_iter
            x_new = D_inv * (b_n - (L + U) * x); % Jacobi update
            rel_error = norm(x_new - x_exact, 2) / norm(x_exact, 2);  % Relative error
            iter = iter + 1;

            %iteration appending
            if c(i) == 2
                iter_c2(j) = iter;
            elseif c(i) == 4
                iter_c4(j) = iter;
            end

            if rel_error < tol
                disp(['c = ', num2str(c(i)), ', n = ', num2str(n(j)), ', Iterations: ', num2str(iter)]);
                disp(['Error: ', num2str(rel_error)])
                disp(" ")
                break
            end
            x = x_new; %Update x
        end
        if iter == max_iter
            warning(['Jacobi method did not converge for c = ', num2str(c(i)), ', n = ', num2str(n(j))]);
            %Safety, meaning max iterations has been hit before convergence
        end
    end
end
disp(iter_c2)
disp(iter_c4)

%Calculating growth rate for C = 2, EMPERICALLY
valid_indices = 1:3; % Keep only the first three elements
ic2 = iter_c2(valid_indices);
n_valid = n(valid_indices);%plotting n vs iter for c=2
r = 0;
for s = 1:length(ic2)
    r = r + ic2(s)/(n_valid(s)^2);
end
c2_r = r/length(ic2);
disp(['Growth Rate for Jacobi c = 2 --> ', num2str(c2_r)])

fig1 = figure;
% no log-log regression as invalid convergence cases
% the failure to converge leads invalid regression line slope
% slope for c = 2 is obtained through imperical methods (ABOVE)
loglog(n, iter_c2,'-s','color', 'r');
xlabel('Matrix size (n)');
ylabel('# of iterations');
title('Iterations require to solve nxn matricies for C = 2');
grid on;
print(fig1, "Jacobi_C2.png", '-dpng'); 


%plotting n vs iter for c=4
fig2 = figure;
loglog(n, iter_c4,'-o', 'color', 'b');
xlabel('Matrix size (n)');
ylabel('# of iterations');
title('Iterations require to solve nxn matricies for C = 4');
grid on;
print(fig2, "Jacobi_C4.png", '-dpng'); 

y2 = polyfit(n, iter_c4, 1);
r2 = y2(1) %Power-Law Big O rep of Jacobi complextity growth for c = 4

