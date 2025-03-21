close all;
clear all;

% Processing:
% x_n+1 = D_inv.*(L+U).*x + D'.*b_n
% Refactoring above: ^^^
% x_new = D_inv * (b_n - (L + U) * x);

tol = 1e-4;
c = [2, 4];  %[Part A, Part B]
n = [10, 50, 100, 500];
%n = [50, 100, 500, 1000, 5000]; %early test
%n = [100, 500, 1000, 5000, 10000, 20000, 30000]; %stress testing
max_iter = 10000;  % Safety limit for iterations

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

        % Declaratives:
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

        % Jacobi Iteration
        while iter < max_iter
            x_new = D_inv * (b_n - (L + U) * x); % Jacobi update
            rel_error = norm(x_new - x_exact, 2) / norm(x_exact, 2);  % Relative error
            iter = iter + 1;

            % Check stopping condition
            if rel_error < tol
                disp(['c = ', num2str(c(i)), ', n = ', num2str(n(j)), ', Iterations: ', num2str(iter)]);
                disp(" ")
                break
            end
            x = x_new; % Update x
        end
        if iter == max_iter
            warning(['Jacobi method did not converge for c = ', num2str(c(i)), ', n = ', num2str(n(j))]);
        end
    end
end