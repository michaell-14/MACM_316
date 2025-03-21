% Construct the tridiagonal matrix and right hand side
n = 50; 
h = 1 / (n + 1);
main_diag = (c/h^2) * ones(n,1);
off_diag = (-1/h^2) * ones(n,1);
A_n = spdiags([off_diag main_diag off_diag], [-1 0 1], n, n);
b_n = ones(n, 1);

% Declaratives:
% D_inv = inv(diag(A_n))
% L = tril(A_n, -1) --> lower tri 
% U = triu(A_n, 1) --> upper tri
% c = [2, 4] --> [Part A, Part B]

% Processing:
% x_n+1 = D_inv.*(L+U).*x + D'.*b_n
% Refactoring above: ^^^
% x_new = D_inv * (b_n - (L + U) * x);



