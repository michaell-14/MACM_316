%(C)Michael Lange; MACM 316; 
close all;
clear;
clc;

tic;

%note about parallel:
%overhead management causes issues with timing and delayed program speed

%parallel computing setup
%pool = parpool('local');
%pool = parpool('local', 8); %on my computer, it  8 cores

%Lower n_val according to computer RAM if you get "Out of Memory"
%n_val = [100, 200, 500, 1000, 5000, 10000, 15000, 20000]; % for n by n matrices
n_val = [100, 500, 1000, 5000, 10000, 20000, 30000]; %stress testing
num_run = 20; % Number of runs for averaging

% Pre-allocate arrays for storing average computation times
avg_time_inv = zeros(size(n_val));
avg_time_chol = zeros(size(n_val));

% Part A: Inversion Method

% Warm-up for Inversion Method
%Because of Just in Time (JIT) in matlab
%JIT is where code is compiled only when it is used, not at all at runtime
%meaning that it can cause delays (false reading as the code wasnt require at that point)

n_warmup_inv = n_val(1);
random_diag_warmup_inv = 4 + rand(n_warmup_inv, 1);
random_offdiag_warmup_inv = 0.5 * rand(n_warmup_inv - 1, 1);
A_warmup_inv = sparse(diag(random_offdiag_warmup_inv, 1) + diag(random_offdiag_warmup_inv, -1) + diag(random_diag_warmup_inv));
A_warmup_inv = A_warmup_inv + n_warmup_inv * speye(n_warmup_inv);  % Ensure positive definiteness
b_warmup_inv = rand(n_warmup_inv, 1);
A_warmup_inv = inv(A_warmup_inv); %A\b is better
x = A_warmup_inv * b_warmup_inv;
clear A_warmup_inv b_warmup_inv x_warmup_inv;

%parfor index = 1:length(n_val) %parallel
for index = 1:length(n_val) %regular
    n = n_val(index);
    random_diag = 4 + rand(n, 1);
    random_offdiag = 0.5 * rand(n-1, 1);
    A_n = sparse(diag(random_offdiag,1) + diag(random_offdiag,-1) + diag(random_diag));
    b_n = rand(n, 1);

    total_time = 0;
    for i = 1:num_run
        tic;
        A_inv = inv(A_n); %A\b is better
        x = A_inv * b_n;
        total_time = total_time + toc;
    end
    avg_time_inv(index) = total_time / num_run;
end
clear A_n b_n
% Part B: Cholesky Factorization

% Warm-up run
n_warmup = n_val(1);
random_diag_warmup = 4 + rand(n_warmup, 1);
random_offdiag_warmup = 0.5 * rand(n_warmup-1, 1);
A_warmup = sparse(diag(random_offdiag_warmup,1) + diag(random_offdiag_warmup,-1) + diag(random_diag_warmup));
A_warmup = A_warmup + n_warmup * speye(n_warmup);  % Make the matrix strictly diagonally dominant
b_warmup = rand(n_warmup, 1);
R_warmup = chol(A_warmup, 'lower');  % warm-up Cholesky Factorization
clear A_warm_up b_warm_up R_warm_up random_diag_warm_up random_offdiag_warm_up


%parfor index = 1:length(n_val) %parallel
for index = 1:length(n_val) %regular
    n = n_val(index);
    random_diag = 4 + rand(n, 1);
    random_offdiag = 0.5 * rand(n-1, 1);
    A_n = sparse(diag(random_offdiag,1) + diag(random_offdiag,-1) + diag(random_diag));
    A_n = A_n + n * speye(n);  %matrix strictly diagonal dominant
    b_n = rand(n, 1);

    total_time_c = 0;
    for i = 1:num_run
        tic;
        R = chol(A_n, 'lower'); %lower triangluar 
        x = R' \ (R \ b_n);
        total_time_c = total_time_c + toc;
    end
    avg_time_chol(index) = total_time_c / num_run;
end

% Overlaying plots onto a single plot
fig1 = figure;
loglog(n_val, avg_time_inv, '-o', 'color', 'b', 'DisplayName', 'Inversion Method');
hold on;
loglog(n_val, avg_time_chol, '-s','color', 'r', 'DisplayName', 'Cholesky Factorization');
xlabel('Matrix Size (n)');
ylabel('Average Computing Time (seconds)');
title('Computing Time vs Matrix Size (Log-Log Scale)');
grid on;
legend('show');
print(fig1, "overlay_CA2.png", '-dpng'); 

large_n = n_val >= 1000; %for polyfit

% For Inversion:
p_inv = polyfit(log(n_val(large_n)), log(avg_time_inv(large_n)), 1);
power_inv = p_inv(1);  % The slope is the power
disp(['Inversion slope:', num2str(power_inv)]);

% For Cholesky:
p_chol = polyfit(log(n_val(large_n)), log(avg_time_chol(large_n)), 1); % Linear fit (degree 1)
power_chol = p_chol(1); % The slope is the power
disp(['Cholesky slope:', num2str(power_chol)]);

% These are only for displaying the time differences between the two
fig2 = figure;
subplot(2,1,1);
semilogy(n_val, avg_time_inv, '-o', 'color', 'b');
title('Inversion Method - Log(y-axis) form');
xlabel('Matrix Size (n)')
ylabel('Average Compute Time (Log)')
grid on;

subplot(2,1,2);
semilogy(n_val, avg_time_chol, '-s','color', 'r');
title('Cholesky Factorization - Log(y-axis) form');
xlabel('Matrix Size (n)')
ylabel('Average Compute Time (Log)')
grid on;
print(fig2, "2xSemilogy_CA2.png", '-dpng'); 

program_run_time = toc;
disp("Program run time: ")
disp(program_run_time)
%

%delete(pool) %end parallel comp stuff