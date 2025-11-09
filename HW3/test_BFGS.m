%% =============================================================
%   TEST SCRIPT : BFGS on Rosenbrock Function
%   Author       : Phubet Suwanno
%   Description  : Run BFGS minimization using two line search settings
% =============================================================

clc; clear;

% ----------------- Problem setup -----------------
x0      = [-1.2; 1];     % standard starting point for Rosenbrock
epsilon = 1e-9;           % convergence tolerance
mu      = 1e-4;           % Armijo parameter
itmax   = 200;            % maximum number of iterations

% ----------------- Case 1 : eta = 0.1 (accurate line search) -----------------
eta1 = 0.1;
fprintf('\n============================================\n');
fprintf('Case 1 : Strong Wolfe Line Search (eta = %.2f)\n', eta1);
fprintf('============================================\n');
[xmin1,fmin1,Xk1,Fk1,Gk1,Lk1,nF1,nG1,IFLAG1] = ...
    BFGS(@Rosenbrock, x0, epsilon, mu, eta1, itmax);

fprintf('Minimum point   : [%f, %f]\n', xmin1(1), xmin1(2));
fprintf('Minimum value   : %.8f\n', fmin1);
fprintf('nF = %d, nG = %d\n', nF1, nG1);
fprintf('Convergence flag: %d\n\n', IFLAG1);

% Print iteration results (first few only)
fprintf('%5s %15s %15s %15s %15s %15s\n', ...
    'Iter', 'x1', 'x2', 'f(x)', 'grad_x1', 'grad_x2');
for k = 1:length(Fk1)
    fprintf('%5d %15.8f %15.8f %15.8f %15.8f %15.8f\n', ...
        k-1, Xk1(1,k), Xk1(2,k), Fk1(k), Gk1(1,k), Gk1(2,k));
end

% ----------------- Case 2 : eta = 0.98 (inaccurate line search) -----------------
eta2 = 0.98;
fprintf('\n============================================\n');
fprintf('Case 2 : Strong Wolfe Line Search (eta = %.2f)\n', eta2);
fprintf('============================================\n');
[xmin2,fmin2,Xk2,Fk2,Gk2,Lk2,nF2,nG2,IFLAG2] = ...
    BFGS(@Rosenbrock, x0, epsilon, mu, eta2, itmax);

fprintf('Minimum point   : [%f, %f]\n', xmin2(1), xmin2(2));
fprintf('Minimum value   : %.8f\n', fmin2);
fprintf('nF = %d, nG = %d\n', nF2, nG2);
fprintf('Convergence flag: %d\n\n', IFLAG2);

fprintf('%5s %15s %15s %15s %15s %15s\n', ...
    'Iter', 'x1', 'x2', 'f(x)', 'grad_x1', 'grad_x2');
for k = 1:length(Fk2)
    fprintf('%5d %15.8f %15.8f %15.8f %15.8f %15.8f\n', ...
        k-1, Xk2(1,k), Xk2(2,k), Fk2(k), Gk2(1,k), Gk2(2,k));
end


% คำนวณ norm ของกราเดียนต์ในแต่ละ iteration
gnorm1 = arrayfun(@(k) norm(Gk1(:,k)), 1:size(Gk1,2));
gnorm2 = arrayfun(@(k) norm(Gk2(:,k)), 1:size(Gk2,2));

iters1 = 0:(numel(gnorm1)-1);
iters2 = 0:(numel(gnorm2)-1);

figure('Name','Convergence: log(norm(grad))','Color','w'); clf;
semilogy(iters1, gnorm1, '-o', 'LineWidth', 1.25, 'MarkerSize', 4); hold on;
semilogy(iters2, gnorm2, '-s', 'LineWidth', 1.25, 'MarkerSize', 4);
grid on; box on;
xlabel('Iteration');
ylabel('||\nabla f(x_k)|| (log scale)');
title('BFGS Convergence on Rosenbrock: log(||\nabla f||)');

% ใส่ legend พร้อมสรุป nF/nG ให้ครบในกราฟ
legend({sprintf('\\eta = 0.10  (nF = %d, nG = %d)', nF1, nG1), ...
        sprintf('\\eta = 0.98  (nF = %d, nG = %d)', nF2, nG2)}, ...
        'Location','northeast');

% บันทึกรูป (เผื่อใส่รายงาน)
set(gcf, 'PaperPositionMode','auto');
print(gcf, 'bfgs_convergence_gradnorm.png', '-dpng', '-r300');
print(gcf, 'bfgs_convergence_gradnorm.pdf', '-dpdf', '-r300');
