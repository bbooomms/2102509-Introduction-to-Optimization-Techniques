clc; clearvars
% Parameters
epsilon = 1e-6; eps_rel = 0.05; eps_abs = 1e-3; itmax = 1000;

% Case 1
disp('--- Case 1: x0 = [-1.2;1] ---');
[xmin1,fmin1,Xk1,Fk1,Gk1,nF1,nG1,nH1,CHN1,IFLAG1] = ...
    Newton(@Fcn,[-1.2;1],epsilon,eps_rel,eps_abs,itmax);

% Case 2
disp('--- Case 2: x0 = [10;10] ---');
[xmin2,fmin2,Xk2,Fk2,Gk2,nF2,nG2,nH2,CHN2,IFLAG2] = ...
    Newton(@Fcn,[10;10],epsilon,eps_rel,eps_abs,itmax);

% Case 3
disp('--- Case 3: x0 = [-1;-1] ---');
[xmin3,fmin3,Xk3,Fk3,Gk3,nF3,nG3,nH3,CHN3,IFLAG3] = ...
    Newton(@Fcn,[-1;-1],epsilon,eps_rel,eps_abs,itmax);

% Print summary with formatted output
fprintf('\n=== Summary of Results ===\n');

fprintf('\nCase 1: x0 = [-1.2; 1]\n');
fprintf('  xmin   = [%.6f, %.6f]\n', xmin1(1), xmin1(2));
fprintf('  fmin   = %.6e\n', fmin1);
fprintf('  Status = %s\n', ternary(IFLAG1==0,'Converged','Not Converged'));

fprintf('\nCase 2: x0 = [10; 10]\n');
fprintf('  xmin   = [%.6f, %.6f]\n', xmin2(1), xmin2(2));
fprintf('  fmin   = %.6e\n', fmin2);
fprintf('  Status = %s\n', ternary(IFLAG2==0,'Converged','Not Converged'));

fprintf('\nCase 3: x0 = [-1; -1]\n');
fprintf('  xmin   = [%.6f, %.6f]\n', xmin3(1), xmin3(2));
fprintf('  fmin   = %.6e\n', fmin3);
fprintf('  Status = %s\n', ternary(IFLAG3==0,'Converged','Not Converged'));



function out = ternary(cond, valTrue, valFalse)
% Simple ternary operator for readability
    if cond
        out = valTrue;
    else
        out = valFalse;
    end
end