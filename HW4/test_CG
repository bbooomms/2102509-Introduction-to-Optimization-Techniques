clear; clc;

% Initial setting
x0      = [-1.2; 1];
epsilon = 1e-6;
mu      = 1e-4;
eta     = 0.9;
itmax   = 500;

% ================= Fletcher–Reeves =================
[xmin1, fmin1, Xk1, Fk1, Gk1, Lk1, nF1, nG1, IFLAG1, rst1] = ...
    CG(@Rosenbrock, x0, epsilon, mu, eta, itmax, 1);

if length(rst1) < length(Fk1)
    rst1 = [rst1; zeros(length(Fk1)-length(rst1),1)];
end

fprintf("\n==== Fletcher-Reeves ====\n");
fprintf("%-6s %-20s %-12s %-12s %-10s %-6s\n", ...
        "Iter","x","f(x)","||g||","lambda","rst");
fprintf('%s\n', repmat('-',1,80));


for k = 1:length(Fk1)
    fprintf("%-6d [%.6f %.6f]  %.4e   %.3e   %.3e   %d\n", ...
        k-1, Xk1(k,1), Xk1(k,2), Fk1(k), norm(Gk1(k,:)), Lk1(k), rst1(k));
end
fprintf("x* = [%.6f %.6f],  f(x*) = %.6e,  IFLAG = %d\n\n", xmin1(1), xmin1(2), fmin1, IFLAG1);


% ================= Polak–Ribiere+ =================
[xmin2, fmin2, Xk2, Fk2, Gk2, Lk2, nF2, nG2, IFLAG2, rst2] = ...
    CG(@Rosenbrock, x0, epsilon, mu, eta, itmax, 2);

if length(rst2) < length(Fk2)
    rst2 = [rst2; zeros(length(Fk2)-length(rst2),1)];
end

fprintf("\n==== Polak-Ribiere+ ====\n");
fprintf("%-6s %-20s %-12s %-12s %-10s %-6s\n", ...
        "Iter","x","f(x)","||g||","lambda","rst");
fprintf('%s\n', repmat('-',1,80));

for k = 1:length(Fk2)
    fprintf("%-6d [%.6f %.6f]  %.4e   %.3e   %.3e   %d\n", ...
        k-1, Xk2(k,1), Xk2(k,2), Fk2(k), norm(Gk2(k,:)), Lk2(k), rst2(k));
end
fprintf("x* = [%.6f %.6f],  f(x*) = %.6e,  IFLAG = %d\n\n", xmin2(1), xmin2(2), fmin2, IFLAG2);
