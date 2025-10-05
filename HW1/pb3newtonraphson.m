clc; clear

f  = @(x) x.^5 - 5*x.^3 - 20*x + 5; % define function
df = @(x) 5*x.^4 - 15*x.^2 - 20; % derivative

x = 0; % initial guess
tol = 1e-10;              
maxIter = 100; % ป้องกัน Infinite Loop
k = 0;
xnew = x; % start same
dx   = Inf; % ใหญ่ ๆ เพื่อเข้า loop ได้
fprintf('iter|         x*          |          error\n')
fprintf('-----------------------------------------------------\n')
while abs(dx) >= tol && k < maxIter
    dx   = -f(x)/df(x);   % Newton step
    xnew = x + dx; 
    k = k + 1;
    x = xnew; % update current
    fprintf('  %d | x* ≈ %.10f | |x%d-x%d| = %.10f \n', k, xnew, k, k-1, abs(dx));
end
