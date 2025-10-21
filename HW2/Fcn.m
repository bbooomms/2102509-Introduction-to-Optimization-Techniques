function [f, grad, H] = Fcn(x, options)
% Generic function interface for Newton's method
% 
% DEVELOPED BY PHUBET SUWANNO 6530316021 (BOOM)
% PRESENTED TO ASSOC. PROF. DR. SUCHIN ARUNSAWATWONG
% 
% INPUTS:
%   x       : n×1 vector
%   options : 1 = f only, 2 = f + gradient, 3 = f + gradient + Hessian
% OUTPUTS:
%   f       : scalar, function value
%   grad    : n×1 vector, gradient
%   H       : n×n matrix, Hessian

    % Rosenbrog fn
    n = length(x);
    f = 0;
    % gradient is a vector
    if options >= 2, grad = zeros(n,1); else, grad = []; end
    % Hessian is a matrix
    if options == 3, H = zeros(n,n); else, H = []; end

    % Objective
    for i = 1:n-1
        f = f + 100*(x(i+1)-x(i)^2)^2 + (1-x(i))^2;
    end

    % Gradient
    if options >= 2
        for i = 1:n-1
            grad(i)   = grad(i) - 400*x(i)*(x(i+1)-x(i)^2) - 2*(1-x(i));
            grad(i+1) = grad(i+1) + 200*(x(i+1)-x(i)^2);
        end
    end

    % Hessian
    if options == 3
        for i = 1:n-1
            H(i,i)     = H(i,i) + 1200*x(i)^2 - 400*x(i+1) + 2;
            H(i,i+1)   = H(i,i+1) - 400*x(i);
            H(i+1,i)   = H(i+1,i) - 400*x(i);
            H(i+1,i+1) = H(i+1,i+1) + 200;
        end
    end
end