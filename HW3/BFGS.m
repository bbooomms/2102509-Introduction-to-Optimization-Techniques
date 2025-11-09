function [xmin,fmin,Xk,Fk,Gk,Lk,nF,nG,IFLAG] = BFGS(FcnName,x0,epsilon,mu,eta,itmax)

% -------------------------------------------------------------
% BFGS METHOD with STRONG WOLFE LINE SEARCH
% -------------------------------------------------------------
% Syntax:
%   [xmin, fmin, Xk, Fk, Gk, Lk, nF, nG, IFLAG] = ...
%       BFGS(FcnName, x0, epsilon, mu, eta, itmax)
%
% Description:
%   Performs unconstrained minimization using the BFGS quasi-Newton method.
%   A line search satisfying the Strong Wolfe conditions is used to determine
%   each step length.
%
% Input Arguments:
%   FcnName  – function handle of the form [f, g] = FcnName(x, mode)
%               mode = 1 → returns f only
%               mode = 2 → returns both f and gradient g
%   x0       – initial guess (column vector)
%   epsilon  – stopping tolerance for convergence (e.g. norm(∇f) < epsilon)
%   mu, eta  – parameters for Armijo and Strong Wolfe conditions (0 < mu < eta < 1)
%   itmax    – maximum number of iterations allowed
%
% Output Arguments:
%   xmin     – estimated minimizer
%   fmin     – function value at xmin
%   Xk       – matrix containing all iterates x_k
%   Fk       – vector of objective values f(x_k)
%   Gk       – matrix of gradients ∇f(x_k)
%   Lk       – vector of step lengths λ_k
%   nF, nG   – counts of function and gradient evaluations
%   IFLAG    – termination flag:
%                 0     → convergence achieved
%                -999   → maximum iteration reached
% -------------------------------------------------------------

Xk = []; Fk = []; Gk = []; Lk = [];
nF = 0; nG = 0;
IFLAG = -999;
B = eye(2);

i = 1;
while i <= itmax

    % evaluate f, g
    [f0, g0] = FcnName(x0, 2); 
    nF = nF + 1; nG = nG + 1;

    % initial step
    a = 1; 
    s = B\(-g0); 

    % line search
    [lambda,nFnew,nGnew] = LineSearch(FcnName, x0, s, a, mu, eta);
    nF = nF + nFnew; nG = nG + nGnew;

    % record data
    Xk(:,i) = x0; 
    Fk(i) = f0; 
    Gk(:,i) = g0; 
    Lk(i) = lambda;

    % update B
    x1 = x0 + lambda*s;
    [f1,g1] = FcnName(x1, 2); 
    nF = nF + 1; nG = nG + 1;
    delta_g = g1 - g0;
    delta_x = lambda*s;
    B = B + delta_g*delta_g'/dot(delta_g,delta_x) - B*(delta_x*delta_x')*B/(delta_x'*B*delta_x);

    % terminate
    if norm(g1) < epsilon
        xmin = x1; fmin = f1; IFLAG = 0;
        disp('search successful.');
        break
    elseif i == itmax
        xmin = 0; fmin = 0; IFLAG = -999;
        disp('search unsuccessful.');
        break
    end

    % update
    x0 = x1; f0 = f1; g0 = g1;
    i = i + 1;
end

end
