function [xmin, fmin, Xk, Fk, Gk, Lk, nF, nG, IFLAG, nReset] = CG(fcn, x0, epsilon, mu, eta, itmax, option)
% ================================================================
% Conjugate Gradient Method (Fletcher–Reeves or Polak–Ribiere)
%
% MODIFIED BY PHUBET SUWANNO
% Syntax:
%   [xmin, fmin, Xk, Fk, Gk, Lk, nF, nG, IFLAG, nReset] = ...
%       CG(fcn, x0, epsilon, mu, eta, itmax, option)
%
% Input Arguments:
%   fcn     : Objective function handle
%             Must support: fcn(x,1) → f(x),   fcn(x,2) → [f(x), ∇f(x)]
%   x0      : Initial point (column vector)
%   epsilon : Convergence tolerance for gradient norm
%   mu, eta : Strong Wolfe condition parameters (0 < mu < eta < 1)
%   itmax   : Maximum number of iterations
%   option  : 1 = Fletcher–Reeves,  2 = Polak–Ribiere+
%
% Output Arguments:
%   xmin    : Approximated minimizer point
%   fmin    : Function value at xmin
%   Xk      : Iterates history (k × dim matrix)
%   Fk      : f(x_k) history
%   Gk      : ∇f(x_k) history
%   Lk      : Step size (lambda) per iteration
%   nF, nG  : Number of function and gradient evaluations per iteration
%   IFLAG   : Exit flag (0 = converged, -999 = max iteration reached)
%   nReset  : Reset indicator (0=no reset, 1=angle reset, 2=loss of descent)
% ================================================================

% storage
Xk = x0';       % iteration 0
Fk = []; Gk = []; Lk = []; nF = []; nG = []; nReset = [];

% Initial evaluation
[fk, gk] = fcn(x0, 2);   
nf_iter = 1; 
ng_iter = 1;
sk = -gk;       % initial steepest descent

i = 0;
while true
    
    % Strong Wolfe line search
    [lambda, nFnew, nGnew] = LineSearch(fcn, x0, sk, 1, mu, eta);
    nf_iter = nf_iter + nFnew;
    ng_iter = ng_iter + nGnew;

    % Step
    x1 = x0 + lambda*sk;
    [fk1, gk1] = fcn(x1, 2);
    nf_iter = nf_iter + 1;
    ng_iter = ng_iter + 1;

    % Record iteration result (k+1)
    Xk = [Xk; x1'];
    Fk = [Fk; fk1];
    Gk = [Gk; gk1'];
    Lk = [Lk; lambda];
    nF = [nF; nf_iter];
    nG = [nG; ng_iter];

    % Check stop
    if norm(gk1) < epsilon
        nReset = [nReset; 0];
        xmin = x1; 
        fmin = fk1; 
        IFLAG = 0;
        return
    end

    % CG beta
    if option == 1
        beta = (gk1'*gk1)/(gk'*gk);           % FR
    else
        beta = max((gk1'*(gk1 - gk))/(gk'*gk), 0); % PR+
    end

    sk1 = -gk1 + beta*sk;

    % Reset conditions
    ct = (sk1'*(-gk1))/(norm(sk1)*norm(gk1));
    ct = max(-1,min(1,ct));   % numerical safety
    angle = acosd(ct);

    if angle > 85
        sk1 = -gk1; rst = 1;    % reset by angle
    elseif sk1'*gk1 >= 0
        sk1 = -gk1; rst = 2;    % lost descent
    else
        rst = 0;
    end
    nReset = [nReset; rst];

    % Prepare next iteration
    x0 = x1;
    fk = fk1;
    gk = gk1;
    sk = sk1;

    % reset counter per iteration
    nf_iter = 0;
    ng_iter = 0;

    i = i + 1;
    if i > itmax
        xmin = x1; fmin = fk1; IFLAG = -999;
        return
    end
end
end
