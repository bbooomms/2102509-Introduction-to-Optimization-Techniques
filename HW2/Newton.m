function [xmin,fmin,Xk,Fk,Gk,nF,nG,nH,CHN,IFLAG] = Newton(Fcn,x0,epsilon,eps_rel,eps_abs,itmax)
% NEWTON  Modified Newton's Method with golden-section line search
%
% DEVELOPED BY PHUBET SUWANNO 6530316021 (BOOM)
% PRESENTED TO ASSOC. PROF. DR. SUCHIN ARUNSAWATWONG
%
%   [xmin,fmin,Xk,Fk,Gk,nF,nG,nH,CHN,IFLAG] = NEWTON(Fcn,x0,epsilon,eps_rel,eps_abs,itmax)
%
%   Minimizes an n-dimensional function using a modified Newton's method.
%   If the Newton direction is not a descent direction, it is replaced
%   by the steepest descent direction.
%   Step lengths are chosen using Boom's golden.m line search.
%
%   INPUTS:
%     Fcn     - function handle: [f,g,H] = Fcn(x,options)
%               options = 1 → f only
%               options = 2 → f and gradient
%               options = 3 → f, gradient, Hessian
%     x0      - starting point (column or row vector)
%     epsilon - tolerance for stopping (based on norm of gradient)
%     eps_rel - relative tolerance for line search (passed to golden)
%     eps_abs - absolute tolerance for line search (passed to golden)
%     itmax   - maximum number of iterations
%
%   OUTPUTS:
%     xmin    - estimated minimizer
%     fmin    - function value at xmin
%     Xk      - array of iterates (rows are x_k)
%     Fk      - values f(x_k)
%     Gk      - norms of gradients ||∇f(x_k)||
%     nF,nG,nH- counters: #evals of f, gradient, Hessian
%     CHN     - 0 if Newton step used, 1 if steepest descent used
%     IFLAG   - termination flag:
%                 0 = converged (||∇f|| <= epsilon)
%                 1 = reached max iterations
%
%   NOTES:
%   - Requires golden.m (Boom's implementation of golden-section search).
%   - Designed for smooth functions with continuous Hessian.
%   - This is derivative-based, but the line search is derivative-free.

    % Initialize
    xk = x0(:);
    [f,g,H] = Fcn(xk,3);

    % History containers
    Xk = xk'; Fk = f; Gk = norm(g);
    nF = 1; nG = 1; nH = 1;
    CHN = NaN;  % sigma undefined at k=0

    k = 0; IFLAG = 1;

    % Print header
    fprintf('\n%-4s | %-20s | %-12s | %-12s | %-6s | %-6s | %-6s | %-6s\n', ...
        'k','x_k','f(x_k)','||grad||','nF','nG','nH','sigma');
    fprintf(repmat('-',1,95)); fprintf('\n');

    % Print row k=0
    fprintf('%4d | %-20s | %12.6e | %12.6e | %6d | %6d | %6d | %6s\n', ...
        k, mat2str(xk',4), f, norm(g), nF, nG, nH, '-');

    % Main loop
    while (norm(g) > epsilon) && (k < itmax)
        k = k + 1;

        % Newton direction
        sk = -H\g;

        % If not descent → steepest descent
        if g'*sk >= 0
            sk = -g;
            sigma = 1;
        else
            sigma = 0;
        end
        CHN = [CHN; sigma];

        % Line search function
        phi = @(alpha) Fcn(xk + alpha*sk, 1);

        % Call Boom's golden2
        [xmin_golden, ~, ~, ~] = golden2(phi, 0, 1, eps_abs, 100);
        lambda = xmin_golden;

        % Update iterate
        xk = xk + lambda*sk;
        [f,g,H] = Fcn(xk,3);

        % Save history
        Xk = [Xk; xk'];
        Fk = [Fk; f];
        Gk = [Gk; norm(g)];
        nF = [nF; nF(end)+1];
        nG = [nG; nG(end)+1];
        nH = [nH; nH(end)+1];

        % Print iteration info
        fprintf('%4d | %-20s | %12.6e | %12.6e | %6d | %6d | %6d | %6d\n', ...
            k, mat2str(xk',4), f, norm(g), nF(end), nG(end), nH(end), sigma);
    end

    % Final outputs
    xmin = xk; fmin = f;
    if norm(g) <= epsilon
        IFLAG = 0;
    end
end