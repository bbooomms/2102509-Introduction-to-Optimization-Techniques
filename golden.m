function [xmin, fmin, IFLAG, IFunc] = golden(a, b, epsilon, itmax)
% GOLDEN  One-dimensional minimization by golden-section search (derivative-free).
%
% DEVELOPED BY PHUBET SUWANNO 6530316021 (BOOM)
% PRESENTED TO ASSOC. PROF. DR. SUCHIN ARUNSAWATWONG
%
%   [xmin, fmin, IFLAG, IFunc] = GOLDEN(f, a, b, epsilon, itmax)
%
%   Finds an approximate minimizer of the real-valued scalar function f on
%   the closed interval [a, b] using the golden-section search. This method
%   is derivative-free and only requires function evaluations at scalar points.
%
%   INPUTS
%     f        - function handle to a real-valued scalar function f(x).
%                The algorithm will evaluate f only at scalar x.
%     a, b     - scalars defining the initial interval with a < b.
%     epsilon  - positive scalar tolerance for the interval width; the loop
%                stops when |b - a| <= epsilon.
%     itmax    - positive integer, maximum number of iterations allowed.
%
%   OUTPUTS
%     xmin     - estimated minimizer (midpoint of the final interval).
%     fmin     - f(xmin), function value at the estimated minimizer.
%     IFLAG    - termination flag:
%                  0 : converged (|b-a| <= epsilon),
%                  1 : reached maximum iterations without meeting tolerance,
%                -999: invalid input (e.g., a >= b, epsilon <= 0, or itmax <= 0).
%     IFunc    - number of iterations actually performed.
%
%   ASSUMPTIONS & NOTES
%   - f should be unimodal on [a, b] (i.e., has a single local minimum).
%     If this assumption is violated, the method may converge to a non-minimal point
%     or stall near a boundary.
%   - This implementation prints a per-iteration debug line showing (a, b, x1, x2,
%     f(x1), f(x2)). To suppress console output, comment out the fprintf lines
%     in the loop (marked as "print this iteration").
%   - No derivatives are used. Only function values are required.
%
%   COMPLEXITY
%   - The interval length shrinks by a constant factor (~0.618) per iteration.
%     Roughly, the number of iterations needed is
%         N ≈ ceil( log(epsilon / (b0 - a0)) / log(0.618...) )
%     where [a0, b0] is the initial interval.
%
%   EXAMPLES
%     % Using a function file f.m:
%     %   function y = f(x), y = x.^5 - 5*x.^3 - 20*x + 5; end
%     [xmin, fmin, flag, k] = golden(@f, 0, 3, 1e-6, 2000);
%
%     % Using an anonymous function:
%     g = @(x) (x-2).^2 + 1;
%     [xmin, fmin, flag, k] = golden(g, 0, 5, 1e-8, 1000);
%
%   EDGE CASES
%   - If any input is invalid (a >= b, epsilon <= 0, itmax <= 0),
%     the function returns NaN outputs and IFLAG = -999.
%   - If the tolerance is very small relative to machine precision, the loop may
%     terminate by itmax; check IFLAG.

% Golden Section Search for 1D minimization
    % IFLAG = 0    → converged
    % IFLAG = 1    → reached max iterations
    % IFLAG = -999 → error (i.e. input are not valid)

    % --- เช็ค input ก่อน ---
    if a >= b || epsilon <= 0 || itmax <= 0
        xmin = NaN; fmin = NaN;
        IFLAG = -999;
        IFunc = 0;
        return;
    end

    golden_ratio = (sqrt(5)-1)/2;
    x1 = b - golden_ratio * (b - a);
    x2 = a + golden_ratio * (b - a);

    k = 0;
    % header print
    fprintf('Iter |        a             b         |         x1            x2        |     f(x1)          f(x2)\n');
    fprintf('-------------------------------------------------------------------------------------------------\n');

    while (abs(b-a) > epsilon && k < itmax)
        k = k + 1;

        if f(x2) > f(x1)
            b = x2;
            x2 = x1;
            x1 = b - golden_ratio * (b - a);
        else
            a = x1;
            x1 = x2;
            x2 = a + golden_ratio * (b - a);
        end
        % print this iteration
        fprintf('%4d | %13.6e  %13.6e | %13.6e  %13.6e | %13.6e  %13.6e\n', ...
        k, a, b, x1, x2, f(x1), f(x2));
    end

    xmin = (a+b)/2;
    fmin = f(xmin);
    IFunc = k;

    if abs(b-a) <= epsilon
        IFLAG = 0;  % converged
    else
        IFLAG = 1;  % reached max iterations
    end

end