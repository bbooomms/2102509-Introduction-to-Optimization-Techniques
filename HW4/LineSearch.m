function [lambda, nF, nG] = LineSearch(FcnName, x0, direction, alpha_init, mu, eta)
% =============================================================
% LineSearch — Strong Wolfe Line Search using Backtracking & Bisection
% =============================================================
% Syntax:
%   [lambda, nF, nG] = LineSearch(FcnName, x0, direction, alpha_init, mu, eta)
%
% Description:
%   Determines a step length λ that satisfies the Strong Wolfe conditions.
%   The algorithm first expands the step size to bracket a minimum, then
%   applies a bisection strategy within that interval.
%
% Input Arguments:
%   FcnName     — Function handle returning f and gradient g
%                 (mode 1 → f only, mode 2 → [f, g])
%   x0          — Current iterate
%   direction   — Search direction
%   alpha_init  — Initial step size
%   mu, eta     — Armijo and Strong Wolfe parameters (0 < mu < eta < 1)
%
% Output Arguments:
%   lambda      — Step length satisfying Strong Wolfe conditions
%   nF, nG      — Number of function and gradient evaluations
% =============================================================

nF = 0; nG = 0;                       % counters
alpha = alpha_init;                   % current step length guess
[f_curr, grad_curr] = FcnName(x0, 2); % evaluate at current point
nF = nF + 1; nG = nG + 1;

% store initial values
f_prev = f_curr;
alpha_prev = 0;

% initialize bracket
alpha_low = NaN; 
alpha_high = NaN;

% -------------------------------------------------------------
% Phase 1: expand alpha until a bracket is found
% -------------------------------------------------------------
while true
    [f_new, grad_new] = FcnName(x0 + alpha * direction, 2);
    nF = nF + 1; nG = nG + 1;

    % Armijo or f increase → bracket found
    if f_new > f_curr + mu * alpha * dot(direction, grad_curr) || f_new >= f_prev
        alpha_low  = alpha_prev;
        alpha_high = alpha;
        break

    % curvature condition satisfied → done
    elseif abs(dot(grad_new, direction)) <= -eta * dot(grad_curr, direction)
        lambda = alpha;
        return

    % gradient sign changed → bracket found
    elseif dot(grad_new, direction) >= 0
        alpha_low  = alpha_prev;
        alpha_high = alpha;
        break
    end

    % expand step length
    f_prev   = f_new;
    alpha_prev = alpha;
    alpha = 2 * alpha;
end

% -------------------------------------------------------------
% Phase 2: zoom (bisection) within [alpha_low, alpha_high]
% -------------------------------------------------------------
[f_low, ~] = FcnName(x0 + alpha_low * direction, 1);
nF = nF + 1;

while true
    alpha = 0.5 * (alpha_low + alpha_high);
    [f_new, grad_new] = FcnName(x0 + alpha * direction, 2);
    nF = nF + 1; nG = nG + 1;

    % Armijo or f still high → move upper bound down
    if f_new > f_curr + mu * alpha * dot(direction, grad_curr) || f_new > f_low
        alpha_high = alpha;

    else
        % strong Wolfe satisfied → return
        if abs(dot(grad_new, direction)) <= -eta * dot(grad_curr, direction)
            lambda = alpha;
            return
        end

        % keep bracketing interval valid
        if dot(grad_new, direction) * (alpha_high - alpha_low) >= 0
            alpha_high = alpha_low;
        end

        alpha_low = alpha;
        f_low = f_new; % update f(alpha_low)
    end
end

end
