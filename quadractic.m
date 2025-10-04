function [xmin, fmin, IFLAG, IFunc] = quadractic(a, b, epsilon, itmax)
    % QUADRRACTIC  One-dimensional minimization by Quadratic Interpolation
    %
    % DEVELOPED BY PHUBET SUWANNO 6530316021 (BOOM)
    % PRESENTED TO ASSOC. PROF. DR. SUCHIN ARUNSAWATWONG
    % 2102509 INTRO OPTIMIZATION TECHNIQUE
    %
    %   [xmin, fmin, IFLAG, IFunc] = quadractic(a, b, epsilon, itmax)
    %
    %   This function finds the minimizer of f(lambda) in the interval [a,b]
    %   using the Quadratic Interpolation Method (Powell style).
    %
    %   Input arguments:
    %     a       - left endpoint of the initial interval
    %     b       - right endpoint of the initial interval
    %     epsilon - stopping tolerance (see eqn (5.2) in lecture notes)
    %     itmax   - maximum number of iterations allowed
    %
    %   Output arguments:
    %     lmin    - estimated minimizer of f
    %     fmin    - function value f(lmin)
    %     IFLAG   - flag (0 if successful, -999 otherwise)
    %     IFunc   - number of function evaluations
    %
    %   Notes:
    %   - The target function f(lambda) must be defined separately in f.m
    %     Example:
    %         function y = f(lambda)
    %             y = lambda.^5 - 5*lambda.^3 - 20*lambda + 5;
    %         end
    %   - Method uses 3 points (A,B,C), fits a parabola h(lambda),
    %     computes vertex lambda*, and updates according to 4 cases.
    %   - Stop criterion: |h(lambda*) - f(lambda*)| / |f(lambda*)| <= epsilon
    %
    %   Example of usage:
    %     [lmin,fmin,IFLAG,IFunc] = quadractic(-5,5,1e-6,100)
  
    % printing format (>= 4 digits; set to 10 for clarity)
    DIG  = 10; % number of digits to print
    fmt  = ['%.', num2str(max(4,DIG)), 'f'];   % e.g., '%.10f'

    A = a; C = b; B = 0.5*(a+b);
    fA = f(A); fB = f(B); fC = f(C);
    IFunc = 3; IFLAG = -999; 
    k = 0;

    % header print
    fprintf('Iter  |      A            B            C       |    lambda*        f(lambda*)         IFunc\n');
    fprintf('------+------------------------------------------------------------------------------------\n');

    while k < itmax
        % use eq (5.1) in EE509 handout to find λ*
        num = fA*(B^2 - C^2) + fB*(C^2 - A^2) + fC*(A^2 - B^2);
        den = 2*( fA*(B - C) + fB*(C - A) + fC*(A - B) );

        % To prevent the case where the denominator = 0, which would cause NaN and an error
        if abs(den) < 1e-14 
            break; 
        end

        lambda_star = num/den;
        fstar = f(lambda_star); 
        IFunc = IFunc+1;

        % define a Parabola h(λ) = a + bλ + cλ²
        denom = (A-B)*(B-C)*(C-A);
        a_coef = (fA*B*C*(C-B) + fB*C*A*(A-C) + fC*A*B*(B-A)) / denom;
        b_coef = (fA*(B^2-C^2) + fB*(C^2-A^2) + fC*(A^2-B^2)) / denom;
        c_coef = -(fA*(B-C) + fB*(C-A) + fC*(A-B)) / denom;

        hstar = a_coef + b_coef*lambda_star + c_coef*lambda_star^2;

        % print this iteration
        fprintf(['%4d | ', fmt, '  ', fmt, '  ', fmt, '  |  ', fmt, '   ',   fmt, '    %5d\n'], ...
                k, A, B, C, lambda_star, fstar, IFunc);

        % use eq (5.2) in EE509 Handout to find a Stop criterion
        if fstar ~= 0
            if abs(hstar - fstar)/abs(fstar) <= epsilon
                xmin = lambda_star; 
                fmin = fstar; 
                IFLAG = 0; 
                return;
            end
        end

        % Update A,B,C follow by Cases
        if (lambda_star < B) && (fstar < fB)
            % Case 1: λ* lies to the left of B and better than B
            C = B; fC = fB; 
            B = lambda_star; fB = fstar;
        
        elseif (lambda_star < B) && (fstar >= fB)
            % Case 2: λ* lies to the left of B and worse than B
            A = lambda_star; fA = fstar;
        
        elseif (lambda_star > B) && (fstar < fB)
            % Case 3: λ* lies to the right of B and better than B
            A = B; fA = fB; 
            B = lambda_star; fB = fstar;
        
        else
            % Case 4: λ* lies to the right of B and worse than B
            C = lambda_star; fC = fstar;
        end
        k = k+1;
    end

    % if not converge
    [fmin, idx] = min([fA, fB, fC]);
    L = [A, B, C];     % save in array
    xmin = L(idx);

    % final line print for visibility
    fprintf('----> Stop by itmax or degenerate parabola. Return best of {A,B,C}.\n');
    fprintf(['Best lmin = ', fmt, ', fmin = ', fmt, ', IFLAG = %d, IFunc = %d\n'], xmin, fmin, IFLAG, IFunc);

end