function [f, gradient] = Rosenbrock(x, options)
f = 100*(x(2)-x(1)^2)^2 + (1 - x(1))^2;
if options == 2
    gradient = [-400*x(1)*(x(2)-x(1)^2) - 2*(1 - x(1));
                 200*(x(2) - x(1)^2)];
else
    gradient = [];
end
end
