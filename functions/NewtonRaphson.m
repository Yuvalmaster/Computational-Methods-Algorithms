function [x, ii] = NewtonRaphson(f, a, b, epsilon, N_iterations)
%   Inputs:
%       a:             min value
%       b:             max value
%       f:             function
%       N_iterations:  Maximum number of iterations
%       epsilon:       Convergence tolerance
%
%   Outputs:
%       x:        Solution
%       ii:       Number of iterations  

    ii = 0;       % Iteration number
    x0 = (a+b)/2; % Initial guess
    xi = x0;      % Set initial movement

    f_der = matlabFunction(diff(sym(f),1)); % Derivative of f
    f_der_zero = f_der(x0);                 % To improve run time the function will not iterate the x0

    while abs(f(xi)) > epsilon && ii < N_iterations
        
        x  = xi - f(xi)/f_der_zero;
        xi = x;
        % Increment the iteration counter
        ii = ii + 1;
    end
end