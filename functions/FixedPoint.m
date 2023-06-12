function [x, ii] = FixedPoint(f, a, b, epsilon, N_iterations)
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

    error = inf;         % Set initial error
    ii    = 0;           % Iteration number
    x0    = (a+b)/2;     % Initial guess 
    g     = @(x) f(x)+x; % set g function for fixedpoint algorithm
    
    while error > epsilon && ii < N_iterations
        x     = g(x0);
        error = abs(x0-x)/abs(x);
        
        x0    = x; 
        % Increment the iteration counter
        ii    = ii + 1;
    end
    if ii == N_iterations
        disp('Bisection algorithm did not converge within the maximum number of iterations');
    end
end