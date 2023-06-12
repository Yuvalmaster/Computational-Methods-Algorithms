function [p, ii] = Bisection(f, a, b, epsilon, N_iterations)
%   Inputs:
%       a:             min value
%       b:             max value
%       f:             function
%       N_iterations:  Maximum number of iterations
%       epsilon:       Convergence tolerance
%
%   Outputs:
%       p:        Solution
%       ii:       Number of iterations

    error = inf;     % Set initial error
    ii    = 0;       % Iteration number
    p     = (a+b)/2; % Initial bisection
    
    while error > epsilon && ii < N_iterations
        if sign(f(p)) == sign(f(b)) % Choose p based on sign (go left/go right)
            b = p;  
            
            new_p = (a+b)/2;
            error = abs(p-new_p)/abs(new_p); % Calculate convergence criterion
            p     = new_p; 
            % Increment the iteration counter
            ii    = ii+1;
        else
            a = p;
            
            new_p = (a+b)/2;
            error = abs(p-new_p)/abs(new_p); % Calculate convergence criterion
            p     = new_p;
            % Increment the iteration counter
            ii    = ii+1;
        end 
    end
    if ii == N_iterations
        disp('Bisection algorithm did not converge within the maximum number of iterations');
    end
end