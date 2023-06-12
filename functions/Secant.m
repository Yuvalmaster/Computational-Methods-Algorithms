function [q, ii] = Secant(f, a, b, epsilon, N_iterations)
%   Inputs:
%   a:             min value
%   b:             max value
%   f:             function
%   N_iterations:  Maximum number of iterations
%   epsilon:       Convergence tolerance
%
%   Outputs:
%   q:        Solution
%   ii:       Number of iterations   
    
    error = inf; % Set initial error
    ii    = 0;   % Iteration number
    xa    = a;   % Starting X_a point
    xb    = b;   % Starting X_b point
    q = xa - f(xa)*(xb-xa)/(f(xb)-f(xa));
    
    while error > epsilon && ii < N_iterations
        if sign(f(q)) == sign(f(xb))
            xb = q;
            
            new_q = xa -f(xa)*(xb-xa)/(f(xb)-f(xa));
            error = abs(q-new_q)/abs(new_q);
            q     = new_q;
            % Increment the iteration counter
            ii    = ii+1;
        else
            xa = q;
            
            new_q = xa -f(xa)*(xb-xa)/(f(xb)-f(xa));
            error = abs(q-new_q)/abs(new_q);
            q     = new_q;
            % Increment the iteration counter
            ii    = ii+1;
        end 
    end
end