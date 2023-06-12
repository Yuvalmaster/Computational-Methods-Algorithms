function [x, iter] = HestenesSteifel(x0, A, b, tol, N_iterations)
%   Inputs:
%   A:             Coefficient matrix
%   b:             Right-hand side vector
%   x0:            Initial guess for the solution
%   N_iterations:  Maximum number of iterations
%   tol:           Convergence tolerance
%   w:             relaxation parameter 
%
%   Outputs:
%   x:        Solution vector
%   iter:     Number of iteration until convergance

    % Initialize the iteration counter and the solution vector
    iter = 0;
    x    = x0;
    err  = inf;
    n    = length(b);
    r    = b-A*x;
    p    = r;

 % Check algorithm
    e = eig(A);
    if size(A,1)~= size(A,2) ||  max(abs(e)) >= 1 || ~issymmetric(A)
        warning('The conditions for the algorithm are not valid')
    end
    
% Loop until the maximum number of iterations is reached or convergence is achieved

    while err > tol && iter < N_iterations
        x_prev = x;

        alpha = (transpose(r)*r)/(transpose(p)*A*p); % length of the step
        x = x+alpha*p; % next approximate solution
        r_prev = r;
        r = r-alpha*A*p; % residual
        beta = (transpose(r)*r)/(transpose(r_prev)*r_prev); % relative improvement of this step
        p = r+beta*p; % new search direction            

        % Check for convergence
        err  = norm(x - x_prev, 2)./... % Calculate error according to given criterion
               norm(x         , 2);     
        if err < tol
            break;
        end
        
        % Increment the iteration counter
        iter = iter + 1;
    end
    
    % Check if the maximum number of iterations was reached
    if iter == N_iterations
        disp('Hestenes-Steifel algorithm did not converge within the maximum number of iterations');
    end
end