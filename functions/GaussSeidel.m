function [x, iter] = GaussSeidel(x0, A, b, tol, N_iterations)
%   Inputs:
%       A:             Coefficient matrix
%       b:             Right-hand side vector
%       x0:            Initial guess for the solution
%       N_iterations:  Maximum number of iterations
%       tol:           Convergence tolerance
%
%   Outputs:
%       x:        Solution vector
%       iter:     Number of iteration until convergance

    % Initialize the iteration counter and the solution vector
    iter = 0;
    x    = x0;
    err  = inf;
    
    % Check algorithm
    a = abs(A)- abs(diag(A));
    if max(max(a > 0)) || size(A,1)~= size(A,2) || det(A) == 0 || rank([A b]) ~= rank(A)
        warning('The conditions for the algorithm are not valid')
    end
    
    % Loop until the maximum number of iterations is reached or convergence is achieved
    while err > tol && iter < N_iterations
        % Save the previous solution
        x_prev = x;
        
        % Loop through each equation
        for i = 1:length(b)
            % Compute the new value for the i-th unknown
            x(i) = (b(i) - A(i,1:i-1)*x(1:i-1) - A(i,i+1:end)*x(i+1:end)) / A(i,i);
        end
        
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
        disp('Gauss-Seidel algorithm did not converge within the maximum number of iterations');
    end

end