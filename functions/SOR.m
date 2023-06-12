function [x, iter] = SOR(x0, A, b, tol, N_iterations, w)
%   Inputs:
%       A:             Coefficient matrix
%       b:             Right-hand side vector
%       x0:            Initial guess for the solution
%       N_iterations:  Maximum number of iterations
%       tol:           Convergence tolerance
%       w:             relaxation parameter 
%
%   Outputs:
%       x:        Solution vector
%       iter:     Number of iteration until convergance

    % Initialize the iteration counter and the solution vector
    iter = 0;
    x    = x0;
    err  = inf;
    n    = length(b);

 % Check algorithm
    a = abs(A)- abs(diag(A));
    if max(max(a > 0)) || size(A,1)~= size(A,2) || w >= 2 || w <= 0
        warning('The conditions for the algorithm are not valid')
    end


    % Loop until the maximum number of iterations is reached or convergence is achieved
    while err > tol && iter < N_iterations
        x_prev = x;
    
        for i = 1:n
            sigma = 0;
            for j = 1:i-1
                sigma = sigma + A(i,j)*x(j);
            end
            for j = i+1:n
                sigma = sigma + A(i,j)*x_prev(j);
            end
            x(i) = (1-w)*x_prev(i) + w*(b(i) - sigma)/A(i,i);
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
        disp('SOR algorithm did not converge within the maximum number of iterations');
    end

end