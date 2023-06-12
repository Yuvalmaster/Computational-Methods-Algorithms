function [x, iter] = Jacobi(x0, A, b, tol, N_iterations)
%   Jacobi method for solving linear equations Ax = b
%   Input:
%       A:             coefficient matrix
%       b:             right hand side vector
%       x0:            initial guess for solution
%       tol:           tolerance for stopping iteration
%       N_interations: maximum number of iterations
% 
%   Output:
%       x:             solution vector
%       iter:          Number of iterations

    % Check algorithm
    a = abs(A)- abs(diag(A));
    e = eig(A);
    if max(max(a > 0)) == 1 || max(diag(A) == 0) == 1 || max(abs(e)) >= 1 || ~isequal(A,A')
        warning('The conditions for the algorithm are not valid')
    end 

    % initialize variables
    n    = length(b);
    x    = x0;
    err  = Inf; % error between successive iterates
    iter = 0;   % iteration counter
    
    % Jacobi iteration
    while err > tol && iter < N_iterations
        x_new = zeros(n, 1); % new iterate
        for i = 1:n
            % Compute the new value for the i-th unknown
            a = (b(i) - A(i,[1:i-1,i+1:n]) * x([1:i-1,i+1:n]))/A(i, i);
            x_new(i) = a;
        end
        err  = norm(x_new - x, 2)./... % Calculate error according to given criterion
               norm(x_new    , 2);          
        x    = x_new;                  % Update X
        iter = iter + 1;               % Iteration counter
    end
    
    % check for convergence
    if err > tol
        fprintf('Jsacobi method did not converge\n');
    end
end