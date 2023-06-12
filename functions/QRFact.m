function [x, iter] = QRFact(x0, A, b, tol, N_iterations)
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

    [n,m] = size(A);
    
    Q = zeros(n,m);
    R = zeros(m,m);
    Z = [];
    
    %initiation (k=1):
    Z(:,1) = A(:,1);
    R(1,1) = sqrt(transpose(Z(:,1))*Z(:,1));
    Q(:,1) = A(:,1)/R(1,1);
    
    %iterations:
    for k=2:m
        x = 0;
        for i=1:k-1
            R(i,k) = transpose(Q(:,i))*A(:,k);
            x = x + R(i,k)*Q(:,i);
        end
        Z(:,k) = A(:,k) - x; 
        R(k,k) = sqrt(transpose(Z(:,k))*Z(:,k));
        Q(:,k) = Z(:,k)/R(k,k);
    end
    
    %check if A=Q*R, with error er:
    er = 1e-14;
    if abs(A-(Q*R))<er
        disp('QR Factorisation: the condition A = Q*R is valid')
    end
    
    %check if transpose(Q)*Q=I, with error er:
    er = 1e-14;
    if abs(eye(n,m)-(transpose(Q)*Q))<er
        disp('QR Factorisation: the condition transpose(Q)*Q = I is valid')
    end
    
    x = inv(R)*transpose(Q)*b; % Estimated solution
    iter = 1;

end