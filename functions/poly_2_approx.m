function [a, b, LS_poly_2] = poly_2_approx(x, y)
    
    n = length(x);
    B = [ones(n, 1), x, x.^2];
    A = transpose(B) * B;

    b = transpose(B) * y;
    a = inv(A) * b;

    LS_poly_2 = a(3)*(x.^2) + a(2)*x + a(1); 

end
