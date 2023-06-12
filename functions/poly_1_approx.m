function [a, b, LS_poly_1] = poly_1_approx(x, y)

    n = length(x);
    B = [ones(n, 1), x];
    A = transpose(B) * B;

    b = transpose(B) * y;
    a = inv(A) * b;

    LS_poly_1 = a(2)*x + a(1); 

end
