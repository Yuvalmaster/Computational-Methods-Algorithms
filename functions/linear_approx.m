function [a, b, LS_linear] = linear_approx(x, y)
    
    n = length(x);
    b = (sum(x.^2)*sum(y) - sum(x.*y)*sum(x))/(n*sum(x.^2)-(sum(x)^2));
    a = (n*sum(x.*y) - sum(x)*sum(y))/(n*sum(x.^2) - (sum(x))^2);

    LS_linear = a.*x + b;

end
