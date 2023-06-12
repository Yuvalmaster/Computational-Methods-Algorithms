function [a, b, LS_log] = log_approx(x, y)

    n = length(x);
    b = (n * sum(y.*log(x)) - sum(y)*sum(log(x)))/(n*sum(log(x).^2) - (sum(log(x)))^2);
    a = (sum(y) - b*sum(log(x)))/n;

    LS_log = a + b*log(x);
    
end
