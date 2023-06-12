function [a, b, LS_exp] = exp_approx(x, y)

    log_a = (sum(x.^2)*sum(log(y)) - sum(x.*log(y))*sum(x))/(length(x)*sum(x.^2)-(sum(x)^2));
    b = (length(x)*sum(x.*log(y)) - sum(x)*sum(log(y)))/(length(x)*sum(x.^2) - (sum(x))^2);
    
    a = exp(log_a);

    log_LS_exp = log_a + b*x;

    LS_exp = exp(log_LS_exp);

end
