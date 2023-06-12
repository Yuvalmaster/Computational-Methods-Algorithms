% First part
data = readtable('Data.xlsx');

x = data.Strain;
y = data.Stress_MPa_;

[a_linear, b_linear, LS_linear] = linear_approx(x, y);
[a_exp, b_exp, LS_exp] = exp_approx(x, y);
[a_log, b_log, LS_log] = log_approx(x, y);
[a_poly_2, b_poly_2, LS_poly_2] = poly_2_approx(x, y);
[a_poly_1, b_poly_1, LS_poly_1] = poly_1_approx(x, y);

LS = [LS_linear, LS_exp, LS_log, LS_poly_2, LS_poly_1];

error = zeros(size(LS, 2), 1);
for i = 1:size(LS, 2)
    error(i) = sum((y - (LS(:, i))).^2);
end

[~, best_approx] = min(error); % index 4 = LS_poly_2 !

%%

plot(x, y, '.k')
hold on
for i = 1:size(LS, 2)
    plot(x, LS(:, i))
end
legend('Real', 'Linear', 'Exponential', 'Logarithmic', 'Polynomial 2', 'Polynomial 1');
xlabel('Strain'); ylabel('Stress')
