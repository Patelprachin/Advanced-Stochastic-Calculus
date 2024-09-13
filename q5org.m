clear all
close all
clc

syms x C G M Y real
C = 1;
G = 4;
M = 3;
Y = 0.5;

f = piecewise(x < 0, C .* exp(-G * abs(x)) / abs(x)^(1 + Y), x > 0, C * exp(-M * x) / x^(1 + Y));

varianceExpr = double(int(x^2 * f, x, -inf, inf));
skewExpr = double(int(x^3 * f, x, -inf, inf));
kurtExpr = double(int(x^4 * f, x, -inf, inf));
meanExpr = double(int(x * f, x, -inf, inf));

results = table(meanExpr, varianceExpr, skewExpr, kurtExpr, 'VariableNames', {'Mean', 'Truncated_Variance', 'Truncated_Skewness', 'Truncated_Kurtosis'});
disp(results);

