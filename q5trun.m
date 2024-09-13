clear all
close all
clc

syms x C G M Y real

C = 1;
G = 4;
M = 3;
Y = 0.5;

f = piecewise(abs(x) >= 0.1, 0, ...
              x < 0 & x > -0.1, C * exp(-G * abs(x)) / abs(x)^(1 + Y), ...
              x > 0 & x < 0.1, C * exp(-M * x) / x^(1 + Y));

% Mean
meanExpr = int(x * f, x, -inf, inf);
currentMeanValue = vpa(0.1-meanExpr, 6);

% Truncated Variance
varExpr = int(x^2 * f, x, -inf, inf);
varianceValue = double(varExpr);

% Truncated Skewness
skewExpr = int(x^3 * f, x, -inf, inf);
skewnessValue = double(skewExpr);

% Truncated Kurtosis
kurtExpr = int(x^4 * f, x, -inf, inf);
kurtosisValue = double(kurtExpr);

% Display results in a table
results = table(currentMeanValue, varianceValue, skewnessValue, kurtosisValue, 'VariableNames', {'Mean', 'Truncated_Variance', 'Truncated_Skewness', 'Truncated_Kurtosis'});
disp(results);
