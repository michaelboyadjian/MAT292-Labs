function [t, y] = heun(f, t0, tN, y0, h)

% Create vector between t0 and tn with step size h
t = t0 : h : tN; 
N = length(t); 
y = NaN(N, 1);

% Initial condition
y(1) = y0;

% Iteration through for loop
for i = 1 : N - 1
    g = y(i) + h * f(t(i), y(i));
    y(i + 1) = y(i) + h / 2 * (f(t(i), y(i)) + f(t(i + 1), g));
end