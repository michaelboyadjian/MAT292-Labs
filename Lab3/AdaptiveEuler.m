function [t, y] = AdaptiveEuler(f, t0, tN, y0, h)

tol = 1e-8;
t = [t0];
y = [y0];

while t(end) + h < tN
    Y = y(end) + h * f(t(end), y(end));
    
    Z = y(end) + h/2 * f(t(end), y(end));
    Z = Z + h/2 * f(t(end) + h/2, Z);
    
    D = Z - Y;
    
    if abs(D) <= tol
        t = [t; t(end) + h];
        y = [y; Y];
    end
    
    h = 0.9 * h * min(max(tol / abs(D), 0.3), 2);
end