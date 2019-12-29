function [t, y] = DE2_boyadj11(f, t0, tN, y0, y1, h)

    t = t0 : h : tN+h;  % Vector with stepsize h
    N = length(t);      % Initialize length of vector
    y = zeros(1,N);     % Create N long vector of 0s

    y(1) = y0;          % Initial condition y(1)
    y(2) = y0 + y1*h;   % Initial condition y(2)

    for i = 2 : tN/h    % Loop through for N = 2, 3, 4 ...
        
        dy = (y(i) - y(i-1)) / h;
        y(i+1) = (h^2)*f(t(i), y(i), dy) + 2*y(i)- y(i-1);
        
    end
    
end

    
