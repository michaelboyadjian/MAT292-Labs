function [t, x] = solvesystem_boyadj11(f1, f2, t0, tN, x0, h)

% Initialize vector with step size h
t = t0:h:tN+h;

% Initialize counter index
j = 1;

% Initial conditions
y1 = [x0(1)];
y2 = [x0(2)];

% For loop to solve ODE
for i = t0:h:tN
    
    % Slopes for IEM and solving for y1
    m0 = f1(i,y1(j), y2(j));
    m1 = f2(i,y1(j), y2(j));
    m2 = f1(i+h,y1(j)+h*m0, y2(j)+m1*h);
    
    y1(j+1) = y1(j)+h*(m0+m2)/2;
    
    % Solving for y2
    m3 = f2(i+h, y1(j)+m0*h, y2(j)+m1*h);
    y2(j+1) = y2(j)+h*(m1+m3)/2;
    
    % Increment counter
    j = j+1;
    
end

% Vectorization of solutions
x=[y1;y2];

end
