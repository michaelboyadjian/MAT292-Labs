%% ODE Lab: Creating your own ODE solver in MATLAB
%
% In this lab, you will write your own ODE solver for the Improved Euler 
% method (also known as the Heun method), and compare its results to those 
% of |ode45|.
%
% You will also learn how to write a function in a separate m-file and 
% execute it.
% 
% Opening the m-file lab3.m in the MATLAB editor, step through each
% part using cell mode to see the results.  Compare the output with the
% PDF, which was generated from this m-file.
%
% There are six (6) exercises in this lab that are to be handed in on the
% due date. Write your solutions in the template, including
% appropriate descriptions in each step. Save the .m files and submit them 
% online on Quercus.
%
% MAT292, Fall 2019, Stinchcombe & Parsch, modified from
% MAT292, Fall 2018, Stinchcombe & Khovanskii, modified from
% MAT292, Fall 2017, Stinchcombe & Sinnamon, modified from
% MAT292, Fall 2015, Sousa, modified from
% MAT292, Fall 2013, Sinnamon & Sousa, modified from
% MAT292, Fall 2011, Hart & Pym

%% Student Information
%
% Student Name: Michael Boyadjian
%
% Student Number: 1005109142
%

%% Creating new functions using m-files.
%  
% Create a new function in a separate m-file:
%
% Specifics:  Create a text file with the file name f.m
% with the following lines of code (text):
%
%  function y = f(a,b,c) 
%  y = a+b+c;
%
% Now MATLAB can call the new function f (which simply accepts 3 numbers
% and adds them together).  
% To see how this works, type the following in the matlab command window:
% sum = f(1,2,3)

%% Exercise 1
%
% Objective: Write your own ODE solver (using the Heun/Improved Euler
% Method).
%
% Details: This m-file should be a function which accepts as variables 
% (t0,tN,y0,h), where t0 and tN are the start and end points of the 
% interval on which to solve the ODE, y0 is the initial condition of the
% ODE, and h is the stepsize.  You may also want to pass the function into
% the ODE the way |ode45| does (check lab 2).
%
% Note: you will need to use a loop to do this exercise.  
% You will also need to recall the Heun/Improved Euler algorithm learned in lectures.  
%

%% Exercise 2
%
% Objective: Compare Heun with |ode45|.
%
% Specifics:  For the following initial-value problems (from lab 2, 
% exercises 1, 4-6), approximate the solutions with your function from
% exercise 1 (Improved Euler Method).
% Plot the graphs of your Improved Euler Approximation with the |ode45| 
% approximation.

% (a) |y' = y tan t + sin t, y(0) = -1/2| from |t = 0| to |t = pi|

% (b) |y' = 1 / y^2 , y(1) = 1| from |t=1| to |t=10|

% (c) |y' =  1 - t y / 2, y(0) = -1| from |t=0| to |t=10|

% (d) |y' = y^3 - t^2, y(0) = 1| from |t=0| to |t=1|

% Comment on any major differences, or the lack thereof. You do not need
% to reproduce all the code here. Simply make note of any differences for
% each of the four IVPs.

% PART A
fa = @(t, y) y*tan(t) + sin(t);
t0a = 0;
y0a = -0.5;
t1a = pi;

[ta, ya] = heun(fa, t0a, t1a, y0a, 0.001);
soln = ode45(fa, [t0a, t1a], y0a);

plot(ta, ya, '-'); 
plot(soln.x, soln.y, '-', 'MarkerSize',10, 'LineWidth', 2);
xlabel('t');
ylabel('y');
legend('Improved Euler Method', 'ODE45', 'Location', 'Best');

figure;
hold on;

% PART B
fb = @(t, y) 1 ./ y^2;
t0b = 1;
y0b = 1;
t1b = 10;

[tb, yb] = heun(fb, t0b, t1b, y0b, 0.001);
soln1 = ode45(fb, [t0b, t1b], y0b);

plot(tb, yb, '-'); 
plot(soln1.x, soln1.y, '-', 'MarkerSize', 10, 'LineWidth', 2);

xlabel('t');
ylabel('y');
legend('Improved Euler Method', 'ODE45', 'Location', 'Best');

figure;
hold on;

% PART C
fc = @(t, y) 1 - (t*y)/2;
t0c = 0;
y0c = -1;
t1c = 10;
hc = 0.1;

[tc, yc] = heun(fc, t0c, t1c, y0c, 0.01);
soln2 = ode45(fc, [t0c, t1c], y0c);

plot(tc, yc, '-'); 
plot(soln2.x, soln2.y, '-', 'MarkerSize',10, 'LineWidth', 2);

xlabel('t');
ylabel('y');
legend('Improved Euler Method', 'ODE45', 'Location', 'Best');

figure;
hold on;

% PART D
fd = @(t, y) y^3 - t^2;
t0d = 0;
y0d = 1;
t1d = 1;
hd = 0.1;

[td, yd] = heun(fd, t0d, t1d, y0d, 0.01);
soln3 = ode45(fd, [t0d, t1d], y0d);

plot(td, yd, '-'); 
plot(soln3.x, soln3.y, '-', 'MarkerSize', 10, 'LineWidth', 2);

xlabel('t');
ylabel('y');

legend('Improved Euler Method', 'ode45', 'Location', 'Best');

% COMMENTS
% When using a step size of h=0.001, there were not many differences
% between ode45 and the part A and B. In part C, there were some 
% differences between 1 < t < 2. For part D, there were noticeable 
% deviations around t = 1. 

%% Exercise 3
%
% Objective: Use Euler's method and verify an estimate for the global error.
%
% Details: 
%
% (a) Use Euler's method (you can use
% euler.m from iode) to solve the IVP
%
% |y' = 2 t sqrt( 1 - y^2 )  ,  y(0) = 0|
%
% from |t=0| to |t=0.5|.
%
% (b) Calculate the solution of the IVP and evaluate it at |t=0.5|.
%
%     Solving as a first-order separable ODE, the solution to the IVP is 
%     y = sin(t^2)
%
% (c) Read the attached derivation of an estimate of the global error for 
%     Euler's method. Type out the resulting bound for En here in
%     a comment. Define each variable.
%
% (d) Compute the error estimate for |t=0.5| and compare with the actual
% error.
%
% (e) Change the time step and compare the new error estimate with the
% actual error. Comment on how it confirms the order of Euler's method.

% SOLUTIONS

% (b) Evaluating as a seperable equation, the solution is
% y(t) = sin(t^2 + c)
% y(0) = 0 -> y(t) = sin(t^2)

% (c) fMAX = 2*0.5*sqrt(1-sin(0.5)^2) = 0.8776
% ftMAX = 2*sqrt(1-sin(0.5)^2) = 1.7552
% fyMAX = 2*0.5*sin(0.5)/sqrt(1-sin(0.5)^2) = 0.5463
% Therefore M is approximately 1.7552

% (d) for dt=0.01, dt*n=0.5; so: En <= (1+M)*dt/2 * (exp(M*dt*n) - 1) ~ 0.0194
% The actual error was 0.0047, which is smaller than
% the upper bound

% (e) For dt=0.01 the error was  approx. 0.0047, and for dt=0.1 the 
% error was  approx. 0.048. This shows that increasing the step size by a
% factor of 10, will also increase the error by a factor of 10, which is a
% linear relationship, as stated in the document

close all; clear; clc;

hold on;
title('Euler Method');

clc;
f = @(t, y) 2 * t * sqrt(1 - y^2);
y0 = 0;

step = 0.01;
td = 0 : step : 0.5;
y = euler(f, y0, td);
fprintf('Error with step=%g is %g\n', step, abs(sin(0.5^2) - y(end)));
plot(td, y, 'r-', 'LineWidth', 1);

step = 0.1;
td = 0 : step : 0.5;
y = euler(f, y0, td);
fprintf('Error with step=%g is %g\n', step, abs(sin(0.5^2) - y(end)));
plot(td, y, 'b-', 'LineWidth', 1);

plot(td, sin(td.^2), 'g--', 'LineWidth', 1);

xlabel('t');
ylabel('y');
legend('Step=0.01', 'Step=0.1', 'Actual', 'Location', 'NorthWest');

%% Adaptive Step Size
%
% As mentioned in lab 2, the step size in |ode45| is adapted to a
% specific error tolerance.
%
% The idea of adaptive step size is to change the step size |h| to a
% smaller number whenever the derivative of the solution changes quickly.
% This is done by evaluating f(t,y) and checking how it changes from one
% iteration to the next.

%% Exercise 4
%
% Objective: Create an Adaptive Euler method, with an adaptive step size |h|.
%
% Details: Create an m-file which accepts the variables |(t0,tN,y0,h)|, as 
% in exercise 1, where |h| is an initial step size. You may also want to 
% pass the function into the ODE the way |ode45| does.
%
% Create an implementation of Euler's method by modifying your solution to 
% exercise 1. Change it to include the following:
%
% (a) On each timestep, make two estimates of the value of the solution at
% the end of the timestep: |Y| from one Euler step of size |h| and |Z| 
% from two successive Euler steps of size |h/2|. The difference in these
% two values is an estimate for the error.
%
% (b) Let |tol=1e-8| and |D=Z-Y|. If |abs(D)<tol|, declare the step to be
% successful and set the new solution value to be |Z+D|. This value has
% local error |O(h^3)|. If |abs(D)>=tol|, reject this step and repeat it 
% with a new step size, from (c).
%
% (c) Update the step size as |h = 0.9*h*min(max(tol/abs(D),0.3),2)|.
%
% Comment on what the formula for updating the step size is attempting to
% achieve.

% COMMENTS
% The Euler method works more accurately the smaller the step size is. So
% by decreasing the step size it produces more precise results.

%% Exercise 5
%
% Objective: Compare Euler to your Adaptive Euler method.
%
% Details: Consider the IVP from exercise 3.
%
% (a) Use Euler method to approximate the solution from |t=0| to |t=0.75|
% with |h=0.025|.
%
% (b) Use your Adaptive Euler method to approximate the solution from |t=0| 
% to |t=0.75| with initial |h=0.025|.
%
% (c) Plot both approximations together with the exact solution.

f = @(t, y) 2 * t * sqrt(1 - y^2);

% Euler Method
te = 0 : 0.025 : 0.75;
ye = euler(f, 0, te);

% Adaptive Euler Method
[ta, ya] = AdaptiveEuler(f, 0, 0.75, 0, 0.025);

% Plotting the solutions
hold on;
title('Adaptive Euler Method');
xlabel('t');
ylabel('y');

plot(te, ye, 'g-', 'LineWidth', 3);
plot(ta, ya, 'r-', 'LineWidth', 3);
plot(te, sin(te.^2), 'k--', 'LineWidth', 3);

legend('Euler Method', 'Adaptive Euler Method', 'Exact Solution', 'Location',...
    'NorthWest');

%% Exercise 6
%
% Objective: Problems with Numerical Methods.
%
% Details: Consider the IVP from exercise 3 (and 5).
% 
% (a) From the two approximations calculated in exercise 5, which one is
% closer to the actual solution (done in 3.b)? Explain why.
% 
% (b) Plot the exact solution (from exercise 3.b), the Euler's 
% approximation (from exercise 3.a) and the adaptive Euler's approximation 
% (from exercise 5) from |t=0| to |t=1.5|.
% 
% PLOT exact, euler's approx, adaptive euler's approximation
%
% (c) Notice how the exact solution and the approximations become very
% different. Why is that? Write your answer as a comment.

f = @(t, y) 2 * t * sqrt(1 - y^2);

% Euler Method
te = 0 : 0.025 : 1.5;
ye = euler(f, 0, te);

% Advanced Euler Method
[ta, ya] = AdaptiveEuler(f, 0, 1.5, 0, 0.025);

% Plotting the solutions
hold on;
title('Adaptive Euler Method');
xlabel('t');
ylabel('y');

plot(te, ye, 'g-', 'LineWidth', 3);
plot(ta, ya, 'r-', 'LineWidth', 3);
plot(te, sin(te.^2), 'k--', 'LineWidth', 3);

legend('Euler Method', 'Adaptive Euler Method', 'Exact Solution', 'Location',...
    'NorthWest');

% COMMENTS:
% 
% The Adaptive Euler Method is closer to the actual solution as the
% stepsize is adjusted. As the stepsize decreases, the accuracy increases.%
%
% The Euler method does not work well and deviates from the actual solution
% as the error becomes too large. Though not as significant, a similar
% problem occurs with the Adapative Euler Method, as the error indicator is
% only an approximation, so it becomes smaller than the actual error and
% the stepsize doesn't decrease enough to stay within the exact solution.