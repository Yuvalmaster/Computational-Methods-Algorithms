function y = RungeKutta4(f, y0, h)
%   Approximate solution to y' = f(t,y) using second-order Runge-Kutta method.
%   Inputs:
%       f       function f(t,y)
%       y0      Initial value y(0)
%       h       Time step
%
%   Output:
%       y       Approximate solution at t = 1

    % Number of time steps
    numberofintervals = (1-y0)/h;

    % Initialize y to a vector of the initial value y0 at each interval
    y = y0;
    
    % Iterate over the intervals
    for i = 1:numberofintervals
        % Compute the next value of y using the Runge-Kutta algorithm
        k1 = h*f(i*h, y);
        k2 = h*f(i*h + h/2, y + k1/2);
        k3 = h*f(i*h + h/2, y + k2/2);
        k4 = h*f(i*h + h  , y + k3  );
        y = y + (1/6)*(k1 + 2*k2 + 2*k3 +  k4);
    end


end