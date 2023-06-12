function y = AdamBash3(f, y0, h)
%   Approximate solution to y' = f(t,y) using second-order Adam Bashforth method.
%   Inputs:
%       f       function f(t,y)
%       y0      Initial value y(0)
%       h       Time step
%
%   Output:
%       y       Approximate solution at t = 1

    % Number of time steps
    numberofintervals = (1-y0)/h;
    % Time vector
    t = 0:h:1;

    % Initialize y to a vector of the initial value y0 at each interval
    y = zeros(size(t));
    y(1) = y0;
    y(2) = y0;
    y(3) = y0;

    if h == 1
        y = y(end);
        return
    end
    
    % Loop through remaining time steps and use Adam Bashforth method to
    % compute approximate solution
    for i = 4:numberofintervals+1
        y(i) = y(i-1) + h/12 * (23*f(t(i-1), y(i-1)) - 16*f(t(i-2), y(i-2)) + 5*f(t(i-3), y(i-3)));
    end
    
    % Return the solution at the final time
    y = y(end);
end