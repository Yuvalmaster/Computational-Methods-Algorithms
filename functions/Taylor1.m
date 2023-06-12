function y = Taylor1(f, y0, h)
%   Approximate solution to y' = f(t,y) using first-order Taylor method.
%   Inputs:
%       f       function f(t,y)
%       y0      Initial value y(0)
%       h       Time step
%
%   Output:
%       y       Approximate solution at t = 1

    % Number of time steps
    numberofintervals = (1-y0)/h;
    
    % Initialize solution
    y = y0;
    
    % Loop over time steps
    for i = 1:numberofintervals
        % Advance solution using first-order Taylor method
        y = y + h * f(i*h, y);
    end
end