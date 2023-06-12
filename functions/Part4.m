%% Set constants & Variables
% Set function
f = @(t,y) (exp(t)+y);

% Set start term
y0 = 0;

% Set real solution for y(t==1) = e
y_real = exp(1);

% Number of intervals vector
h = [1, 0.25, 0.125, 0.01];

% Number of time steps
numberofintervals = (1-y0)./h;

% Methods vector
Methods = {'Taylor n=1', 'Runge Kutta n=2', 'Runge Kutta n=4' ,...
           'Adams-Bashforth m=2'      , 'Adams-Bashforth m=3' ,...
           'Predictor-corrector m=2','Predictor-corrector m=3'};

funcs   = {@Taylor1, @RungeKutta2, @RungeKutta4, @AdamBash2, @AdamBash3, @PredCorr2, @PredCorr3};

Y         = zeros(size(Methods,2),numel(h));  % Vector of y results for y(t=1) in each method for each h
diff_list = zeros(size(Methods,2),numel(h));  % Difference vector between real result and approximate result in each method for each h
time      = zeros(size(Methods,2),numel(h));  % Run time vector for each method for each number of intervals
%% Run functions
%   This section loads each method and gives the approximate y(t=1) 
%   solution, number of iterations, and output lists for analysis table.
%   Inputs:
%       f              function f
%       h              list of intervals coefficient
%       y0             start term
%       funcs          set of algorithms' functions to run
%       Methods        list of names for each function
%
%   Output:
%       Y             list of solutions for each method for each number of intervals
%       time          list of runtime for each method for each number of intervals
%       diff_list     list of differences between method's solution and real solution for each number of intervals

for method_i=1:numel(funcs)
    disp(append('<strong>', Methods{method_i} ,' Method</strong>'));
    for ii=1:length(h)
        % Choose current interval
        h_current = h(ii);  
        
        % Run algorithm & calculate run time for current interval
        tic
        Y(method_i,ii) = funcs{method_i}(f, y0, h_current);
        time(method_i,ii) = toc;

        % Calculate difference between real solution
        diff_list(method_i,ii) = abs(y_real - Y(method_i,ii));

        % Display stats
        disp(append(num2str(numberofintervals(ii)), ' intervals: Y(t=1) = ' , num2str(round(Y(method_i,ii),4))));
        disp(append('Run time: ', num2str(time(method_i,ii)), ' [Sec]'));
        disp(append('Difference from real solution: ', num2str(diff_list(method_i,ii)))); 
        disp(' ');
    end
end


%% Export data to table
columns = strcat('h=',sprintfc('%.3f',h));

% Solutions table
T_1 = splitvars(table(Y,'RowNames',Methods));
T_1.Properties.VariableNames = columns;
writetable(T_1,'MethodsAnalysis4_solutions.csv','WriteRowNames',true)

% differences table
T_2 = splitvars(table(diff_list,'RowNames',Methods));
T_2.Properties.VariableNames = columns;
writetable(T_2,'MethodsAnalysis4_differences.csv','WriteRowNames',true)

% runtimes table
T_3 = splitvars(table(time,'RowNames',Methods));
T_3.Properties.VariableNames = columns;
writetable(T_3,'MethodsAnalysis4_runtime.csv','WriteRowNames',true)