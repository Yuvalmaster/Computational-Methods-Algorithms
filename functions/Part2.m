%% Set constants & Variables
f=@(x) cos(x);                  % Set function
a = 1; b = 2; real_sol = pi/2;  % Limits & real solution
epsilon       = 0.001;          % Threshold
N_interations = 1000;           % Max iteration

Methods = {'Bisection', 'Secant', 'Fixed-Point', 'Newton-Raphson'};
funcs   = {@Bisection, @Secant, @FixedPoint, @NewtonRaphson};

X  = zeros(size(Methods,2),1);        % Vector of x results for f(x)==0 in each method
ii = zeros(size(Methods,2),1);        % Vecor of numbers of iteration for each method
diff_list = zeros(size(Methods,2),1); % Difference vector between real result and approximate result in each method

%% Run functions
%   This section loads each method and gives the approximate x solution,
%   number of iterations, and output lists for analysis table.
%   Inputs:
%       f              function f
%       epsilon        Method's threshold for convergance
%       N_iterations   maximum number of iterations
%       a              bottom limit
%       b              top limit
%       funcs          set of algorithms' functions to run
%       Methods        list of names for each function
%
%   Output:
%       X             list of solutions for each method
%       ii            list of number of iterations for each method
%       diff_list     list of differences between method's solution and real solution

for method_i=1:length(funcs)
    disp(append('<strong>', Methods{method_i} ,' Method</strong>'));
    tic
    [X(method_i), ii(method_i)] = funcs{method_i}(f, a, b, epsilon, N_interations);
    toc
    % Create differences list from real solution
    diff_list(method_i) = abs(real_sol - X(method_i));
    
    % Display stats 
    disp(append('F(X) ≅ 0 for X = '     ,num2str(round(X(method_i),4))));
    disp(append('Number of iterations: ', num2str(ii(method_i)))       );
    disp(append('Difference from π/2: ' , num2str(diff_list(method_i)))); 
    disp(' ');

end

%% Export data to table
T = table(X,ii,diff_list,...
         'VariableNames',{'X','Number of Iterations','Difference from π/2'},...
         'RowNames',Methods);

% Adding real solution
T = [T; table(real_sol, 0, 0,...
    'VariableNames',{'X','Number of Iterations','Difference from π/2'},...
    'RowNames',{'Real solution'})];

writetable(T,'MethodsAnalysis2.csv','WriteRowNames',true)
