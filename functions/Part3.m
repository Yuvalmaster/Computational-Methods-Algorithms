%% Set constants & Variables
S=...
[-0.40,+0.20,+0.00,+0.00,+0.20,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00;...
 +0.20,-0.53,+0.20,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00;...
 +0.00,+0.20,-0.53,+0.20,+0.00,+0.13,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00;...
 +0.00,+0.00,+0.20,-0.40,+0.00,+0.00,+0.20,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00;...
 +0.20,+0.00,+0.00,+0.00,-0.53,+0.00,+0.00,+0.20,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00;...
 +0.00,+0.00,+0.13,+0.00,+0.00,-0.46,+0.13,+0.00,+0.00,+0.10,+0.00,+0.00,+0.00,+0.00,+0.00;...
 +0.00,+0.00,+0.00,+0.20,+0.00,+0.13,-0.53,+0.00,+0.00,+0.00,+0.20,+0.00,+0.00,+0.00,+0.00;...
 +0.00,+0.00,+0.00,+0.00,+0.20,+0.00,+0.00,-0.53,+0.13,+0.00,+0.00,+0.20,+0.00,+0.00,+0.00;...
 +0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.13,-0.46,+0.10,+0.00,+0.00,+0.13,+0.00,+0.00;...
 +0.00,+0.00,+0.00,+0.00,+0.00,+0.10,+0.00,+0.00,+0.10,-0.46,+0.13,+0.00,+0.00,+0.13,+0.00;...
 +0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.20,+0.00,+0.00,+0.13,-0.53,+0.00,+0.00,+0.00,+0.20;...
 +0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.20,+0.00,+0.00,+0.00,-0.40,+0.20,+0.00,+0.00;...
 +0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.13,+0.00,+0.00,+0.20,-0.53,+0.20,+0.00;...
 +0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.13,+0.00,+0.00,+0.20,-0.53,+0.20;...
 +0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.00,+0.20,+0.00,+0.00,+0.20,-0.40];

Is = [0;0.5;0.5;0;0;0;0;0;0;0;0;0;-0.5;-0.5;0];
x0 = zeros(15,1);                               % Initial guess  (φ'0)

f = @(x) S*x - Is;                              % Set function
real_sol      = linsolve(S,Is);                 % Find the real solution

tol           = 0.001;                          % Threshold
N_iterations  = 1000;                           % Max iteration

Methods = {'Jacobi', 'Gauss-Seidel', 'SOR w = 1' ,...
           'SOR w = 1.5', 'Binconjugate gradient',...
           'QR Factorisation','Real Solution'       };

funcs   = {@Jacobi, @GaussSeidel, @SOR, @SOR, @HestenesSteifel, @QRFact};

X  = zeros(size(funcs,2),size(x0,1));  % Vector of x results of each method
ii = zeros(size(funcs,2),1);           % Vecor of numbers of iteration for each method
diff_list = zeros(size(funcs,2),1);    % Difference vector between real result and approximate result in each method

%% Run functions
%   This section loads each method and gives the approximate x solution,
%   number of iterations, and output lists for analysis table.
%   Inputs:
%       S              function S'
%       tol            Method's threshold for convergance
%       N_iterations   maximum number of iterations
%       x0             Initial guess
%       Is             vector Is' of the function S'φ'=Is'
%       w              coefficient for SOR method
%       funcs          set of algorithms' functions to run
%       Methods        list of names for each function
%
%   Output:
%       X             list of solutions for each method
%       ii            list of number of iterations for each method
%       diff_list     list of differences between method's solution and real solution

for method_i=1:numel(funcs)
    disp(append('<strong>', Methods{method_i} ,' Method</strong>'));
    % Set different w for SOR algorithm
    if contains(Methods{method_i},'SOR')
        % Choose w from Method's name
        w = str2double(regexp(Methods{method_i},'\d+\.?\d*','Match'));
        tic
        [X(method_i,:), ii(method_i)] = funcs{method_i}(x0, S, Is, tol, N_iterations, w);
        toc
    else
        tic
        [X(method_i,:), ii(method_i)] = funcs{method_i}(x0, S, Is, tol, N_iterations);
        toc
    end
    
    % Create differences list from real solution
    difference            = setdiff(X(method_i,:),real_sol);
    diff_list(method_i,1) = mean(difference);
    
    % Display stats
    disp(append('A*X-b ≅ 0 for X = '            , num2str(round(X(method_i,:),4))));
    disp(append('Number of iterations: '        , num2str(ii(method_i)))           );
    disp(append('Mean difference from real X: ' , num2str(diff_list(method_i)))    ); 
    disp(' ');

end
ii(end+1,1)        = NaN; % Add for real Solution
diff_list(end+1,1) = NaN; % Add for real Solution

%% Export data to table
row_names        = cell(size(real_sol));
for i=1:numel(real_sol)
    row_names{i} = append('X',num2str(i));
end

row_names{end+1} = 'Number of Iterations for convergance';
row_names{end+1} = 'Difference from real Solution'       ;
table_array      = [X' real_sol; ii'; diff_list']        ;

T = array2table(table_array, 'RowNames', row_names, 'VariableNames', Methods);
writetable(T, 'MethodsAnalysis3.csv', 'WriteRowNames', true)
