%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DESCRIPTION:
%   main function for Munk & Dette (1998) simulation for Type 1 error using gp sample  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ntrain = 500; %number of training samples
n_simulation = 1000; %number of simualtions
n_type1 = 0 ; %variable to store the number of simulations with type 1 error
confLevel = 0.95; %significance level

fprintf('Munk and Dette experiment for Type 1 error using gpSample function\n')
fprintf('Size of train data: %d\n', ntrain);
fprintf('Nominal level of the test: %.7f\n', confLevel);
fprintf("Number of simulations: %d\n",n_simulation);

logfile = "MD_gpSample_Type1.log";
logID = fopen(logfile,'w');
fprintf('Storing computation log in file: %s\n',logfile);

%model specification
sigma_f = 5; %true value of sigma_f
sigma_n = 0.5; %true value of sigma_n
theta = 0.2; %true value of theta
n_function_eval = 5000; %number of points to evaluate the gp sample function

rng(1); %setting the seed 

for j = 1:n_simulation
    x = linspace(0,1,n_function_eval)'; %points to evaluate the gp functions 
    mu = zeros(n_function_eval,1); %mean function value
    covMat = (sigma_f^2)*computeCorrelMat(x,x,theta); %covariance matrix for the function f(.)
    F = mvnrnd(mu,covMat,1)' ; %vector of function evaluations
    index1 = randsample(1:n_function_eval,ntrain)'; %sampling index of dataset1
    index2 = randsample(1:n_function_eval,ntrain)'; %sampling index for dataset2
    x1 = x(index1); %input points for dataset1
    x2 = x(index2); %input points for dataset2
    y1 = F(index1) + normrnd(0,sigma_n,ntrain,1); %response for dataset1
    y2 = F(index2)+normrnd(0,sigma_n,ntrain,1); %response for dataset2
    
    out = MunkDetteTest(x1,y1,x2,y2,confLevel, [0,1]);
    
    if out.differ == true
        n_type1 = n_type1 + 1;        
    end
    fprintf(logID,'Number of replications with Type 1 error: %d out of %d simulations\n',n_type1,j);
end
fprintf('Estimated Type 1 error: %.3f\n',n_type1/n_simulation);
