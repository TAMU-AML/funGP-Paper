%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DESCRIPTION:
%   main function for funGP simulation for Type 2 error using gp sample  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ntrain = 500; %number of training samples
ntest = 500; %number of test samples
n_simulation = 1000; %number of simualtions
n_type2 = 0 ; %variable to store the number of simulations with type 2 error
confLevel = 0.95; %significance level

fprintf("Simulation for function: %s\n","GP sample");
fprintf("Type of test: %s\n","Type 2");
fprintf("Sample size for each dataset: %d\n",ntrain);
fprintf("Size of test grid: %d\n",ntest);
fprintf("Confidence level: %.7f\n",confLevel);
fprintf("Number of simulations: %d\n",n_simulation);

logfile = "gpSample_Type2.log";
logID = fopen(logfile,'w');
fprintf('Storing computation log in file: %s\n',logfile);

%model specification
sigma_f = 5; %true value of sigma_f
sigma_n = 0.5; %true value of sigma_n
theta = 0.2; %true value of theta
n_function_eval = 5000; %number of points to evaluate the gp sample function
l2dist_vector = zeros(1000,1); %array to store percentage l2 distance for each simulation
rng(2); %setting the seed 

for j = 1:n_simulation
    x = linspace(0,1,n_function_eval)'; %points to evaluate the gp functions 
    mu = zeros(n_function_eval,1); %mean function value
    covMat = (sigma_f^2)*computeCorrelMat(x,x,theta); %covariance matrix for the function f(.)
    F = mvnrnd(mu,covMat,1)' ; %vector of function evaluations
    FSq = F.^2; %square of the function
    l2norm = sqrt(trapz(x,FSq)); %l2 norm of the function calculated using trapezoidal rule
    perturbationRegion = [0.2,0.8]; %perturbation region
    %function handle for the perturbation function
    fun = @(x) 0.33*sin(pi*(x-perturbationRegion(1))/(perturbationRegion(2)-perturbationRegion(1))); 
    l2dist = ComputeL2(fun,0.2,0.8); %computing l2 distance
    percentDiff = l2dist*100/l2norm; %computing percent l2 distance
    fprintf(logID,'Percent L2 distance for the current simulation: %.2f\n',percentDiff);
    l2dist_vector(j) = percentDiff; %storing percent l2 distance 
    index1 = randsample(1:n_function_eval,ntrain)'; %sampling index of dataset1
    index2 = randsample(1:n_function_eval,ntrain)'; %sampling index for dataset2
    x1 = x(index1); %input points for dataset1
    x2 = x(index2); %input points for dataset2
    y1 = F(index1) + normrnd(0,sigma_n,ntrain,1); %response for dataset1
    y2 = F(index2)+normrnd(0,sigma_n,ntrain,1); %response for dataset2
    %finding the index of input points for dataset2 to perturb
    perturb_index = find(x2>=perturbationRegion(1) & x2<=perturbationRegion(2)); 
    y2(perturb_index)= y2(perturb_index) + fun(x2(perturb_index)); %perturbing response for dataset2
    xtest = linspace(0,1,ntest)'; %creating test dataset 
    
    out = funGP(x1, y1, x2, y2, xtest, confLevel);
    if out.differ == false
        n_type2 = n_type2 + 1;
    end
    fprintf(logID,'Number of replications with Type 2 error: %d out of %d simulations\n',n_type2,j);  
end
fprintf('Estimated Type 2 error: %.3f\n',n_type2/n_simulation);
fprintf('Average percent L2 distance between the functions: %.2f\n',mean(l2dist_vector));