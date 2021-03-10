%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DESCRIPTION:
%   main function for Cox & Lee (2008) simulation for Type 2 error using gp sample  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_simulation = 1000; %number of simualtions
n_type2 = 0 ; %variable to store the number of simulations with type 1 error
n_type2_funGP = 0;
confLevel = 0.95; %significance level
logfile = "CL_gpSample_Type2.log";
logID = fopen(logfile,'w');

fprintf('Cox and Lee experiment for Type 2 error using GP sample function\n')
fprintf('Storing computation log in file: %s\n',logfile);

%model specification
sigma_f = 5; %true value of sigma_f
sigma_n = 0.5; %true value of sigma_n
theta = 0.2; %true value of theta
n_function_eval = 50; %number of points to evaluate the gp sample function
n_replications = 10;
l2dist_vector = zeros(n_simulation,1); %array to store percentage l2 distance for each simulation
rng(2); %setting the seed 
n_CL_SigPoints = zeros(n_simulation,1); %vector to store number of stat sig points for CL test
n_funGP_SigPoints = zeros(n_simulation,1); %vector to store number of stat sig points for funGP test
corrected_p_values_mat = zeros(n_function_eval,n_simulation); % matrix to store corrected p-values for CL test

for j = 1:n_simulation
    X = linspace(0,1,n_function_eval)'; %points to evaluate the gp functions 
    mu = zeros(n_function_eval,1); %mean function value
    covMat = (sigma_f^2)*computeCorrelMat(X,X,theta); %covariance matrix for the function f(.)
    F = mvnrnd(mu,covMat,1)' ; %vector of function evaluations
    FSq = F.^2; %square of the function
    l2norm = sqrt(trapz(X,FSq)); %l2 norm of the function calculated using trapezoidal rule
    perturbationRegion = [0.2,0.8]; %perturbation region
    %function handle for the perturbation function
    fun = @(x) 0.33*sin(pi*(x-perturbationRegion(1))/(perturbationRegion(2)-perturbationRegion(1))); 
    l2dist = ComputeL2(fun,0.2,0.8); %computing l2 distance
    percentDiff = l2dist*100/l2norm; %computing percent l2 distance
    fprintf(logID,'Percent L2 distance for the current simulation: %2.2f\n',percentDiff);
    l2dist_vector(j) = percentDiff; %storing percent l2 distance 
    perturb_index = find(X>=perturbationRegion(1) & X<=perturbationRegion(2)); 
    G = F;
    G(perturb_index)= G(perturb_index) + fun(X(perturb_index)); %perturbing response for dataset2
    Y1_noisy = zeros(n_function_eval,n_replications);
    Y2_noisy = zeros(n_function_eval,n_replications);
    Y1_smooth = zeros(n_function_eval,n_replications);
    Y2_smooth = zeros(n_function_eval,n_replications);

    for i = 1:n_replications
        Y1_noisy(:,i) = F + normrnd(0,sigma_n,n_function_eval,1);
        gpMdl1 = fitrgp(X,Y1_noisy(:,i));
        Y1_smooth(:,i) = predict(gpMdl1,X);
        Y2_noisy(:,i) = G + normrnd(0,sigma_n,n_function_eval,1);
        gpMdl2 = fitrgp(X,Y2_noisy(:,i));
        Y2_smooth(:,i) = predict(gpMdl2,X);
    end
    
    nRandomize = 1000;
    out_CL = CoxLeeTest(Y1_smooth,Y2_smooth,X,confLevel, nRandomize);
    if out_CL.differ == false
        n_type2 = n_type2 + 1;        
    end
    n_CL_SigPoints(j) = length(out_CL.statSigPointIndices);
    fprintf(logID,'Number of statistically significant test points for CL: %d\n',n_CL_SigPoints(j));
    corrected_p_values_mat(:,j) = out_CL.corrected_p_values; 
    fprintf(logID,'Number of replications with Type 2 error for CL: %d out of %d simulations\n',n_type2,j);
  
    %%funGP algorithm

    %vectorizing the response matrix
    Y1_noisy = Y1_noisy(:);
    Y2_noisy = Y2_noisy(:);
    
    %replicating X 
    X_rep = repmat(X,n_replications,1);
    
    out_funGP = funGP(X_rep,Y1_noisy,X_rep,Y2_noisy,X,confLevel);
    if out_funGP.differ == false
        n_type2_funGP = n_type2_funGP + 1;        
    end
    n_funGP_SigPoints(j) = out_funGP.nPoints;
    fprintf(logID,'Number of statistically significant test points for funGP: %d\n',n_funGP_SigPoints(j));
    fprintf(logID,'Number of replications with Type 2 error for funGP: %d out of %d simulations\n',n_type2_funGP,j);
end

fprintf('Estimated Type 2 error for CL: %.3f\n',n_type2/n_simulation);
fprintf('Estimated Type 2 error for funGP: %.3f\n',n_type2_funGP/n_simulation);