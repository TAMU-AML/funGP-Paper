%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DESCRIPTION:
%   main function for Cox & Lee (2008) simulation for Type 1 error using gp sample  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_simulation = 1000; %number of simualtions
n_type1 = 0 ; %variable to store the number of simulations with type 1 error
n_type1_funGP = 0;
confLevel = 0.95; %significance level
logfile = "CL_gpSample_Type1.log";
logID = fopen(logfile,'w');

fprintf('Cox and Lee experiment for Type 1 error using GP sample function\n')
fprintf('Storing computation log in file: %s\n',logfile);

%model specification
sigma_f = 5; %true value of sigma_f
sigma_n = 0.5; %true value of sigma_n
theta = 0.2; %true value of theta
n_function_eval = 50; %number of points to evaluate the gp sample function
n_replications = 10;
rng(1); %setting the seed 

n_CL_SigPoints = zeros(n_simulation,1); %vector to store number of stat sig points for CL test
n_funGP_SigPoints = zeros(n_simulation,1); %vector to store number of stat sig points for funGP test
corrected_p_values_mat = zeros(n_function_eval,n_simulation); % matrix to store corrected p-values for CL test

for j = 1:n_simulation
    X = linspace(0,1,n_function_eval)'; %points to evaluate the gp functions 
    mu = zeros(n_function_eval,1); %mean function value
    covMat = (sigma_f^2)*computeCorrelMat(X,X,theta); %covariance matrix for the function f(.)
    F = mvnrnd(mu,covMat,1)' ; %vector of function evaluations
    Y1_noisy = zeros(n_function_eval,n_replications);
    Y2_noisy = zeros(n_function_eval,n_replications);
    Y1_smooth = zeros(n_function_eval,n_replications);
    Y2_smooth = zeros(n_function_eval,n_replications);
    for i = 1:n_replications
        Y1_noisy(:,i) = F + normrnd(0,sigma_n,n_function_eval,1);
        gpMdl1 = fitrgp(X,Y1_noisy(:,i));
        Y1_smooth(:,i) = predict(gpMdl1,X);
        Y2_noisy(:,i) = F + normrnd(0,sigma_n,n_function_eval,1);
        gpMdl2 = fitrgp(X,Y2_noisy(:,i));
        Y2_smooth(:,i) = predict(gpMdl2,X);
    end
    
    nRandomize = 1000;
    out_CL = CoxLeeTest(Y1_smooth,Y2_smooth,X,confLevel, nRandomize);
    if out_CL.differ == true
        n_type1 = n_type1 + 1;        
    end
    n_CL_SigPoints(j) = length(out_CL.statSigPointIndices);
    fprintf(logID,'Number of statistically significant test points for CL: %d\n',n_CL_SigPoints(j));
    corrected_p_values_mat(:,j) = out_CL.corrected_p_values; 
    fprintf(logID,'Number of replications with Type 1 error for CL: %d out of %d simulations\n',n_type1,j);
    
    %%funGP algorithm

    %vectorizing the response matrix
    Y1_noisy = Y1_noisy(:);
    Y2_noisy = Y2_noisy(:);
    
    %replicating X 
    X_rep = repmat(X,n_replications,1);
    
    out_funGP = funGP(X_rep,Y1_noisy,X_rep,Y2_noisy,X,confLevel);
    if out_funGP.differ == true
        n_type1_funGP = n_type1_funGP + 1;        
    end
    n_funGP_SigPoints(j) = out_funGP.nPoints;
    fprintf(logID,'Number of statistically significant test points for funGP: %d\n',n_funGP_SigPoints(j));
    fprintf(logID,'Number of replications with Type 1 error for funGP: %d out of %d simulations\n',n_type1_funGP,j);
end

fprintf('Estimated Type 1 error for CL: %.3f\n',n_type1/n_simulation);
fprintf('Estimated Type 1 error for funGP: %.3f\n',n_type1_funGP/n_simulation);
