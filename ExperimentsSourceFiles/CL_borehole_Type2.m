%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DESCRIPTION:
%   main function for Cox and Lee (2008) simulation for Type 2 error for borehole function 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../funGP-Matlab'); %adding path for the subroutines
n_simulation = 1000; %number of simualtions
n_type2 = 0 ; %variable to store the number of simulations with type 2 error
n_type2_funGP = 0; %variable to store the number of simulations with type 2 error for funGP
confLevel = 0.95; %significance level
n_eval_points = 10; %number of points at which the functions are evaluated
n_replications = 10; %number of replications at each eval point.
n_CL_SigPoints = zeros(n_simulation,1); %vector to store number of stat sig points for CL test
n_funGP_SigPoints = zeros(n_simulation,1); %vector to store number of stat sig points for funGP test
corrected_p_values_mat = zeros(n_eval_points^2,n_simulation); % matrix to store corrected p-values for CL test

logfile = "CL_borehole_Type2.log";
logID = fopen(logfile,'w');

fprintf('Cox and Lee experiment for Type 2 error using borehole function\n')
fprintf('Storing computation log in file: %s\n',logfile);

%Fixed input values 
Tu = 78000;
Hu = 1050;
Tl = 84;
Hl = 760;
L  = 1400;
Kw = 11000;
Lmod = L + 50;

rng(2); %setting the seed 
j = 1;
while j <= n_simulation
    %Variable input values for two datasets
    rw0 = linspace(0.05,0.15,n_eval_points);
    r0 = linspace(100,50000,n_eval_points);

    
    %constructing the datasets
    X = combvec(rw0, r0)';
    length_X = size(X,1);
    
    Y1_noisy = zeros(length_X,n_replications);  
    Y2_noisy = zeros(length_X,n_replications);  
    Y1_smooth = zeros(length_X,n_replications);
    Y2_smooth = zeros(length_X,n_replications);
    
    for nrows = 1:length_X
        Y1_noisy(nrows,:) = borehole([X(nrows,1),X(nrows,2),Tu,Hu,Tl,Hl,L,Kw])+ normrnd(0,10,1,n_replications);
        Y2_noisy(nrows,:) = borehole([X(nrows,1),X(nrows,2),Tu,Hu,Tl,Hl,Lmod,Kw])+ normrnd(0,10,1,n_replications);
    end

    for replicates = 1:n_replications
        gpMdl1 = fitrgp(X,Y1_noisy(:,replicates));
        Y1_smooth(:,replicates) = predict(gpMdl1,X);
        gpMdl2 = fitrgp(X,Y2_noisy(:,replicates));
        Y2_smooth(:,replicates) = predict(gpMdl2,X);
    end

    nRandomize = 1000;
    out = CoxLeeTest(Y1_smooth,Y2_smooth,X,confLevel,nRandomize);
    if out.differ == false
        n_type2 = n_type2 + 1;        
    end
    n_CL_SigPoints(j) = length(out.statSigPointIndices);
    fprintf(logID,'Number of statistically significant test points for CL: %d\n',n_CL_SigPoints(j));
    corrected_p_values_mat(:,j) = out.corrected_p_values; 
    fprintf(logID,'Number of replications with Type 2 error for CL: %d out of %d simulations\n',n_type2,j);
    
    %%funGP algorithm

    %vectorizing the response matrix
    Y1_noisy = Y1_noisy(:);
    Y2_noisy = Y2_noisy(:);

    %replicating X 
    X_rep = repmat(X,n_replications,1);
    
    try
        out_funGP = funGP(X_rep, Y1_noisy, X_rep, Y2_noisy, X, confLevel);
    catch
        fprintf(logID,"Error in parameter estimation. Skipping this iteration\n");
        continue;
    end
    if out_funGP.differ == false
        n_type2_funGP = n_type2_funGP + 1;        
    end
    n_funGP_SigPoints(j) = out_funGP.nPoints;
    fprintf(logID,'Number of statistically significant test points for funGP: %d\n',n_funGP_SigPoints(j));
    fprintf(logID,'Number of replications with Type 2 error for funGP: %d out of %d simulations\n',n_type2_funGP,j);
    j = j+1;
end

fprintf('Estimated Type 2 error for CL: %.3f\n',n_type2/n_simulation);
fprintf('Estimated Type 2 error for funGP: %.3f\n',n_type2_funGP/n_simulation);
