%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   AUTHOR: ABHINAV PRAKASH
%   DESCRIPTION:
%   main function for funGP sample size experiment simulation for Type 1 error for piston function 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[] = SampleSizeExperiment_piston_Type1(ntrain)
   
    logfile = strcat('SampleSize_',num2str(ntrain),'_piston_Type1.log');
    logID = fopen(logfile,'w');
    n_simulation = 1000; %number of simualtions
    n_type1 = 0 ; %variable to store the number of simulations with type 1 error
    confLevel = 0.95; %significance level
    ntest = 50*50;

    fprintf('Sample size experiment for Type 1 error using piston function\n')
    fprintf('Size of train data: %d\n', ntrain);
    fprintf('Size of test data: %d\n', ntest);
    fprintf('Nominal level of the test: %.7f\n', confLevel);

    fprintf('Storing computation log in file: %s\n',logfile);

    %Fixed input values 
    M  = 45;
    S  = 0.010;
    k  = 2000;
    P0 = 100000;
    Ta = 292;

    rng(1); %setting the seed 
    j = 1;
    while j <= n_simulation
        %Variable input values for two datasets
        V01 = 0.002 + (0.010 - 0.002).*rand(ntrain,1);
        T01 = 340 + (360 - 340).*rand(ntrain,1);
        V02 = 0.002 + (0.010 - 0.002).*rand(ntrain,1);
        T02 = 340 + (360 - 340).*rand(ntrain,1);

        %constructing first dataset
        x1 = [V01,T01];
        y1 = zeros(length(V01),1);
        for i = 1:length(V01)
            y1(i) = piston([M,S,V01(i),k,P0,Ta,T01(i)])+ normrnd(0,0.05,1,1);
        end

        %constructing second dataset
        x2 = [V02,T02];
        y2 = zeros(length(V01),1);
        for i = 1:length(V01)
            y2(i) = piston([M,S,V02(i),k,P0,Ta,T02(i)])+ normrnd(0,0.05,1,1);
        end

        %constructing test dataset
        V0ts = linspace(0.002,0.010,sqrt(ntest));
        T0ts = linspace(340,360,sqrt(ntest));
        xtest = combvec(V0ts,T0ts)';
        try
            out = funGP(x1, y1, x2, y2, xtest, confLevel);
        catch
            fprintf(logID,"Error in parameter estimation. Skipping this iteration\n");
            continue;
        end
        if out.differ == true
            n_type1 = n_type1 + 1;        
        end
        fprintf(logID,'Number of replications with Type 1 error: %d out of %d simulations\n',n_type1,j);
        j = j+1;
    end
    fprintf('Estimated Type 1 error: %.3f\n',n_type1/n_simulation);
end