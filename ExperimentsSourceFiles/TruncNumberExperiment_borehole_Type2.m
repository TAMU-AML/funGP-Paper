%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DESCRIPTION:
%   main function for funGP truncation number experiment for Type 2 error for borehole function 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[] = TruncNumberExperiment_borehole_Type2(trunc_number)
   
    logfile = strcat('TruncNumber_',num2str(trunc_number),'_borehole_Type2.log');
    logID = fopen(logfile,'w');
    n_simulation = 1000; %number of simualtions
    n_type2 = 0 ; %variable to store the number of simulations with type 2 error
    confLevel = 0.95; %significance level
    ntrain = 1000;
    ntest = 50*50;
    option.truncation.trunc_number = trunc_number; %specifying the truncation number for calculating the confBand

    fprintf('Truncation number experiment for Type 2 error using borehole function\n')
    fprintf('Size of train data: %d\n', ntrain);
    fprintf('Size of test data: %d\n', ntest);
    fprintf('Truncation Number: %d\n', trunc_number);
    fprintf('Nominal level of the test: %.7f\n', confLevel);

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
        rw01 = 0.05 + (0.15 - 0.05).*rand(ntrain,1);
        r01  = 100 + (50000 - 100).*rand(ntrain,1);
        rw02 = 0.05 + (0.15 - 0.05).*rand(ntrain,1);
        r02  = 100 + (50000 - 100).*rand(ntrain,1);

        %constructing first dataset
        x1 = [rw01,r01];
        y1 = zeros(length(rw01),1);
        for i = 1:length(rw01)
            y1(i) = borehole([rw01(i),r01(i),Tu,Hu,Tl,Hl,L,Kw])+ normrnd(0,10,1,1);
        end

        %constructing second dataset
        x2 = [rw02,r02];
        y2 = zeros(length(rw01),1);
        for i = 1:length(rw01)
            y2(i) = borehole([rw02(i),r02(i),Tu,Hu,Tl,Hl,Lmod,Kw])+ normrnd(0,10,1,1);
        end

        %constructing test dataset
        rw0ts = linspace(0.05,0.15,sqrt(ntest));
        r0ts = linspace(100,50000,sqrt(ntest));
        xtest = combvec(rw0ts,r0ts)';
        try
            out = funGP(x1, y1, x2, y2, xtest, confLevel, option);
        catch
            fprintf(logID,"Error in parameter estimation. Skipping this iteration\n");
            continue;
        end
        if out.differ == false
            n_type2 = n_type2 + 1;        
        end
        fprintf(logID,'Number of replications with Type 2 error: %d out of %d simulations\n',n_type2,j);
        j = j+1;
    end
    fprintf('Estimated Type 2 error: %.3f\n',n_type2/n_simulation);
end