%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DESCRIPTION:
%   main function for funGP simulation for Type 2 error for piston function 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_simulation = 1000; %number of simualtions
n_type2 = 0 ; %variable to store the number of simulations with type 2 error
confLevel = 0.95; %significance level
ntrain = 1000; %number of training points
ntest = 50*50; %number of testpoints

fprintf("Simulation for function: %s\n","piston");
fprintf("Type of test: %s\n","Type 2");
fprintf("Sample size for each dataset: %d\n",ntrain);
fprintf("Size of test grid: %d\n",ntest);
fprintf("Confidence level: %.7f\n",confLevel);
fprintf("Number of simulations: %d\n",n_simulation);

logfile = "piston_Type2.log";
logID = fopen(logfile,'w');
fprintf('Storing computation log in file: %s\n',logfile);

%Fixed input values 
M  = 45;
S  = 0.010;
k  = 2000;
P0 = 100000;
Ta = 292;
kmod = k+500;

rng(2); %setting the seed 
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
        y2(i) = piston([M,S,V02(i),kmod,P0,Ta,T02(i)])+ normrnd(0,0.05,1,1);
    end

    %constructing test dataset
    V0ts = linspace(0.002,0.010,floor(sqrt(ntest)));
    T0ts = linspace(340,360,floor(sqrt(ntest)));
    xtest = combvec(V0ts,T0ts)';
    try
        out = funGP(x1, y1, x2, y2, xtest, confLevel);
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
