%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DESCRIPTION:
%   Function to compute the confidence band on difference of two GPs over a finite grid.
%   INPUTS:
%       K: Covariance matrix of the difference process over a finite grid.
%       confLevel: Significance level
%   OUTPUT:
%       band: positive vector of upper bounds for the confidence band. 
%              lower band can be computed as negative of the upper band.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[band] = computeConfBand(K,confLevel,option)
   
    K = (K + K')/2; %to make sure that the matrix is numerically symmetric
    [V,D] = eig(K); %eigen decomposition
    [d, sorted_index] = sort(diag(D),'descend'); %sort eigenvalues in decreasing order
    Ds = D(sorted_index,sorted_index); %sorted eigenvalue diagonal matrix
    Vs = V(:,sorted_index); %sorted eigenvectors
    %print largest and smallest eigenvalues
    %fprintf('Largest eigenvalue: %8.10f \n',max(d));
    %fprintf('Smallest eigenvalue: %8.10f \n',min(d));

    if ismember('rel_threshold',fieldnames(option.truncation))
        rel_threshold = option.truncation.rel_threshold; %threshold to compute truncation number
        %find the number of eigenvalues greater than threshold*max(d)
        if (rel_threshold*max(d)) > min(d)
            m = find(d<rel_threshold*max(d),1);
        else
            m = size(K,1);
        end
    elseif ismember('abs_threshold',fieldnames(option.truncation))
        abs_threshold = option.truncation.abs_threshold;
        if abs_threshold > min(d)
            m = find(d<abs_threshold,1);
        else
            m = size(K,1);
        end
    else
        if option.truncation.trunc_number <= size(K,1)
            m = option.truncation.trunc_number;
        else
            m = size(K,1);
        end
    end



    %print the tructation number and the smallest eigenvalue selected
    %fprintf('Number of eigenvalues selected: %d \n',m);
    %fprintf('Smallest eigenvalue selected: %f \n',d(m));


    R = sqrt(chi2inv(confLevel,m)); %compute the radius of coverage region with probability confLevel
    t_samples = option.zsamples; %number of samples to be used in simulating the band at each point
    Z = zeros(t_samples,m); %matrix to strore the simulated standard normal random variables. 
    
    %sample random vectors till the number of samples reach t_samples
    t = 1;
    while t <= t_samples
        Z(t,:) = normrnd(0,1,m,1);
        if norm(Z(t,:))<=R
         t = t+1;
        end
    end

    Gx = Vs(:,1:m)*sqrt(Ds(1:m,1:m))*(Z'); %compute the simulated sample paths 
    band =  max(abs(Gx),[],2); %take the absolute maximum as the upper band.
end