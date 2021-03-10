%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DESCRIPTION:
%   Function to compute the correlation matrix for ARD squared exponential kernel
%   INPUTS:
%       x1: Input variable matrix or vector for first set of points
%       x2: Input variable matrix or vector for second set of points
%       theta: lengthscale for each input variable
%   OUTPUT:
%       correlMat: correlation matrix between x1 and x2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[correlMat] = computeCorrelMat(x1,x2, theta)
    nrows = size(x1,1); %number of rows for the correlation matrix
    ncols = size(x2,1); %number of columns for the correlation matrix
    ncov = size(x1,2); %number of input variables
    correlMat = zeros(nrows,ncols); %initialize the correlation matrix
    if ncov > 1
        %compute the squared distance matrix by adding one input variable at a time
        for i = 1:ncov
            correlMat = correlMat + (((x1(:,i) - (x2(:,i)')).^2)/(theta(i)^2));
        end
    else         
        correlMat =  (((x1 - (x2')).^2)/(theta^2));         
    end
    correlMat = exp(-0.5*correlMat); %compute correlation mat from squared dist mat
end