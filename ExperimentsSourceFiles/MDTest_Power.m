function[out] = MDTest_Power(MSq,confLevel,sigmaSq1,sigmaSq2,n1, n2)
    q = n2/n1;
    xiSq0 = ((q+1)*(sigmaSq1^2)) + ((1+(1/q))*(sigmaSq2^2)) + (2*sigmaSq1*sigmaSq2*(1+(1/q)));
    xiSqM = ((q+1)*(sigmaSq1^2 + (4*sigmaSq1*MSq))) + ((1+(1/q))*(sigmaSq2^2 + (4*sigmaSq2*MSq))) + (2*sigmaSq1*sigmaSq2*(1+(1/q)));
    out = 1-normcdf(((norminv(confLevel)*sqrt(xiSq0)) - (sqrt(n1+n2)*MSq))/sqrt(xiSqM));
end

