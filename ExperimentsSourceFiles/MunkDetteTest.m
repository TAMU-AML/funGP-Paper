function[out] = MunkDetteTest(x1, y1, x2, y2, confLevel,rangeX)
    if rangeX(1) ~= 0 && rangeX(2) ~= 1 
        error('The domain of the function should be 0 to 1');
    end
   [x1_sorted, sort_index1] = sort(x1);
   y1_sorted = y1(sort_index1);
   length_y1 = length(y1_sorted);
   
   [x2_sorted, sort_index2] = sort(x2);
   y2_sorted = y2(sort_index2);
   length_y2 = length(y2_sorted);
   
   q = length_y2/length_y1; %now handling different lengths of y's
   
   sigmaSq1 = sum((y1_sorted(2:end) - y1_sorted(1:(end-1))).^2)/(2*(length_y1 -1));
   sigmaSq2 = sum((y2_sorted(2:end) - y2_sorted(1:(end-1))).^2)/(2*(length_y2 -1));
   xiSq = ((1+q)*(sigmaSq1^2))+((1+(1/q))*(sigmaSq2^2))+(2*sigmaSq1*sigmaSq2*(1+(1/q)));  
   
   y1_0 = y1(1);
   x1_0 = rangeX(1);
   y1_max = y1(end);
   x1_max = rangeX(2);
   y1 = [y1_0;y1_sorted;y1_max];
   x1 = [x1_0;x1_sorted;x1_max];
   
   y2_0 = y2(1);
   x2_0 = rangeX(1);
   y2_max = y2(end);
   x2_max = rangeX(2);
   y2 = [y2_0;y2_sorted;y2_max];
   x2 = [x2_0;x2_sorted;x2_max];
   
   
   mHatSq = 0;
   for i = 1:(length_y1+1)
       for j = 1:(length_y2+1)
           lambda_ij =  max(min(x1(i+1),x2(j+1)) - max(x1(i),x2(j)), 0);
           mHatSq = mHatSq + (lambda_ij*(y1(i+1) - y2(j+1))*(y1(i)-y2(j)));
       end
   end
   
   test_stat = sqrt(length_y1 + length_y2)*mHatSq/sqrt(xiSq);
   
   threshold = norminv(confLevel);
  
   if test_stat > threshold
       out.differ = true;
   else
       out.differ = false;
   end
   
   out.test_stat = test_stat;
   out.threshold = threshold;
  
end

   