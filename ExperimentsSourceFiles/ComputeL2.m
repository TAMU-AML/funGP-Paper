function[value] = ComputeL2(fun,fmin,fmax)
    if length(fmin) == 2
        integrand = @(x,y)sqFun2(fun,x,y);
        value = sqrt(integral2(integrand,fmin(1),fmax(1),fmin(2),fmax(2)));
    elseif length(fmin) == 1
        integrand = @(x)sqFun1(fun,x);
        value = sqrt(integral(integrand,fmin,fmax));
    else
        error("The method works only for 1-dim or 2-dim functions");
    end
             
     function[vector] = sqFun1(inputfun, x)
        vector = zeros(1,length(x));
        for i = 1:length(x)
            vector(i) = inputfun(x(i))^2;
        end
     end
     
     function[vector] = sqFun2(inputfun, x,y)
        vector = zeros(size(x));
        for i = 1:length(x)
            vector(i) = inputfun(x(i),y(i))^2;
        end
     end
end
