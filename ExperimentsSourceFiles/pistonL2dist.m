%L2 distance for piston function
warning('off') % turning off tolerance warning
M  = 45;
S  = 0.010;
%x = V0;
k  = 2000;
P0 = 100000;
Ta = 292;
%y = T0;

fmin = [0.002, 340];
fmax = [0.01, 360];

funDiff = @(x,y) 100*(piston([M,S,x,k,P0,Ta,y]) - piston([M,S,x,k+500,P0,Ta,y]));
fun = @(x,y) 100*piston([M,S,x,k,P0,Ta,y]);
L2Dist = ComputeL2(funDiff,fmin,fmax);
L2NormF = ComputeL2(fun,fmin,fmax);
percentDiff = L2Dist*100/L2NormF;
fprintf('l2 distance for piston: %f \n',percentDiff);
clear