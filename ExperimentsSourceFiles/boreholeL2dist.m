%L2 distance for borehole function
warning('off') % turning off tolerance warning
%x = rw;
%y  = r;
Tu = 78000;
Hu = 1050;
Tl = 84;
Hl = 760;
L  = 1400;
Kw = 11000;

fmin = [0.05, 100];
fmax = [0.15, 50000];

funDiff = @(x,y) borehole([x,y,Tu,Hu,Tl,Hl,L,Kw]) - borehole([x,y,Tu,Hu,Tl,Hl,L+50,Kw]);
fun = @(x,y) borehole([x,y,Tu,Hu,Tl,Hl,L,Kw]);
L2Dist = ComputeL2(funDiff,fmin,fmax);
L2NormF = ComputeL2(fun,fmin,fmax);
percentDiff = L2Dist*100/L2NormF;
fprintf('l2 distance for borehole: %f \n',percentDiff);
clear