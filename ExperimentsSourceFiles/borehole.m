function [y] = borehole(xx)
rw = xx(1);
r  = xx(2);
Tu = xx(3);
Hu = xx(4);
Tl = xx(5);
Hl = xx(6);
L  = xx(7);
Kw = xx(8);

frac1 = 2 * pi * Tu * (Hu-Hl);

frac2a = 2*L*Tu / (log(r/rw)*rw^2*Kw);
frac2b = Tu / Tl;
frac2 = log(r/rw) * (1+frac2a+frac2b);

y = frac1 / frac2;

end