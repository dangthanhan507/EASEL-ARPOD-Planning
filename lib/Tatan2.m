function [thrad] = Tatan2(y,x)
% computes atan2 for tensor variables

pp = pi*Tones(1);
hf = 0.5*Tones(1);

nrmy = sqrt(y'*y);
nrmx = sqrt(x'*x);

thrad = atan(y./x);
temp1 = (y./nrmy);
temp2 = (nrmx - x);
thrad = thrad + temp1.*temp2;

