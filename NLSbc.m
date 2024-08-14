function [pl,ql,pr,qr] = NLSbc(xl,ul,xr,ur,t)
pl = ul-2^0.5*exp(1i*t)*sech(xl); %ignored by solver since m=1
ql = 0; %ignored by solver since m=1
pr = ur-2^0.5*exp(1i*t)*sech(xr);
qr = 0; 
end