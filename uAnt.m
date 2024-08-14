function out = uAnt(t,x)
    out = 2^0.5 * exp(1i*(0.5*x+t)).*sech(x);
end