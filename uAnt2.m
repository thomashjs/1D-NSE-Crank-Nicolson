function out = uAnt2(t,x)
    out = 2^0.5 * exp(1i*(t)).*sech(x);
end