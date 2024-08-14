function u0 = NLSic(x)
u0 = 2^0.5 * exp(0.5i*x).*sech(x);
end