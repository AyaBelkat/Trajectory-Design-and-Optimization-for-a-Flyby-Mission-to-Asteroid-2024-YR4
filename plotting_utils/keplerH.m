function H = keplerH(M,e)
H = asinh(M/e);                  % good initial guess
for k = 1:30
    f  = e*sinh(H) - H - M;
    fp = e*cosh(H) - 1;
    dH = -f/fp;
    H  = H + dH;
    if abs(dH) < 1e-12, break; end
end
end
