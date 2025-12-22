% hypergeometricF  Compute 2F1(3,1;2.5;z) via series until term < tol

function Q = hypergeometricF(z, tol)

Sj = 1;
Cj = 1;
err = 1;

j = 0;

while err > tol
    Cj1 = Cj * (3 + j) * (1 + j)/(2.5 + j) * z/(j + 1);
    Sj1 = Sj + Cj1;
    err = abs(Cj1);
    Sj = Sj1;
    Cj = Cj1;
    j = j + 1;
end
Q = Sj;
end
