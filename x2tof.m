% x to time of flight transformation

function TOF = x2tof(x, N, lambda)

battin  = 0.01;
lagrange = 0.2;
dist = abs(x - 1);

% Langrange region |x - 1| < 0.2 but not very small

if dist < lagrange && dist > battin
    TOF = x2tof2(x, N, lambda);
end
K = lambda * lambda;
E = x * x - 1;
rho = abs(E);
z = sqrt(1 + K * E);

% Battin region |x - 1| < 0.01

if dist < battin 
    eta = z - lambda * x;
    S1 = 0.5 * (1 - lambda - x * eta);
    Q = hypergeometricF(S1, 1e-11);
    Q = 4/3 * Q;
    TOF = (eta * eta * eta * Q + 4 * lambda * eta)/2 + (N*pi)/(rho^1.5);
   
else           % Lancaster region for all other values of |x - 1|
    y = sqrt(rho);
    g = x * z - lambda * E;
   
    if E < 0
        l = acos(g);
        d = N * pi + l;
    else
        f = y * (z - lambda * x);
        d = log(f + g);
    end
    TOF = (x - lambda * z - d / y)/E;
    
end
end
