% Lancaster-Blanchard time of flight formula for ellipse/hyperbola

function TOF = x2tof2(x, N, lambda)

a = 1/(1 - x * x);

if a > 0   % ellipse

    alfa = 2 * acos(x);
    beta = 2 * asin(sqrt(lambda * lambda/a));

    if lambda < 0
        beta = -beta;
        
    end
    TOF = (a * sqrt(a) * ((alfa - sin(alfa)) -(beta - sin(beta)) + 2 * pi * N))/2;
    
else       % hyperbola
    alfa = 2 * acosh(x);
    beta = 2 * asinh(sqrt(-lambda * lambda / a));

    if lambda < 0
        beta = -beta;
    end
    TOF = (-a * sqrt(-a) * ((beta - sinh(beta)) - (alfa - sinh(alfa))))/2;

end
end
