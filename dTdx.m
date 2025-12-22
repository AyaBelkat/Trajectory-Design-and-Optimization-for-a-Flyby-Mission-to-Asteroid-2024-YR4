% dTdx funtion to find the time of flight derivatives

function [DTdx, DTdx2, DTdx3] = dTdx(x, T, lambda2, lambda3)

tmp = 1/(1 - x^2);
y = sqrt(1 - lambda2*(1 - x^2));

DTdx = tmp * (3 * T * x - 2 + 2 * lambda3 * (x/y));
DTdx2 = tmp * (3 * T + 5 * x * DTdx + 2 * (1 - lambda2) * (lambda3/y^3));
DTdx3 = tmp * (7 * x * DTdx2 + 8 * DTdx - 6 * (1 - lambda2) * lambda2 * lambda3 * (x/y^5));

end