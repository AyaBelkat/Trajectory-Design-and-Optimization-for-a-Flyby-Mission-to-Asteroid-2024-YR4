% Householder iterations

function [x, it] = householder(T, x0, M, eps, max_iter, lambda)

it = 0;
lambda2 = lambda * lambda;
lambda3 = lambda2 * lambda;
buffer = 1e-8;

% Setting both ellipical and hyperbolic solutions for 0 revolutions
% and only elliptical for multiple revolutions
if M > 0
    x0 = max(-1+buffer, min(1-buffer, x0));    % Multi-rev
else
    if x0 < 1
        x0 = max(-1+buffer, min(1-buffer, x0));   % 0 rev elliptical 
    else
        x0 = max(-1+buffer, min(x0, inf));                    % 0 rev hyperbolic
    end
end

err = inf;
while(err > eps && it < max_iter)
    TOF = x2tof(x0, M, lambda);
    [DTdx, DTdx2, DTdx3] = dTdx(x0, TOF, lambda2, lambda3);   % Time of flight derivatives
    delta = TOF - T;
    DT2 = DTdx.^2;
    x_new = x0 - delta * (DT2 - delta * DTdx2/2) / (DTdx * (DT2 - delta * DTdx2) + DTdx3 * delta * delta / 6);
    
    if M > 0
    x_new = max(-1+eps, min(1-eps, x_new));    % Multi-rev
    else
        if x_new < 1
        x_new = max(-1+eps, min(1-eps, x_new));   % 0 rev elliptical 
        else
        x_new = min(x_new, inf);                    % 0 rev hyperbolic
        end
    end

    err = abs(x_new - x0);
    x0 = x_new;
    it = it + 1;
    
end
x = x0; 

end



