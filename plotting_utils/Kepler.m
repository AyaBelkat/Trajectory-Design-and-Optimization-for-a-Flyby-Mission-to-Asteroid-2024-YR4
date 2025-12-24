%% Kepler function to implement Newton's method for Kepler's equation%%

function [E] = Kepler(e, M, tol)

% retrieve initial guess
E_n = M;

for ii = 1:10

    fE_n = E_n - e*sin(E_n) - M;   % evaluate function in E_n
    fpE_n = 1 - e*cos(E_n);        % evaluate first derivative in E_n

    % Print output
    %fprintf('E_n: %f,  |f(E_n)|: %e\n', E_n, abs(fE_n));

    % check for convergence using IF statement
    if(abs(fE_n) < tol)
        break;
    else
        %Newton's method
        dE_n = -fE_n/fpE_n;        % delta E_n
        E_n = E_n + dE_n;          % update value of E_n
    end
end

% output variables
E = E_n;


  
