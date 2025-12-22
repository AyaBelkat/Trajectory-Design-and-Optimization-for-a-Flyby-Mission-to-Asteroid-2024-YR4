% === Lambert solver main function =========================================
% Purpose:
%   - Given two position vectors r1, r2, a time of flight (TOF) and the
%     gravitational parameter mu, solve Lambert's orbital boundary value
%     problem.
%   - Returns the possible departure and arrival velocity vectors (v1, v2).
% References: References: Lancaster–Blanchard formulation, Izzo, D. (2015). "Revised Lambert's problem solver".
% =========================================================================


function [v1, v2] = lambert_solver(r1, r2, TOF, mu)

M = 1;      % Number of revolutions considered
cw = 0;     % Motion type: prograde cw=0, retrograde cw=1

% --- Geometry setup -------------------------------------------------------
c = r1 - r2;                   % Chord vector
R1 = norm(r1);                 % Magnitude of departure position
R2 = norm(r2);                 % Magnitude of arrival position
C = norm(c);                   % Chord length
s = 0.5 * (R1 + R2 + C);       % Semi-perimeter of triangle (r1, r2, chord)

% --- Local reference frame ------------------------------------------------
ir1 = r1/R1;              % Unit vector along r1
ir2 = r2/R2;              % Unit vector along r2
ih = cross(ir1, ir2);     % Angular momentum direction
ih = ih/norm(ih);         % Normalized

if ih(3) == 0
    % Guard against singular geometry: trajectory undefined if ih has no z-component
    error('The angular momentum vector has no z component -> impossible to define cw or ccw');
end

% --- Define lambda parameter (Lancaster–Blanchard variable) ---------------
lambda2 = 1 - C/s;          % Square of lambda
lambda = sqrt(lambda2);


% --- Correct for transfer angles > 180° -----------------------------------
% --- Setting up the frame where the velocities will be constructed -------
if ih(3) < 0
    lambda = -lambda;       % Adjust lambda for transfer angle greater than 180°
    it1 = cross(ir1, ih);  
    it2 = cross(ir2, ih);  
    
else
    it1 = cross(ih, ir1);
    it2 = cross(ih, ir2);
    
end

% --- Retrograde motion correction -----------------------------------------
if cw == 1
    it1 = -it1;         % Adjust it1 for retrograde motion
    it2 = -it2;         % Adjust it2 for retrograde motion
    lambda = -lambda;   % Adjust lambda for retrograde motion
end

% --- Non-dimensionalize time of flight ------------------------------------
T = sqrt((2*mu)/s^3) * TOF;    % Dimensionless TOF (scaled by s and mu)

% --- Auxiliary values for iteration ---------------------------------------
lambda3 = lambda * lambda2;
M_max = floor(T/pi);     % Maximum feasible revolutions

% Bounds for TOF
T00 = acos(lambda) + lambda * sqrt(1 - lambda2);   % Minimum TOF, zero rev
T0 = T00 + M_max * pi;           % Multi-rev TOF bound
T1 = (2/3) * (1 - lambda3);      % Upper bound for small TOF

% --- Pre-check to cap multi-rev -------------------------------
err = 1;  
if M_max > 0
    if T < T0
        x_old = 0;
        T_min = T0;
        for i = 1:13
            [DTdx, DTdx2, DTdx3] = dTdx(x_old, T_min, lambda2, lambda3);
            if DTdx ~= 0
                x_new = x_old - (DTdx * DTdx2)/(DTdx2 * DTdx2 - (DTdx * DTdx3)/2);
            end
           
            err = abs(x_old - x_new);
            if err < 1e-13
                break;
            end

            T_min = x2tof(T_min, x_new, M_max);

            x_old = x_new;
            
        end
        if T_min > T
            M_max = M_max - 1;
        end

    end
end

M_max = min(M, M_max);     % Enforcing the upper limit of the number of revolutions

% --- Initial guess generation -------------------------------
if T >= T00
    x0 = power((T00/T), 2/3) - 1;    
elseif T <= T1
    x0 = (5/2) * (T1 * (T1 -T))/(T*(1 - lambda2 * lambda3)) + 1;
else
    x0 = power(T00/T, log2(T1/T00)) - 1;   
end

% --- Storage for solutions ------------------------------------------------
n_solutions = 2*M_max+1;              
x = zeros(n_solutions, 1);
it = zeros(n_solutions, 1);

% --- Householder iterations for 0 rev  -----------------
[x(1), it(1)] = householder(T, x0, 0, 1e-5, 15, lambda);


% --- Householder iterations for multiple revolutions----------------
for i = 1:M_max
    % Left householder iterations
    tmp = power((i * pi + pi)/(8 * T), 2/3);
    x0r = (tmp - 1)/(tmp + 1);
    [x(2 * i), it(2 * i)] = householder(T, x0r, i, 1e-8, 15, lambda);

    % Right householder iterations
    tmp = power((8 * T)/(i * pi), 2/3);
    x0l = (tmp - 1)/(tmp + 1);
    [x(2 * i + 1), it(2 * i + 1)] = householder(T, x0l, i, 1e-8, 15, lambda);

end

% --- Velocity reconstruction ----------------------------------------------
gamma = sqrt((mu * s)/2);
rho = (R1 - R2)/C;
sigma = sqrt(1 - rho^2);

for i = 1:length(x)
    y = sqrt(1 - lambda2 + lambda2 * x(i) * x(i));      % Lancaster–Blanchard relation
    vr1 = gamma * ((lambda * y - x(i)) - rho * (lambda * y + x(i)))/R1;
    vr2 = -gamma * ((lambda * y - x(i)) + rho * (lambda * y + x(i)))/R2;
    vt = gamma * sigma * (y +lambda * x(i));      
    vt1 = vt/R1;                  % Tangential components
    vt2 = vt/R2;
    
    v1(i, :) = vr1 * ir1 + vt1 * it1;      % Departure velocity solution(s)
    v2(i, :) = vr2 * ir2 + vt2 * it2;      % Arrival velocity solution(s)
end
