%% COE2RV function to convert from classical orbit elements to SCI position and velocity coordinates %%

function [X] = COE2RV(coe, mu)

% Retrieving input
a = coe(1,:);
e = coe(2,:);
i = coe(3,:);
RAAN = coe(4,:);
argP = coe(5,:);
theta = coe(6,:);

% Rotation Matrix R3(RAAN)

R3_RAAN = [cos(RAAN), sin(RAAN), 0;
    -sin(RAAN), cos(RAAN), 0;
    0, 0, 1];

% Rotation Matrix R1(i)

R1_i = [1, 0, 0;
    0, cos(i), sin(i);
    0, -sin(i), cos(i)];

% Rotation Matrix R3(omega)

R3_argP = [cos(argP), sin(argP), 0;
    -sin(argP), cos(argP), 0;
    0, 0, 1];

% Direction cosine matrix IP

IP = R3_argP*R1_i*R3_RAAN;

r = (a*(1-e^2))/(1+e.*cos(theta));    % calculating position using orbit elements

h = sqrt(mu*a*(1-e^2));                          % calculating angular momentum

% Poition and velocity coordinates

r_coe = [r*cos(theta); r*sin(theta); 0];
v_coe = [-(mu / h) * sin(theta); (mu / h) * (e + cos(theta)); 0];

r_ECI = IP'*r_coe;      % converting to SCI position
v_ECI = IP'*v_coe;      % converting to SCI velocity


% Output

X = [r_ECI; v_ECI];
