% === Optimize a trajectory with one Gravity Assist (GA) around Venus and one DSM ===

% High-level:
%   - Decision variables x = [t0, TOF, eta1, eta2, theta_b, Rperi]
%   - Objective     : cost_GA_1DSM(x)  (user-defined)
%   - Constraints   : circlcon_GA_Venus_dsm(x) (user-defined)
%   - Local solver  : fmincon (SQP), wrapped in a simple MBH loop
%   - Ephemerides   : SPICE (cspice_*), kernels loaded via meta files

clc; clear;

% --- SPICE setup (adjust paths for your environment) ---------------------
addpath("C:\Users\wbook\Desktop\Dissertation files\Dissertation_codes\mice\mice\src\mice");
addpath("C:\Users\wbook\Desktop\Dissertation files\Dissertation_codes\mice\mice\lib");

cspice_kclear;
cspice_furnsh('C:\Users\wbook\Desktop\Dissertation files\Dissertation_codes\mission_meta.tm');
cspice_furnsh('C:\Users\wbook\Desktop\Dissertation files\Dissertation_codes\mission_equa_data.tm');


% --- Initial epoch guesses ------------------------------
t0_guess = cspice_str2et(datestr(datetime(2030, 1, 17)));     % Departure date guess
tf_guess = cspice_str2et(datestr(datetime(2032, 11, 17)));    % Arrival date guess

Re = 6051.8;   % Venus equatorial radius [km]

% --- Initial decision vector (TOF, eta1, eta2, theta_b, Rperi) ----------
TOF0   = tf_guess - t0_guess;              % time of flight [s]
eta1   = 0.6;                              % fraction at which GA occurs splitting leg 1 and 2
eta2   = 0.4;                              % fraction along leg 2 at which DSM occurs 
t0     = t0_guess;                         % departure epoch [ET seconds]
theta_b = 0.0;                             % B-plane angle [rad]
Rperi   = 2 * Re;                          % flyby periapsis radius [km] (>= safety margin)

x0 = [t0, TOF0, eta1, eta2, theta_b, Rperi];   % Initial guesses vector

% --- Bounds and linear constraints ---------------------------------------

% Time window
et_lb_dep   = cspice_str2et('2028 JAN 01 TDB');    % Earliest departure date from Earth
et_ub_arr   = cspice_str2et('2032 DEC 15 TDB');    % Latest arrival date at the asteroid

% TOF limits (days → seconds)
minTOF_days = 100;                       % days
maxTOF_days = 1200;                      % days
minTOF_sec  = minTOF_days * 86400;       % Convert to second
maxTOF_sec  = maxTOF_days * 86400;       % Convert to second


% Linear inequality: [1 1 0 0 0 0] * x <= b  →  t0 + TOF <= et_ub_arr
A  = [1 1 0 0 0 0]; b = et_ub_arr;   % Arrival date no later than 15 December 2032

% No linear equalities
Aeq = []; beq = [];

% Variable-wise lower/upper bounds:
%  x = [ t0, TOF, eta1, eta2, theta_b, Rperi ]

lb = [ et_lb_dep,               minTOF_sec, 0.05, 0.05, -pi, 1.2*Re ];
ub = [ et_ub_arr - minTOF_sec, et_ub_arr + t0, 0.95, 0.95,  pi, 10.0*Re ];

% --- Objective and nonlinear constraints --------
nonlcon = @circlcon_GA_Venus_dsm;       % Enforces nonlinear constraints
J       = @cost_GA_Venus_1DSM;          % Returns scalar cost 

% Quick constraint check at the initial guess
[c0, ceq0] = circlcon_GA_dsm(x0);
disp('Initial constraint values:'); disp(c0(:).');

% --- Local optimizer settings (SQP) --------------------------------------
options = optimoptions('fmincon','Display','iter','Algorithm','sqp', ...
    'UseParallel',false,'ConstraintTolerance',1e-7,'ScaleProblem',true);

% --- Local solve (starting from x0) --------------------------------------
[X, Jmin] = fmincon(J, x0, A, b, Aeq, beq, lb, ub, nonlcon, options); 

check_more(X);    % optional routine for testing

% Report local results
optimized_departure_date = cspice_et2utc(X(1), 'C', 3);
optimized_arrival_date  = cspice_et2utc(X(1)+X(2), 'C', 3);
eta1_opt  = X(3);
eta2_opt  = X(4);
theta_b_opt = X(5);
Rperi_opt   = X(6);

disp('optimized_departure_date'); disp(optimized_departure_date)
disp('optimized_arrival_date');  disp(optimized_arrival_date)
disp('eta1');  disp(eta1_opt)
disp('eta2');  disp(eta2_opt)
disp('theta_b'); disp(theta_b_opt)
disp('Rperi [km]'); disp(Rperi_opt)

% === MONOTONIC BASIN HOPPING (MBH) =======================================
% Simple global-local hybrid: random perturbation → local SQP → keep best

numIterations = 60;    % Number of iterations
bestSolution = X;      % Initial guess from local optimizer results
bestCost     = Jmin;   % Best cost function from local optimizer results

% Utility: wrap any angle to (-pi, pi]
wrapPi = @(theta_b) atan2(sin(theta_b), cos(theta_b));   % maps any angle to (−π, π]

%rng(42,'twister');   % make the run deterministic
for iter = 1:numIterations
    % --- Random perturbation around current best ----------
    perturbVec = [ 86400*randn, ...         % t0 shift ~ ±1 day
                   5*86400*randn, ...       % TOF shift ~ ±5 days
                   0.05*randn, ...          % eta1 shift
                   0.05*randn, ...          % eta2 shift
                   0.2*randn,  ...          % theta_b shift [rad]
                   0.2*Re*randn ];          % Rperi shift [km] ~ ±0.2 Re

    newGuess = bestSolution + perturbVec;   % Introducing perturbations
    newGuess = max(newGuess, lb);           % Enforcing bounds
    newGuess = min(newGuess, ub);
    newGuess(5) = wrapPi(newGuess(5));

    % Ensure arrival cap not violated while keeping ≥ min TOF
    arr_cap = et_ub_arr;
    if newGuess(1) + newGuess(2) > arr_cap
    newGuess(2) = max(minTOF_sec, arr_cap - newGuess(1));  % keep ≥ min TOF
    end

    % --- Local solve from perturbed seed --------------------------------
    [X_trial, J_trial] = fmincon(J, newGuess, A, b, Aeq, beq, lb, ub, nonlcon, options);

    % --- Accept improvement ---------------------------------------------
    if J_trial < bestCost
        bestSolution = X_trial;   % <- use optimized local result
        bestCost     = J_trial;
        fprintf('Iteration %d improved: cost = %.3f\n', iter, bestCost);
        check_more(bestSolution);   % optional diagnostics for the new best
    end
end


% --- Final report --------------------------------------------------------
optimized_departure_date = cspice_et2utc(bestSolution(1), 'C', 3);
optimized_arrival_date = cspice_et2utc(bestSolution(1) + bestSolution(2), 'C', 3);


disp('optimized_departure_date');  disp(optimized_departure_date)
disp('optimized_arrival_date');   disp(optimized_arrival_date)

disp('Best eta1 achieved');    disp(bestSolution(3));
disp('Best eta2 achieved:'); disp(bestSolution(4));
disp('Best theta_b achieved:'); disp(bestSolution(5));
disp('Best Rperi achieved:'); disp(bestSolution(6));

disp('Best cost achieved:');
disp(bestCost);

check_more(bestSolution);