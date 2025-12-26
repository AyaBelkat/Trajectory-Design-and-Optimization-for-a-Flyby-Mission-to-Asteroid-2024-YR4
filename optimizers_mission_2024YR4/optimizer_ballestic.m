% =========================================================================
%   OPTIMIZER_BALlISTIC.M
%
%   Purpose:
%       This script sets up and solves a preliminary Earth-to-NEO transfer 
%       optimisation problem using a purely ballistic trajectory model. 
%       It frames the mission design as a constrained nonlinear program 
%       (NLP) with decision variables representing departure and arrival 
%       epochs.
%
%   Approach:
%       - The spacecraft dynamics are modelled with a two-body Sun-centred 
%         approximation (patched-conics at departure and arrival).
%       - The cost function (cost_function_ballistic_velocity) evaluates the 
%         characteristic energy (C3) or equivalent v∞ velocity at departure, 
%         to be minimized for a feasible trajectory.
%       - Nonlinear constraints (circlcon) enforce mission requirements such 
%         as minimum time of flight, planetary geometry conditions, and 
%         arrival deadline compatibility.
%       - The optimisation is handled by fmincon (SQP), with bounds ensuring 
%         realistic departure and arrival dates within the mission window.
%
%   Decision variables:
%       x = [ t0, tf ]
%           t0  : Departure epoch [ET seconds past J2000]
%           tf  : Arrival epoch   [ET seconds past J2000]
%
%   Constraints:
%       - t0 ≥ lower bound (2028-01-01)
%       - tf ≤ upper bound (2032-12-15)
%       - tf - t0 ≥ minimum time of flight (100 days)
%       - circlcon(x) imposes additional dynamical and geometry conditions
%
%   Output:
%       The optimiser returns optimal departure and arrival epochs, printed 
%       both in ET seconds and as UTC calendar dates. These serve as 
%       baseline solutions for subsequent higher-fidelity trajectory design.
%
%   Context:
%       This ballistic optimiser provides a first-cut solution for the 
%       mission trajectory design process. It is deliberately simplified to 
%       highlight feasibility regions before extending to models that include 
%       midcourse DSMs or gravity assists ( MGA/MGADSM frameworks) 
clc; clear;
startup;

% === MICE (SPICE) setup ====================================================
cspice_kclear;
load_kernels;

% === Initial guesses (calendar → SPICE ET) =================================
% Seed departure/arrival epochs for the local optimiser (fmincon).
t0_guess = cspice_str2et(datestr(datetime(2032, 2, 17)));     % First guess for departure date
tf_guess = cspice_str2et(datestr(datetime(2032, 11, 17)));    % First guess for arrival date

% === Hard window bounds (design constraints) ==============================
% Lower bound on departure, upper bound on arrival.
et_lb_dep = cspice_str2et('2028 JAN 01 TDB');   % lower bound for departure
et_ub_arr = cspice_str2et('2032 DEC 15 TDB');   % upper bound for arrival

% Minimum TOF requirement (robustness/mission feasibility)
minTOF_days = 1;                                % minimum time of flight 
minTOF_sec  = minTOF_days*86400;

% Decision vector (ballistic case): x = [t0, tf]
x0 = [t0_guess, tf_guess];

% === Linear constraints ====================================================
% Inequality: A*x <= b  → here ensures arrival before the deadline (tf ≤ et_ub_arr)
A = [0, 1];
b = et_ub_arr;

% No linear equalities
Aeq = [];
beq = [];

% Variable bounds:
%  t0 ∈ [et_lb_dep, et_ub_arr - minTOF]
%  tf ∈ [t0 + minTOF, et_ub_arr - 15 d]  (keep 15 d headroom before hard deadline)
lb = [et_lb_dep,  et_lb_dep + minTOF_sec];    % tf can't be before t0+minTOF
ub = [et_ub_arr - minTOF_sec,  et_ub_arr - 15*86400];    % leave room for minTOF before deadline

% === Nonlinear constraints and objective ==================================
% circlcon : user-defined nonlinear feasibility (e.g., geometry, segment minima)
% J        : user-defined scalar cost (e.g., launch energy or Δv proxy)
nonlcon = @circlcon;
J = @cost_function_ballistic_velocity;

% Quick check of initial constraint vector at x0
[c0, ceq0] = circlcon(x0);
disp('Initial constraint values:')
disp(c0)

% === Local optimiser (SQP) =================================================
% Display iterations; use moderate constraint tolerance; single-threaded.
options = optimoptions('fmincon','Display','iter','Algorithm','sqp', 'UseParallel',false, 'ConstraintTolerance', 1e-7);

% Solve the constrained optimisation starting from x0
[X, Jmin] = fmincon(J, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);

% === Report optimised epochs ==============================================
% Convert optimised ET back to human-readable UTC strings.
optimized_departure_date = cspice_et2utc(X(1), 'C', 3);
optimized_arrival_date = cspice_et2utc(X(2), 'C', 3);

disp('optimized departure date');   disp(optimized_departure_date)
disp('optimized arrival date');     disp(optimized_arrival_date)

disp('Best cost');       disp(Jmin)




