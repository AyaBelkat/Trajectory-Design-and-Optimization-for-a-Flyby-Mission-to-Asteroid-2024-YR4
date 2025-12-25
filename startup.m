% =========================================================================
% Startup Script: 2024 YR4 Mission Analysis & Optimization
% =========================================================================
% Purpose:
%   - Automatically configures the MATLAB search path for all subfolders.
%   - Enables modular access to /core_solvers from any directory.
%   - Verifies the availability of the SPICE (Mice) Toolkit.
% =========================================================================

fprintf('Initializing Mission Design Environment...\n');

% 1. Add ALL folders and subfolders to the MATLAB path
% genpath(pwd) recursively finds every folder in your project structure.
% This allows main scripts to "see" functions inside /core_solvers,
% /direct_transfer, and /mgadsm_optimization.
addpath(genpath(pwd));
fprintf('-> Recursive path integration complete.\n');

% 2. Check for SPICE (Mice) Toolkit
% Essential for high-fidelity planetary states and J2000 ephemeris data.
% If this fails, cspice_str2et or cspice_spkezr calls will error out.
% 2. Check for SPICE (Mice) Toolkit
if exist('cspice_spkezr', 'file') ~= 3
    % If not found, try adding a standard relative path 
    % (Useful if you keep 'mice' in your project root)
    if isfolder('mice')
        addpath(genpath('mice'));
        fprintf('-> SPICE Toolkit: Found in local /mice folder.\n');
    else
        warning('SPICE (Mice) Toolkit not detected.');
        fprintf('   Action: Ensure the Mice toolkit is in your MATLAB path.\n');
    end
else
    fprintf('-> SPICE Toolkit: OK (Detected in Path)\n');
end

% 3. Verify Core Solver and MGADSM Engine Existence
if exist('lambert_solver.m', 'file') && exist('cost_GA_1DSM.m', 'file')
    fprintf('-> Core Solvers & MGADSM Engines: Detected\n');
else
    warning('One or more mission engines are missing. Check your folder structure.');
end

fprintf('------------------------------------------------------------\n');
fprintf('Environment Ready.\n');
fprintf('Entry point: Run /validation/earth_mars_test.m to verify solver.\n');
