% =========================================================================
%   PLOT_BALLISTIC_E2AST.M
%
%   Purpose:
%       Plot a purely ballistic (single Lambert) heliocentric trajectory
%       from Earth at t0 to asteroid 2024 YR4 at tf. Mirrors the visual
%       style of the MGADSM plotter (continuous real-scale arc + key arrows).
%
%   Dependencies:
%       - MICE (SPICE) installed and meta-kernel furnished
%       - lambert_solver(r1, r2, TOF, muS) available (Izzo-based)
%
%   Notes:
%       - Branch selection: chooses the Lambert branch that minimizes
%         Earth-escape injection Δv from a reference GTO (rp=Re+250 km,
%         ra=Re+22500 km). Switch strategy easily to min |v1| if desired.
% =========================================================================
clear; clc;
startup;

% ---------- USER INPUT: epochs (UTC calendar) ----------
t0_cal = datetime(2032,  2, 17);   % example departure
tf_cal = datetime(2032,  11, 17);   % example arrival

% ---------- SPICE setup ----------
cspice_kclear;
load_kernels;

% ---------- Constants ----------
muS = 1.32712440018e11;   % Sun GM [km^3/s^2]
muE = 398600;             % Earth GM [km^3/s^2]
Re  = 6378.137;           % Earth radius [km]

% GTO reference (Ariane-6-like) for Δv_inj proxy
rp  = Re + 250;           % km
ra  = Re + 22500;         % km
aGTO = 0.5*(rp + ra);     % km

% ---------- Convert epochs to ET ----------
t0 = cspice_str2et(datestr(t0_cal));
tf = cspice_str2et(datestr(tf_cal));
TOF = tf - t0;

% ---------- Ephemerides (Sun/J2000) ----------
E0  = cspice_spkezr('Earth',        t0, 'J2000','NONE','Sun');
AST = cspice_spkezr('54509621',     tf, 'J2000','NONE','Sun');  % 2024 YR4

rE0  = E0(1:3);  vE0  = E0(4:6);
rAST = AST(1:3); vAST = AST(4:6);

% ---------- Lambert solve (may have multiple branches) ----------
[v1, v2] = lambert_solver(rE0, rAST, TOF, muS);

% If multiple Lambert branches exist, select the lowest |v1| branch
if size(v1,1) > 1
    v0 = norm(v1(1, :));
    v11 = norm(v1(2, :));
    v22 = norm(v1(3, :));
    if v0 < v11 && v0 < v22
    v1 = v1(1, :);
    v2 = v2(1, :);
    end
    if v11 < v0 && v11 < v22
    v1 = v1(2, :);
    v2 = v2(2, :);
    end
    if v22 < v0 && v22 < v11 
    v1 = v1(3, :);
    v2 = v2(3, :);
    end
end

% ---------- Key relative velocities ----------
v_inf_dep  = v1(:) - vE0(:);
v_inf_depN = norm(v_inf_dep);

v_inf_arr  = v2(:) - vAST(:);
v_inf_arrN = norm(v_inf_arr);

% ---------- Injection Δv from GTO perigee ----------
dvinj = sqrt(v_inf_depN^2 + 2*muE/rp) - sqrt(muE*(2/rp - 1/aGTO));

% ---------- Sample the Lambert arc for plotting ----------
N = 1600;
[Rarc, Varc] = sample_conic_by_time(rE0, v1, TOF, muS, N);

% Verification: end-point miss to asteroid
miss_km = norm(Rarc(:,end) - rAST(:));
fprintf('Endpoint miss to asteroid: %.3f km\n', miss_km);
fprintf('Injection Δv (GTO parking orbit):  %.3f km/s\n', dvinj);
fprintf('Arrival v_inf:             %.3f km/s\n', v_inf_arrN);

% ---------- Prepare 2D plot (heliocentric XY, actual scale) ----------
xT = Rarc(1,:); yT = Rarc(2,:);
xE = rE0(1);    yE = rE0(2);
xA = rAST(1);   yA = rAST(2);
xSun = 0;       ySun = 0;

figure('Color','w','units','normalized','outerposition',[0 0 1 1]); hold on;
% Sun
scatter(xSun,ySun,240,'filled','MarkerFaceColor',[1 0.84 0],'DisplayName','Sun');

% Trajectory
pT = plot(xT,yT,'k','LineWidth',2.4,'DisplayName','Ballistic Lambert arc');

% Key points
scatter(xE,yE,60,'b','filled'); 
text(xE,yE,'  Earth @ t_0','FontWeight','bold', 'VerticalAlignment','bottom');

scatter(xA,yA,60,[0.2 0.2 0.2],'filled'); 
text(xA,yA,'  Asteroid @ t_f','FontWeight','bold','VerticalAlignment','top');

% ---------- Arrow overlays (scaled to figure extent) ----------
rngX = max(xT) - min(xT);
rngY = max(yT) - min(yT);
arrow_scale = 0.05 * min(rngX, rngY);

% Departure ΔV_inj arrow direction (use v_inf_dep direction)
u_dep = (v_inf_dep / norm(v_inf_dep)) * arrow_scale;
quiver(xE, yE, u_dep(1), u_dep(2), 0, 'g','LineWidth',2,'MaxHeadSize',0.8);
text(xE + 1.2*u_dep(1), yE + 1.2*u_dep(2), '\DeltaV_{inj} (dep.)', ...
    'FontWeight','bold','HorizontalAlignment','left', 'VerticalAlignment','bottom');

% Arrival v_inf arrow
u_arr = (v_inf_arr / norm(v_inf_arr)) * arrow_scale;
quiver(xA, yA, u_arr(1), u_arr(2), 0, 'r','LineWidth',2,'MaxHeadSize',0.8);
text(xA + 1.2*u_arr(1), yA + 1.2*u_arr(2), 'v_{\infty,arr}', ...
    'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top');

axis equal; grid on; box on;
xlabel('X [km]'); ylabel('Y [km]');
title(sprintf('Ballistic Earth \\rightarrow 2024 YR4 (t_0=%s, t_f=%s)', ...
    datestr(t0_cal,'yyyy-mm-dd'), datestr(tf_cal,'yyyy-mm-dd')));
h = legend(pT, 'Location','best');
set(h,'FontSize',14);

% ========================== helpers =======================================
function [R, V] = sample_conic_by_time(r0, v0, TOF, mu, N)
    if nargin < 5, N = 1200; end
    t = linspace(0, TOF, N);
    R = zeros(3,N); V = zeros(3,N);
    state0 = [r0(:); v0(:)];
    for k = 1:N
        st = cspice_prop2b(mu, state0, t(k));
        R(:,k) = st(1:3);
        V(:,k) = st(4:6);
    end
end


