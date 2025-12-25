clear;clc;
startup;

% --- SPICE setup (adjust paths for your environment) ---------------------
cspice_kclear;
load_kernels;

% ---------- plot the optimized trajectory MGADSM Earth ----------
t0 = datetime(2028, 10, 14);
tf = datetime(2032, 11, 7);
t0 = cspice_str2et(datestr(t0));
tf = cspice_str2et(datestr(tf));
TOF = tf - t0;
eta1 = 0.712;
eta2 = 0.773;
theta_b = -1.75279421621431;
Rperi = 36463.1375370251;

% -------------------- Constants & GTO orbit info ---------------------------
muS = 1.32712440018e11;   % Sun GM [km^3/s^2]
muE = 398600;             % Earth GM [km^3/s^2]
Re  = 6378.137;           % Earth radius [km]
rp  = Re + 250;           % GTO perigee radius [km]
ra  = Re + 22500;         % GTO apogee radius [km]
a   = 0.5*(ra + rp);      % GTO semi-major axis [km]

% -------------------- Event times from decision variables ------------------
T1       = eta1 * TOF;
tGA      = t0 + T1;
tArr     = t0 + TOF;
T2       = (1 - eta1) * TOF;
tGA2DSM  = eta2 * T2;
tDSM2Ast = (1 - eta2) * T2;

% -------------------- Planetary states (Sun/J2000) -------------------------
D_earth1 = cspice_spkezr('Earth', t0,  'J2000','NONE','Sun');
D_earth2 = cspice_spkezr('Earth', tGA, 'J2000','NONE','Sun');
D_ast = cspice_spkezr('54509621', tArr, 'J2000','NONE','Sun');    % 2024 YR4

rE1 = D_earth1(1:3); vE1 = D_earth1(4:6);
rE2 = D_earth2(1:3); vE2 = D_earth2(4:6);
rAst = D_ast(1:3); vAst = D_ast(4:6);

% -------------------- Lambert arc: Earth(t0) → Earth(tGA) ------------------
[v1, v2] = lambert_solver(rE1, rE2, T1, muS);   % may return multiple branches

if size(v1,1) > 1
    v0 = norm(v1(1, :));
    v11 = norm(v1(2, :));
    v22 = norm(v1(3, :));
    if v0 < v11 && v0 < v22
    v1 = v1(1, :);
    v2 = v2(1, :);

    elseif v11 < v0 && v11 < v22
    v1 = v1(2, :);
    v2 = v2(2, :);
    else 
    v1 = v1(3, :);
    v2 = v2(3, :);
    end
end

% -------------------- Departure injection Δv from GTO ----------------------
v_inf1    = v1(:) - vE1(:);
v_inf1_n  = norm(v_inf1);
delta_v_inj = sqrt(v_inf1_n^2 + 2*muE/rp) - sqrt(muE*(2/rp - 1/a));

% -------------------- Unpowered GA at Earth --------------------------------
v_inf2 = v2(:) - vE2(:);
v_out  = GA_Earth_calculations(Rperi, theta_b, v_inf2, vE2, tGA);  % heliocentric

% -------------------- Propagate GA → DSM (Sun two-body) --------------------
dt_GA2DSM = tGA2DSM;                                  % seconds from tGA to DSM
stateGA   = [rE2(:); v_out(:)];                       % 6x1 at GA epoch
stateDSM  = cspice_prop2b(muS, stateGA, dt_GA2DSM);   % propagate by Δt

r_pre_DSM = stateDSM(1:3);
v_pre_DSM = stateDSM(4:6);

% -------------------- DSM → asteroid Lambert arc ---------------------------
[v_post_DSM, v_sc_ast] = lambert_solver(r_pre_DSM(:), rAst(:), tDSM2Ast, muS);
if size(v_post_DSM) > 1
    v0 = norm(v_post_DSM(1, :));
    v11 = norm(v_post_DSM(2, :));
    v12 = norm(v_post_DSM(3, :));
    if v0 < v11 && v0 < v12
    v_post_DSM = v_post_DSM(1, :);
    v_sc_ast = v_sc_ast(1, :);
    end
    if v11 < v0 && v11 < v12
    v_post_DSM = v_post_DSM(2, :);
    v_sc_ast = v_sc_ast(2, :);
    end
    if v12 < v0 && v12 < v11
    v_post_DSM = v_post_DSM(3, :);
    v_sc_ast = v_sc_ast(3, :);
    end
end

% -------------------- DSM impulse and arrival v∞ ---------------------------
delta_v_DSM   = v_post_DSM(:) - v_pre_DSM(:);
delta_v_DSM_n = norm(delta_v_DSM);
v_inf_arr     = v_sc_ast(:) - vAst(:);
v_inf_arr_n   = norm(v_inf_arr);

% -- sampling densities --
N1 = 900;           % Earth->Earth (resonant) leg
N2 = 1200;          % GA->DSM 
N3 = 900;           % DSM->Arrival

% --- Leg 1: Earth -> Earth using Lambert solution ---
[R1s, V1s] = sample_conic_by_time(rE1, v1, T1, muS, N1);

% Verification: endpoint should match Earth position at tGA
st_end = R1s(:,end);
fprintf('Leg1 miss to Earth at tGA: %.1f km\n', norm(st_end - rE2(:)));

% --- Leg 2: Earth GA -> DSM using post-GA velocity ---
[R2s, V2s] = sample_conic_by_time(rE2, v_out, tGA2DSM, muS, N2);
R_DSM_computed = R2s(:,end);
v_pre_DSM_computed = V2s(:,end);

% Use the consistently computed DSM position and velocity
R_DSM = R_DSM_computed;
v_pre_DSM = v_pre_DSM_computed;

% --- Leg 3: DSM -> Asteroid using Lambert solution ---
[R3s, V3s] = sample_conic_by_time(R_DSM, v_post_DSM, tDSM2Ast, muS, N3);

% Verification: endpoint should match asteroid position
st_end3 = R3s(:,end);
fprintf('Leg3 miss to Asteroid: %.1f km\n', norm(st_end3 - rAst(:)));

% --- key points for labels (no coordinate transformation) ---
RE_t0  = rE1(:);  
RE_tGA = rE2(:);  
R_AST = rAst(:);

% -------- Remove distortion - use actual coordinates --------
x1 = R1s(1,:); y1 = R1s(2,:);
x2 = R2s(1,:); y2 = R2s(2,:);  
x3 = R3s(1,:); y3 = R3s(2,:);

xe0 = RE_t0(1);   ye0 = RE_t0(2);
xeGA = RE_tGA(1); yeGA = RE_tGA(2);
xDSM = R_DSM(1);  yDSM = R_DSM(2);
xAST = R_AST(1);  yAST = R_AST(2);
xSun = 0;         ySun = 0;

% -------- plot --------
figure('Color','w','units','normalized','outerposition',[0 0 1 1]); hold on;

% Plot Sun
scatter(xSun,ySun,240,'filled','MarkerFaceColor',[1 0.84 0],'DisplayName','Sun');

% Plot trajectory legs with proper continuity
p1 = plot(x1,y1,'r','LineWidth',2.4, ...
    'DisplayName','Leg 1: Earth \rightarrow Earth (resonant)');
p2 = plot(x2,y2,'b','LineWidth',2.0, ...
    'DisplayName','Leg 2 segment 1: GA \rightarrow DSM');
p3 = plot(x3,y3,'k','LineWidth',2.4, ...
    'DisplayName','Leg 2 segment 2: DSM \rightarrow Asteroid');

% Plot key points
scatter(xe0, ye0, 60,'b','filled');  
text(xe0, ye0,'  Earth @ t_0','FontWeight','bold', 'VerticalAlignment','bottom');

scatter(xeGA, yeGA, 60, 'r', 'filled');
text(xeGA, yeGA, '  Earth @ t_{GA}', 'FontWeight','bold', ...
     'VerticalAlignment','top','HorizontalAlignment','left');

scatter(xDSM, yDSM, 55,[0.6 0 0.8],'filled'); 
text(xDSM, yDSM,'  DSM','FontWeight','bold','Color',[0.4 0 0.6], ...
     'HorizontalAlignment','right', 'VerticalAlignment','bottom');

scatter(xAST, yAST, 55,[0.2 0.2 0.2],'filled'); 
text(xAST, yAST,'2024 YR4','FontWeight','bold', 'VerticalAlignment','bottom', ...
     'HorizontalAlignment', 'right');

% Add velocity arrows at key points with proper scaling
% Calculate range for arrow scaling
rngX = max([x1,x2,x3]) - min([x1,x2,x3]);
rngY = max([y1,y2,y3]) - min([y1,y2,y3]);
arrow_scale = min(rngX, rngY) * 0.05; % 5% of plot range

% Departure velocity arrow (normalised and scaled)
v_inf1_norm = v_inf1 / norm(v_inf1) * arrow_scale;
quiver(xe0, ye0, v_inf1_norm(1), v_inf1_norm(2), 0, 'g','LineWidth',2,'MaxHeadSize',0.8);
text(xe0 + v_inf1_norm(1)*1.2, ye0 + v_inf1_norm(2)*1.2, '\DeltaV_1 (departure)', ...
     'FontWeight','bold', 'HorizontalAlignment','left');

% DSM velocity arrow (normalised and scaled)
delta_v_DSM_norm = delta_v_DSM / norm(delta_v_DSM) * arrow_scale;
quiver(xDSM, yDSM, delta_v_DSM_norm(1), delta_v_DSM_norm(2), 0, 'm','LineWidth',2,'MaxHeadSize',0.8);
text(xDSM + delta_v_DSM_norm(1)*1.2, yDSM + delta_v_DSM_norm(2)*1.2, '\DeltaV_{DSM}', ...
     'FontWeight','bold', 'HorizontalAlignment','left','VerticalAlignment','bottom');

% Arrival velocity arrow (normalised and scaled)
v_inf_arr_norm = v_inf_arr / norm(v_inf_arr) * arrow_scale;
quiver(xAST, yAST, v_inf_arr_norm(1), v_inf_arr_norm(2), 0, 'r','LineWidth',2,'MaxHeadSize',0.8);
text(xAST + v_inf_arr_norm(1)*2, yAST + v_inf_arr_norm(2)*2, 'v_{\infty,arr}', ...
     'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','bottom');

axis equal; grid on; box on;
xlabel('X [km]'); ylabel('Y [km]');
title('MGADSM – Resonant Earth Flyby Trajectory');
h = legend([p1 p2 p3],'Location','best');
set(h,'FontSize',14);
print(gcf,'-dpdf','MGADSM_resonant_trajectory_fixed.pdf');

disp('DeltaV_inj');   disp(delta_v_inj)
disp('DeltaV_DSM');   disp(delta_v_DSM_n)
disp('DeltaV_arr');   disp(v_inf_arr_n)

% ========================== helper function ================================
% --- helper: sample conic by time (no branch logic, no clipping) ---
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


