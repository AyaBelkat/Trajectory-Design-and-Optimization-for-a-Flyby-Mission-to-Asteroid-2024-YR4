clear;clc;
% --- SPICE setup (adjust paths for your environment) ---------------------
addpath("C:\Users\wbook\Desktop\Dissertation files\Dissertation_codes\mice\mice\src\mice");
addpath("C:\Users\wbook\Desktop\Dissertation files\Dissertation_codes\mice\mice\lib");

cspice_kclear;
cspice_furnsh('C:\Users\wbook\Desktop\Dissertation files\Dissertation_codes\mission_meta.tm');
cspice_furnsh('C:\Users\wbook\Desktop\Dissertation files\Dissertation_codes\mission_equa_data.tm');
% ---------- plot the optimized trajectory MGADSM Venus ----------
t0 = datetime(2031, 5, 27);
tf = datetime(2032, 10, 13);
t0 = cspice_str2et(datestr(t0));
tf = cspice_str2et(datestr(tf));
TOF = tf - t0;
eta1 = 0.289;
eta2 = 0.130;
theta_b = -1.29904089014891;
Rperi = 17420.1897023278;
% -------------------- Constants & GTO orbit info ----------
muS = 1.32712440018e11;   % Sun GM [km^3/s^2]
muE = 398600;             % Earth GM [km^3/s^2]
Re  = 6378.137;           % Earth radius [km]
rp  = Re + 250;           % GTO perigee radius [km]
ra  = Re + 22500;         % GTO apogee radius [km]
a   = 0.5*(ra + rp);      % GTO semi-major axis [km]

% -------------------- Event times from decision variables ------------------
T1   = eta1 * TOF;                  % Earth→GA leg duration
tGA  = t0 + T1;                     % GA epoch
tArr = t0 + TOF;                    % Arrival epoch
T2   = (1 - eta1) * TOF;            % GA→arrival total
tGA2DSM   = eta2 * T2;              % GA→DSM duration
tDSM2Ast  = (1 - eta2) * T2;        % DSM→arrival duration

% -------------------- Planetary states (Sun/J2000) ------------------------
D_earth1 = cspice_spkezr('Earth', t0,  'J2000','NONE','Sun');   
D_venus2 = cspice_spkezr('Venus', tGA, 'J2000','NONE','Sun');
D_ast = cspice_spkezr('54509621', tArr, 'J2000','NONE','Sun');   % 2024 YR4

rE1 = D_earth1(1:3); vE1 = D_earth1(4:6);
rV2 = D_venus2(1:3); vV2 = D_venus2(4:6);
rAst = D_ast(1:3); vAst = D_ast(4:6);

% -------------------- Lambert arc: Earth(t0) → Venus(tGA) -----------------
% Heliocentric velocities at endpoints (may return multiple branches)
[v1, v2] = lambert_solver(rE1, rV2, T1, muS);
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

% -------------------- Departure injection Δv from GTO ---------------------
% Required Earth-relative v∞ at departure:
v_inf1    = v1(:) - vE1;                    % [km/s]
v_inf1_n  = norm(v_inf1);

% Tangential perigee burn to achieve that v∞ from GTO:
% Δv0 = v_peri(hyper) − v_peri(ellipse) 
delta_v_inj = sqrt(v_inf1_n^2 + 2*muE/rp) - sqrt(muE*(2/rp - 1/a));

% -------------------- Unpowered GA at Venus -------------------------------
% Inbound v∞ at GA epoch (heliocentric → planet-relative):
v_inf2 = v2(:) - vV2;                       % inbound v∞ at GA (not charged)
v_out  = GA_Venus_calculations(Rperi, theta_b, v_inf2, vV2, tGA);    % GA calculations

% -------------------- Propagate GA → DSM (Sun two-body) -----------------------
dt_GA2DSM = tGA2DSM;                                  % seconds from tGA to DSM
stateGA   = [rV2(:); v_out(:)];                       % 6x1 at GA epoch
stateDSM  = cspice_prop2b(muS, stateGA, dt_GA2DSM);   % propagate by Δt

r_pre_DSM = stateDSM(1:3);
v_pre_DSM = stateDSM(4:6);

% -------------------- DSM → asteroid Lambert arc --------------------------
[v_post_DSM, v_sc_ast] = lambert_solver(r_pre_DSM, rAst, tDSM2Ast, muS);

% Calculate the velocity for each revolution and select the one with least
% energy (keeping your Lambert selection logic unchanged)
if size(v_post_DSM,1) > 1
    v0 = norm(v_post_DSM(1,:));
    v11 = norm(v_post_DSM(2,:));
    v22 = norm(v_post_DSM(3,:));
    if v0 < v11 && v0 < v22
    v_post_DSM = v_post_DSM(1, :);
    v_sc_ast = v_sc_ast(1, :);

    elseif v11 < v0 && v11 < v22
    v_post_DSM = v_post_DSM(2, :);
    v_sc_ast = v_sc_ast(2, :);

    else
    v_post_DSM = v_post_DSM(3, :);
    v_sc_ast = v_sc_ast(3, :);
    end
end

% -------------------- DSM impulse and arrival v∞ -------------
delta_v_DSM   = v_post_DSM(:) - v_pre_DSM(:);   % DSM vector Δv [km/s]
delta_v_DSM_n = norm(delta_v_DSM);

% Arrival hyperbolic excess (not charged in J, but useful for reporting)
v_inf_arr = v_sc_ast(:) - vAst;
v_inf_arr_n = norm(v_inf_arr);

% -- sampling densities --
N1 = 1200;           % Earth->Venus leg
N2 = 1200;          % GA->DSM 
N3 = 900;           % DSM->Arrival

% --- Leg 1: Earth -> Venus using Lambert solution ---
[R1s, V1s] = sample_conic_by_time(rE1, v1, T1, muS, N1);

% Verification: endpoint should match Venus position
st_end = R1s(:,end);
fprintf('Leg1 miss to Venus: %.1f km\n', norm(st_end - rV2(:)));

% --- Leg 2: Venus GA -> DSM using post-GA velocity ---
[R2s, V2s] = sample_conic_by_time(rV2, v_out, tGA2DSM, muS, N2);
R_DSM_computed = R2s(:,end);
v_pre_DSM_computed = V2s(:,end);

% Use the consistently computed DSM position
R_DSM = R_DSM_computed;
v_pre_DSM = v_pre_DSM_computed;

% --- Leg 3: DSM -> Asteroid using Lambert solution ---
[R3s, V3s] = sample_conic_by_time(R_DSM, v_post_DSM, tDSM2Ast, muS, N3);

% Verification: endpoint should match asteroid position
st_end3 = R3s(:,end);
fprintf('Leg3 miss to 2024 YR4: %.1f km\n', norm(st_end3 - rAst(:)));

% --- key points for labels (no coordinate transformation) ---
RE_t0  = rE1(:); 

RE_tGA = rV2(:);  
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
    'DisplayName','Leg 1: Earth \rightarrow Venus');
p2 = plot(x2,y2,'b','LineWidth',2.0, ...
    'DisplayName','Leg 2 segment 1: Venus \rightarrow DSM');
p3 = plot(x3,y3,'k','LineWidth',2.4, ...
    'DisplayName','Leg 2 segment 2: DSM \rightarrow 2024 YR4');

% Plot key points
scatter(xe0, ye0, 60,'b','filled');  
text(xe0, ye0,'  Earth @ t_0','FontWeight','bold', 'VerticalAlignment','bottom');

scatter(xeGA, yeGA, 60, 'r', 'filled');
text(xeGA, yeGA, '  Venus @ t_{GA}', 'FontWeight','bold', ...
     'VerticalAlignment','top','HorizontalAlignment','right');

scatter(xDSM, yDSM, 55,[0.6 0 0.8],'filled'); 
text(xDSM, yDSM,'  DSM','FontWeight','bold','Color',[0.4 0 0.6], 'VerticalAlignment','bottom');

scatter(xAST, yAST, 55,[0.2 0.2 0.2],'filled'); 
text(xAST, yAST,'  2024 YR4','FontWeight','bold', 'VerticalAlignment','top');

% Add velocity arrows at key points (simplified, no complex tangent calculation)
% Departure velocity arrow
quiver(xe0, ye0, -v_inf1(1)*1e7, v_inf1(2)*1e7, 0, 'g','LineWidth',2,'MaxHeadSize',0.5);
text(xe0, ye0-2e7, '\DeltaV_1 (departure)', 'FontWeight','bold', 'HorizontalAlignment','center');

% DSM velocity arrow  
quiver(xDSM, yDSM, delta_v_DSM(1)*2e7, delta_v_DSM(2)*2e7, 0, 'm','LineWidth',2,'MaxHeadSize',0.5);
text(xDSM, yDSM-1.5e7, '\DeltaV_{DSM}', 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','bottom');

% Arrival velocity arrow 
quiver(xAST, yAST, -v_inf_arr(1)*1e7,- v_inf_arr(2)*1e7, 0, 'm','LineWidth',2,'MaxHeadSize',0.5);
text(xAST, yAST-2e7, '\DeltaV_(arrival)', 'FontWeight','bold', 'HorizontalAlignment','right', 'VerticalAlignment','bottom');


axis equal; grid on; box on;
xlabel('X [km]'); ylabel('Y [km]');
title('MGADSM – Venus Flyby Trajectory');
h = legend([p1 p2 p3],'Location','best');
set(h,'FontSize',14);
print(gcf,'-dpdf','MGADSM_trajectory_fixed.pdf');

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
