% AAE 537 Hypersonic Propulsion
% Project Sounding Rocket Sizing (Delta-V method w/ coast to hit gamma_t)
% Tyler Nord

clear; clc; close all;

%% Rugved look here
gamma_t_deg = 3;                               % [deg] FIXED target angle
M_target = 4;
H_target = 20000; %m

%% ---------------- constants/requirements ----------------
g0 = 9.80665;                                 % [m/s^2]
m_pl = 158.757;                                   % [kg] payload mass

H = linspace(10, 30, 100).*1000;              % target altitude grid [m]
M = linspace(2, 4, 100);                      % target Mach grid [-]

% Get speed of sound at each target altitude
[~, aH, ~, ~] = atmosisa(H, 'extended', true); % a(H) [m/s]

%% ---------------- propulsion assumptions ----------------
tb     = 30;                                   % [s] burn time (you can flex later)
Ivac   = 258.29;                               % [s] vacuum Isp
lambda = 0.8;                                  % [-] propellant mass fraction

%% ---------------- aero / drag assumptions ----------------
D      = 0.30;                                 % [m] diameter (set to your vehicle)
Aref   = pi*(D/2)^2;                           % [m^2] frontal area

% Burn drag model (average)
Cd_burn     = 0.35;                             % [-]
hEffFrac    = 0.30;                             % rho evaluated at hEffFrac*Hb for burn loss
k_drag_burn = 1.2;                              % [-] conservatism factor

% Coast drag model (ballistic coefficient average)
Cd_coast    = 0.25;                             % [-] typically smaller than powered AoA case
Hs          = 7200;                             % [m] scale height for rho-averaging
k_drag_coast= 1.3;                              % [-] conservatism

% Target flight-path angle at the target point
gamma_t = deg2rad(gamma_t_deg);

% Choose where burnout happens as a fraction of target altitude
hb_frac = 0.55;                                 % e.g., burn out around 55% of target altitude
hb_min  = 6000;                                 % [m]
hb_margin = 500;                                % [m] keep hb < Ht

% Iteration control (keeps m used in drag consistent with solved m_prop)
m_prop_seed = 400;                              % [kg] initial seed
nIter = 4;                                      % fixed-point iterations per grid point

%% ---------------- prealloc ----------------
v_t      = zeros(numel(H), numel(M));          % target speed at (Ht,M)
Vb_req   = zeros(numel(H), numel(M));          % required burnout speed
gammab   = zeros(numel(H), numel(M));          % required burnout gamma
t_coast  = zeros(numel(H), numel(M));          % coast time
dV_drag_b= zeros(numel(H), numel(M));          % burn drag loss
dV_drag_c= zeros(numel(H), numel(M));          % coast drag loss
dV_design= zeros(numel(H), numel(M));          % total ΔV used in rocket eq
m_prop   = zeros(numel(H), numel(M));          % solved propellant mass

%% ---------------- main sweep ----------------
for i = 1:length(M)
    for j = 1:length(H)

        Ht = H(j);
        at = aH(j);
        Vt = M(i) * at;                     % target speed at target point
        v_t(j,i) = Vt;

        % Choose burnout altitude for this target
        hb = hb_frac * Ht;
        hb = max(hb, hb_min);
        hb = min(hb, Ht - hb_margin);       % ensure hb < Ht

        % Fixed-point iteration for consistent masses in drag terms
        mp = m_prop_seed;

        for it = 1:nIter
            % Mass bookkeeping from lambda
            m_inert = mp*(1-lambda)/lambda;
            m0 = m_pl + m_inert + mp;                 % lift-off mass
            mb = m_pl + m_inert;                      % burnout (coast) mass

            % --- COAST (gravity-only kinematics) to hit (Ht, Vt, gamma_t) ---
            [Vb, gamma_b, tc] = coast_requirements(Ht, hb, Vt, gamma_t, g0);
            Vb_req(j,i) = Vb;
            gammab(j,i) = gamma_b;
            t_coast(j,i)= tc;

            % --- Burn drag loss (average, no stepping) ---
            % Use RMS speed during burn based on burnout speed requirement
            Vrms2_burn = Vb^2/3;

            % Evaluate density at an effective altitude during burn
            hEff_b = hEffFrac * hb;
            [~, ~, ~, rhoEff_b] = atmosisa(max(hEff_b,0),'extended',true);

            % Effective drag force (burn)
            Deff_b = 0.5 * rhoEff_b * Vrms2_burn * Cd_burn * Aref;

            % Use an effective mass during burn (average of m0 to mb)
            mEff_b = 0.5*(m0 + mb);

            dVb_drag = k_drag_burn * (Deff_b/mEff_b) * tb;

            % --- Coast drag loss (ballistic coefficient average) ---
            dh = Ht - hb;
            [~, ~, ~, rho_b] = atmosisa(max(hb,0),'extended',true);

            % altitude-avg density for exponential model from hb to Ht
            rho_eff_c = rho_b * (Hs/dh) * (1 - exp(-dh/Hs));

            % choose V^2 average during coast from Vb down to Vt
            Veff2_c = 0.5*(Vb^2 + Vt^2);

            beta = mb / (Cd_coast*Aref);   % ballistic coefficient at coast mass

            dVc_drag = k_drag_coast * (0.5*rho_eff_c*Veff2_c/beta) * tc;

            % --- Gravity loss during burn (conservative) ---
            dV_grav = g0 * tb;

            % Total ΔV for rocket equation (ΔV-level model)
            dVtot = Vb + dV_grav + dVb_drag + dVc_drag;

            % Rocket equation -> MR -> prop mass with lambda
            MR = exp(dVtot/(g0*Ivac));
            mp_new = m_pl*(MR-1)/(MR - (MR-1)/lambda);

            % Relax for stability
            mp = 0.6*mp + 0.4*mp_new;
        end

        % Store final iteration results
        m_prop(j,i)    = mp;
        dV_drag_b(j,i) = dVb_drag;
        dV_drag_c(j,i) = dVc_drag;
        dV_design(j,i) = dVtot;
    end
end

%% ---------------- plots ----------------
figure;
imagesc(M, H/1000, dV_design);
cb = colorbar; title(cb, "Design \DeltaV [m/s]");
axis xy; xlabel("Target Mach # at H_t"); ylabel("Target Altitude H_t [km]");
title(sprintf('Design \\DeltaV (coast to \\gamma_t = %g^\\circ)', gamma_t_deg));

figure;
imagesc(M, H/1000, m_prop);
cb = colorbar; title(cb, "Propellant Mass [kg]");
axis xy; xlabel("Target Mach # at H_t"); ylabel("Target Altitude H_t [km]");
title(sprintf('Propellant mass (coast to \\gamma_t = %g^\\circ)', gamma_t_deg));

figure;
imagesc(M, H/1000, rad2deg(gammab));
cb = colorbar; title(cb, "Required burnout \gamma_b [deg]");
axis xy; xlabel("Target Mach # at H_t"); ylabel("Target Altitude H_t [km]");
title('Required burnout flight-path angle to meet target');

figure;
imagesc(M, H/1000, t_coast);
cb = colorbar; title(cb, "Coast time t_c [s]");
axis xy; xlabel("Target Mach # at H_t"); ylabel("Target Altitude H_t [km]");
title('Coast time required');

%% ---------------- report one design point ----------------
iM = find(abs(M-M_target) == min(abs(M-M_target)),1);
jH = find(abs(H-H_target) == min(abs(H-H_target)),1);

fprintf('\n--- Example point ---\n');
fprintf('Target: H=%.0f m, M=%.2f, gamma_t=%.1f deg\n', H(jH), M(iM), gamma_t_deg);
fprintf('Propellant mass      : %.2f kg\n', m_prop(jH,iM));
fprintf('Inert mass           : %.2f kg\n', m_prop(jH,iM)*(1-lambda)/lambda);
fprintf('GLOW:                : %.2f kg\n', m_pl+m_prop(jH,iM)/lambda);
fprintf('Design ΔV            : %.1f m/s\n', dV_design(jH,iM));
fprintf('Burn drag ΔV         : %.1f m/s\n', dV_drag_b(jH,iM));
fprintf('Coast drag ΔV        : %.1f m/s\n', dV_drag_c(jH,iM));
fprintf('Required burnout Vb  : %.1f m/s\n', Vb_req(jH,iM));
fprintf('Required burnout γb  : %.2f deg\n', rad2deg(gammab(jH,iM)));
fprintf('Coast time tc        : %.2f s\n', t_coast(jH,iM));

%% --------- EXAMPLE POINT TRAJECTORY + q(t) (kinematic, no sim) ----------
Ht_ex   = H(jH);
Mt_ex   = M(iM);
Vt_ex   = v_t(jH,iM);
tc_ex   = t_coast(jH,iM);

% Burnout altitude consistent with policy
hb_ex = hb_frac*Ht_ex;
hb_ex = max(hb_ex, hb_min);
hb_ex = min(hb_ex, Ht_ex - hb_margin);

% Target velocity components at target point
gam_t_ex = deg2rad(gamma_t_deg);
vx_t = Vt_ex*cos(gam_t_ex);
vy_t = Vt_ex*sin(gam_t_ex);

% Burnout velocity components from coast kinematics
vx_b = vx_t;
vy_b = vy_t + g0*tc_ex;

fprintf('\n--- Kinematic trajectory build (example point) ---\n');
fprintf('hb used              : %.1f m\n', hb_ex);
fprintf('Vb (reconstructed)   : %.2f m/s\n', hypot(vx_b,vy_b));
fprintf('gamma_b (reconstruct): %.2f deg\n', rad2deg(atan2(vy_b,vx_b)));

% ---- Build time vector for burn+coast ----
dt_plot = 0.05;
t_all = (0:dt_plot:(tb+tc_ex)).';
Nall  = numel(t_all);

x_tr = zeros(Nall,1);
y_tr = zeros(Nall,1);
vx_tr = zeros(Nall,1);
vy_tr = zeros(Nall,1);
V_tr  = zeros(Nall,1);
q_tr  = zeros(Nall,1);
Mach_tr = zeros(Nall,1);

% ---- Burn vertical profile that guarantees continuity at burnout ----
% Enforce y(tb)=hb_ex AND vy(tb)=vy_b with vy(0)=0 using vy(t)=a t^2 + b t
A = [tb^2,      tb;
     tb^3/3,    tb^2/2];
ab = A \ [vy_b; hb_ex];
ay2 = ab(1);   % 'a' coefficient
by1 = ab(2);   % 'b' coefficient

% ---- NEW: Burn horizontal profile (NO linear term) so gamma starts ~90 deg ----
% Use vx(t)=c3*t^3 + c2*t^2 so vx ~ t^2 near liftoff (much smaller than vy ~ t)
% Enforce vx(tb)=vx_b and x(tb)=x_burn_des
x_burn_des = 0.5*vx_b*tb;  % keep prior convention for burnout downrange (just a kinematic choice)

Ax = [tb^3,      tb^2;
      tb^4/4,    tb^3/3];
c32 = Ax \ [vx_b; x_burn_des];
cx3 = c32(1);   % cubic coefficient
cx2 = c32(2);   % quadratic coefficient

for k = 1:Nall
    tnow = t_all(k);

    if tnow <= tb
        % ---- Burn segment ----
        % Horizontal: higher-order start -> gamma starts near 90 deg smoothly
        vx_tr(k) = cx3*tnow^3 + cx2*tnow^2;
        x_tr(k)  = (cx3/4)*tnow^4 + (cx2/3)*tnow^3;

        % Vertical: quadratic vy(t) and cubic y(t) (continuous at burnout)
        vy_tr(k) = ay2*tnow^2 + by1*tnow;
        y_tr(k)  = (ay2/3)*tnow^3 + (by1/2)*tnow^2;

    else
        % ---- Coast segment: ballistic from burnout (continuous start) ----
        tc = tnow - tb;

        x_burn = x_burn_des;
        y_burn = hb_ex;

        vx_tr(k) = vx_b;
        vy_tr(k) = vy_b - g0*tc;

        x_tr(k) = x_burn + vx_b*tc;
        y_tr(k) = y_burn + vy_b*tc - 0.5*g0*tc^2;
    end

    % Atmosphere at current altitude
    yk = max(y_tr(k),0);
    [~, a_k, ~, rho_k] = atmosisa(yk,'extended',true);

    V_tr(k) = hypot(vx_tr(k), vy_tr(k));
    q_tr(k) = 0.5*rho_k*V_tr(k)^2;
    Mach_tr(k) = V_tr(k)/max(a_k,1e-9);
end



% Plot trajectory
figure; plot(x_tr/1000, y_tr/1000,'LineWidth',1.8); grid on;
xlabel('Downrange [km]'); ylabel('Altitude [km]');
title(sprintf('Example trajectory (kinematic, continuous at burnout) — target: H=%.1f km, M=%.1f, \\gamma_t=%g^\\circ', ...
    Ht_ex/1000, Mt_ex, gamma_t_deg));

% Plot dynamic pressure vs time
figure; plot(t_all, q_tr/1000,'LineWidth',1.8); grid on;
xlabel('Time [s]'); ylabel('q [kPa]');
title('Dynamic pressure vs time (example point)');

gamma_tr = atan2(vy_tr, vx_tr);   % [rad]

figure;
subplot(2,2,1)
plot(t_all, y_tr/1000,'LineWidth',1.6); grid on;
ylabel('Altitude [km]')

subplot(2,2,2)
plot(t_all, x_tr/1000,'LineWidth',1.6); grid on;
ylabel('Downrange [km]')

subplot(2,2,3)
plot(t_all, Mach_tr,'LineWidth',1.6); grid on;
xlabel('Time [s]'); ylabel('Mach')

subplot(2,2,4)
plot(t_all, rad2deg(gamma_tr),'LineWidth',1.6); grid on;
xlabel('Time [s]'); ylabel('\gamma [deg]')

%% --------- SAVE EXAMPLE TRAJECTORY DATA ----------
t     = t_all;            % [s] time
x     = x_tr;             % [m] downrange
h     = y_tr;             % [m] altitude
M     = Mach_tr;          % [-] Mach number
angle = rad2deg(gamma_tr);   % [deg] flight-path angle
q     = q_tr;             % [Pa] dynamic pressure

save('example_trajectory.mat', 't', 'x', 'h', 'M', 'angle', 'q');

%% ===================== local function =====================
function [Vb, gamma_b, tc] = coast_requirements(Ht, hb, Vt, gamma_t, g)
% Gravity-only coast from burnout (hb, Vb, gamma_b) to target (Ht, Vt, gamma_t).
% Assumes vx constant during coast, vy decreases linearly with time.

    dh = Ht - hb;
    if dh <= 0
        error('hb must be below Ht');
    end

    vx_t = Vt*cos(gamma_t);
    vy_t = Vt*sin(gamma_t);

    % dh = vy_t*tc + 0.5*g*tc^2  (derived by eliminating vy_b)
    tc = (-vy_t + sqrt(vy_t^2 + 2*g*dh))/g;

    vy_b = vy_t + g*tc;
    vx_b = vx_t;

    Vb = hypot(vx_b, vy_b);
    gamma_b = atan2(vy_b, vx_b);
end

