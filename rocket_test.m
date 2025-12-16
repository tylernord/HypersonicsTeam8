%% rocket_prescribed_pitch_sweep_fpitch.m
% Sweep f_pitch values for a prescribed BODY PITCH program.
% Plots trajectory, Mach(t), q(t), gamma(t), AoA(t) for each case.
% Uses a CELL array for results to avoid "dissimilar structures" errors.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% THIS CODE IS OUTDATED: DO NOT USE
% USE rocket_test_2.m INSTEAD



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% ---------------- USER INPUTS ----------------
g0 = 9.80665;                 % gravity [m/s^2]

% Vehicle / propulsion
m_pl   = 160;                 % payload mass [kg]
m_prop = 900;                 % propellant mass [kg] (chosen)
lambda = 0.80;                % propellant mass fraction
Isp    = 258.29;              % vacuum Isp [s]
tb     = 26.53;               % burn time [s]

% Geometry / aero
D      = 0.30;                                % diameter [m]
Aref   = pi*(D/2)^2;                          % frontal area [m^2]
Cd0    = 0.25;                                % zero-AoA Cd
Cd_max = 1.0;                                 % cap
kCd    = 10.0;                                % Cd(alpha) strength

% AoA → projected area amplification
A_side_factor = 8.0;                          % side-on penalty strength
alpha_max = deg2rad(60);                      % AoA limit

% Integration
dt = 0.01;                                    % time step [s]

% Initial conditions
x0 = 0; y0 = 0;
V0 = 1;                                       % [m/s]
gamma0 = deg2rad(90);                         % vertical

% -------- Pitch program knobs (fixed across sweep) --------
theta0_deg = 90;      % initial body pitch
thetaf_deg = 10;      % final body pitch
p_shape    = 1.5;     % <1 = aggressive early pitch, >1 = later pitch

% >>> Sweep these <<<
f_pitch_vec = [0.2 0.4 0.6 0.8 0.9 0.95];     % fraction of burn before pitching begins

%% ---------------- LOOP OVER f_pitch ----------------
results = cell(numel(f_pitch_vec),1);

fprintf('\n===== SWEEP RESULTS =====\n');

for ii = 1:numel(f_pitch_vec)
    f_pitch = f_pitch_vec(ii);

    out = run_case( ...
        g0, m_pl, m_prop, lambda, Isp, tb, ...
        Aref, Cd0, Cd_max, kCd, A_side_factor, alpha_max, ...
        dt, x0, y0, V0, gamma0, ...
        theta0_deg, thetaf_deg, f_pitch, p_shape);

    results{ii} = out;

    fprintf('f_pitch = %.3f | y(tb)=%.0f m | V(tb)=%.0f m/s | Mach(tb)=%.3f | gamma(tb)=%.2f deg | max q=%.0f kPa\n', ...
        f_pitch, out.y(end), out.Vend, out.Mach(end), rad2deg(out.gamma(end)), max(out.qbar)/1000);
end

%% ---------------- PLOTS ----------------

% 1) Trajectory y vs x
figure; hold on; grid on;
for ii = 1:numel(results)
    plot(results{ii}.x/1000, results{ii}.y/1000, 'LineWidth', 1.4, ...
        'DisplayName', sprintf('f=%.2f', results{ii}.f_pitch));
end
xlabel('Downrange [km]'); ylabel('Altitude [km]');
title('Trajectory (prescribed pitch) — sweep f\_pitch');
legend('Location','best');

% 2) Mach vs time
figure; hold on; grid on;
for ii = 1:numel(results)
    plot(results{ii}.t, results{ii}.Mach, 'LineWidth', 1.4, ...
        'DisplayName', sprintf('f=%.2f', results{ii}.f_pitch));
end
xlabel('Time [s]'); ylabel('Mach');
title('Mach vs time — sweep f\_pitch');
legend('Location','best');

% 3) qbar vs time
figure; hold on; grid on;
for ii = 1:numel(results)
    plot(results{ii}.t, results{ii}.qbar/1000, 'LineWidth', 1.4, ...
        'DisplayName', sprintf('f=%.2f', results{ii}.f_pitch));
end
xlabel('Time [s]'); ylabel('q [kPa]');
title('Dynamic pressure vs time — sweep f\_pitch');
legend('Location','best');

% 4) Velocity flight-path angle gamma vs time
figure; hold on; grid on;
for ii = 1:numel(results)
    plot(results{ii}.t, rad2deg(results{ii}.gamma), 'LineWidth', 1.4, ...
        'DisplayName', sprintf('f=%.2f', results{ii}.f_pitch));
end
xlabel('Time [s]'); ylabel('\gamma_v [deg]');
title('Velocity flight-path angle vs time — sweep f\_pitch');
legend('Location','best');

% 5) Angle of attack alpha vs time
figure; hold on; grid on;
for ii = 1:numel(results)
    plot(results{ii}.t, rad2deg(results{ii}.alpha), 'LineWidth', 1.4, ...
        'DisplayName', sprintf('f=%.2f', results{ii}.f_pitch));
end
xlabel('Time [s]'); ylabel('\alpha [deg]');
title('Angle of attack vs time — sweep f\_pitch');
legend('Location','best');

% 6) Summary plot: burnout altitude vs f_pitch
fp = zeros(numel(results),1);
yb = zeros(numel(results),1);
Mb = zeros(numel(results),1);
qb = zeros(numel(results),1);
gb = zeros(numel(results),1);

for ii = 1:numel(results)
    fp(ii) = results{ii}.f_pitch;
    yb(ii) = results{ii}.y(end);
    Mb(ii) = results{ii}.Mach(end);
    qb(ii) = max(results{ii}.qbar);
    gb(ii) = rad2deg(results{ii}.gamma(end));
end

figure; grid on; hold on;
plot(fp, yb/1000, '-o', 'LineWidth', 1.4);
xlabel('f\_pitch [-]'); ylabel('Burnout altitude [km]');
title('Burnout altitude vs f\_pitch');

figure; grid on; hold on;
plot(fp, Mb, '-o', 'LineWidth', 1.4);
xlabel('f\_pitch [-]'); ylabel('Burnout Mach [-]');
title('Burnout Mach vs f\_pitch');

figure; grid on; hold on;
plot(fp, qb/1000, '-o', 'LineWidth', 1.4);
xlabel('f\_pitch [-]'); ylabel('Max q [kPa]');
title('Max dynamic pressure vs f\_pitch');

figure; grid on; hold on;
plot(fp, gb, '-o', 'LineWidth', 1.4);
xlabel('f\_pitch [-]'); ylabel('Burnout flight angle \gamma_v [deg]');
title('Burnout flight angle vs f\_pitch');

%% ======================= LOCAL FUNCTION =======================
function out = run_case( ...
    g0, m_pl, m_prop, lambda, Isp, tb, ...
    Aref, Cd0, Cd_max, kCd, A_side_factor, alpha_max, ...
    dt, x0, y0, V0, gamma0, ...
    theta0_deg, thetaf_deg, f_pitch, p_shape)

    t = (0:dt:tb).';  N = numel(t);

    % derived masses and thrust
    m_inert = m_prop*(1-lambda)/lambda;
    m0 = m_pl + m_inert + m_prop;

    mdot = m_prop / tb;
    T = Isp * g0 * mdot;

    theta0 = deg2rad(theta0_deg);
    thetaf = deg2rad(thetaf_deg);

    % prealloc
    x  = zeros(N,1);  y  = zeros(N,1);
    vx = zeros(N,1);  vy = zeros(N,1);
    m  = zeros(N,1);

    theta = zeros(N,1);
    gamma = zeros(N,1);
    alpha = zeros(N,1);
    Cd    = zeros(N,1);
    Ddrag = zeros(N,1);
    Mach  = zeros(N,1);
    qbar  = zeros(N,1);

    % init
    x(1)=x0; y(1)=y0;
    vx(1)=V0*cos(gamma0);
    vy(1)=V0*sin(gamma0);

    % pitch schedule
    t_pitch = f_pitch * tb;
    for k = 1:N
        if t(k) <= t_pitch
            theta(k) = theta0;
        else
            tau = (t(k)-t_pitch)/max(tb-t_pitch,eps);
            tau = min(max(tau,0),1);
            tau = tau^p_shape;
            theta(k) = theta0 + (thetaf-theta0)*tau;
        end
    end
    theta(end) = thetaf;

    % integrate
    for k = 1:N-1
        tk = t(k);

        m(k) = m0 - mdot*tk;

        yk = max(y(k),0);
        [~, a, ~, rho] = atmosisa(yk,'extended',true);

        V = hypot(vx(k), vy(k));
        gamma(k) = atan2(vy(k), vx(k));

        % vhat always defined
        if V < 1e-6
            vhat = [cos(theta(k)); sin(theta(k))];
        else
            vhat = [vx(k); vy(k)]/V;
        end

        alpha(k) = wrapToPi(theta(k) - gamma(k));
        alpha_eff = max(min(alpha(k), alpha_max), -alpha_max);

        Cd(k) = min(Cd0 + kCd*alpha_eff^2, Cd_max);
        Aproj = Aref * (1 + A_side_factor*abs(sin(alpha_eff)));

        Ddrag(k) = 0.5*rho*V^2 * Cd(k) * Aproj;

        qbar(k) = 0.5*rho*V^2;
        Mach(k) = V/max(a,1e-9);

        Tx = T*cos(theta(k));
        Ty = T*sin(theta(k));

        Fx = Tx - Ddrag(k)*vhat(1);
        Fy = Ty - Ddrag(k)*vhat(2) - m(k)*g0;

        ax = Fx/m(k);
        ay = Fy/m(k);

        vx(k+1) = vx(k) + ax*dt;
        vy(k+1) = vy(k) + ay*dt;

        x(k+1)  = x(k) + vx(k+1)*dt;
        y(k+1)  = y(k) + vy(k+1)*dt;

        if y(k+1) < 0
            y(k+1) = 0;
            vy(k+1) = max(vy(k+1),0);
        end
    end

    % final bookkeeping
    m(end) = m0 - mdot*tb;
    [~, a_end, ~, rho_end] = atmosisa(max(y(end),0),'extended',true);
    Vend = hypot(vx(end),vy(end));
    gamma(end) = atan2(vy(end),vx(end));
    alpha(end) = wrapToPi(theta(end)-gamma(end));
    Cd(end) = min(Cd0 + kCd*alpha(end)^2, Cd_max);
    qbar(end) = 0.5*rho_end*Vend^2;
    Mach(end) = Vend/max(a_end,1e-9);
    Ddrag(end) = 0.5*rho_end*Vend^2 * Cd(end) * Aref;

    % pack outputs
    out = struct();
    out.f_pitch = f_pitch;
    out.t = t;
    out.x = x; out.y = y;
    out.vx = vx; out.vy = vy;
    out.theta = theta;
    out.gamma = gamma;
    out.alpha = alpha;
    out.Cd = Cd;
    out.Ddrag = Ddrag;
    out.Mach = Mach;
    out.qbar = qbar;
    out.Vend = Vend;
end
