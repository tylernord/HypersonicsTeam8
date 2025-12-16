% AAE 537 Hypersonic Propulsion
% Project Sounding Rocket Sizing 
% Tyler Nord
clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% THIS CODE IS OUTDATED: DO NOT USE
% USE rocket_test_2.m INSTEAD



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%constants/requirements
g = 9.81;                                                   %gravitational acceleration, m/s^2
m_pl = 160;                                                 %payload mass, kg
H = linspace(10, 30, 100).*1000;                            %take over altitude, m
[~, a, ~, ~, ~, ~] = atmosisa(H, 'extended', true);         %take over speed of sound, m/s
M = linspace(2, 4, 100);                                    %take over Mach number

%assumptions taken from Black Brant V
tb = 30;                                                 %burn time, s
Ivac = 258.29;                                              %specific impulse, sec
lambda = 0.8;                                               %propellant mass fraction
Ae = 0.0773;                                                %nozzle exit area, m^2

% drag stuff
Cd   = 0.35;           % [-] assumed average Cd over burn
Aref = 0.0707;         % [m^2] ref area (D=0.30 m => 0.0707)
hEffFrac = 0.30;       % [-] effective altitude fraction for rho_eff (0.2–0.4 typical)
k_drag = 1.2;          % [-] fudge factor (>=1 conservative)
m_prop_guess = 240;    % kg
m0 = m_pl + (m_prop_guess)/lambda; % needs m_prop -> see note below
m1 = m0 - m_prop_guess;
mEff = 0.5*(m0+m1);

for i = 1:length(M)
    for j = 1:length(a)
        v(j,i) = M(i)*a(j);                                 %take over velocity, m/s
        % dV_ideal(j,i) = sqrt(v(j,i)^2+2*g*H(j));            %from conservation of energy, neglects drag, steering, and grav losses
        Vrms2 = v(j,i)^2/3;                % <V^2> assuming linear ramp in V
        hEff = hEffFrac * H(j);            % effective altitude for density
        [~, ~, ~, rhoEff] = atmosisa(hEff,'extended',true);
        Deff = 0.5*rhoEff*Vrms2*Cd*Aref;
        dV_drag(j,i) = k_drag*(Deff/mEff)*tb;
        %end of drag stuff
        dV_grav(j,i) = g*tb;
        dV_design(j,i) = v(j,i)+dV_grav(j,i)+dV_drag(j,i);
        Isp = Ivac;
        MR = exp(dV_design(j,i)/(g*Isp));                   %Mass Ratio
        m_prop(j,i) = m_pl*(MR-1)/(MR-(MR-1)/lambda);       %mass of propellant, kg
        m_lv(j,i) = m_prop(j,i)/lambda;                     %mass of unladen launch vehicle (prop+inert) 
        GLOW(j,i) = m_lv(j,i)+m_pl;                         %mass of vehicle with payload at lift-off, kg
    end 
end

figure;
imagesc(M, H/1000, dV_design);
cb = colorbar;
title(cb, "Design \DeltaV [m/s]");
axis xy;
xlabel("Take-over Mach #");
ylabel("Take-over Altitude [km]");

figure;
imagesc(M, H/1000, m_prop);
cb = colorbar;
title(cb, "Propellant Mass [kg]");
axis xy;
xlabel("Take-over Mach #");
ylabel("Take-over Altitude [km]");

i = find(abs(M-4) == min(abs(M-4)));
j = find(abs(H-25000) == min(abs(H-25000)));
disp(m_prop(j,i))
