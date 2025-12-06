% AAE 537 Hypersonic Propulsion
% Project Sounding Rocket Sizing 
% Tyler Nord
clear; clc; close all;

%constants/requirements
g = 9.81;                                                   %gravitational acceleration, m/s^2
m_pl = 600;                                                 %payload mass, lbf
H = linspace(10, 30, 100).*1000;                            %take over altitude, m
[~, a, ~, ~, ~, ~] = atmosisa(H, 'extended', true);         %take over speed of sound, m/s
M = linspace(2, 4, 100);                                    %take over Mach number

%assumptions taken from Black Brant V
tb = 26.53;                                                 %burn time, s
Ivac = 258.29;                                              %specific impulse, sec
lambda = 0.8;                                               %propellant mass fraction
Ae = 0.0773;                                                %nozzle exit area, m^2

for i = 1:length(M)
    for j = 1:length(a)
        v(j,i) = M(i)*a(j);                                 %take over velocity, m/s
        dV_ideal(j,i) = sqrt(v(j,i)^2+2*g*H(j));            %from conservation of energy, neglects drag, steering, and grav losses
        f = 0;                                              %time at which linear pitching begins, normalized by burn time
        dV_drag(j,i) = dV_ideal(j,i)*0.1;                   %rough assumption for drag losses
        dV_grav(j,i) = g*tb*(f+(1-f)*(2/pi));
        dV_design(j,i) = dV_ideal(j,i)+dV_grav(j,i)+dV_drag(j,i);
        Isp = Ivac;
        MR = exp(dV_design(j,i)/(g*Isp));                   %Mass Ratio
        m_prop(j,i) = m_pl*(MR-1)/(MR-(MR-1)/lambda);       %mass of propellant, kg
        m_lv(j,i) = m_prop(j,i)/lambda;                     %mass of unladen launch vehicle (prop+inert) 
        GLOW(j,i) = m_lv(j,i)+m_pl;                         %mass of vehicle with payload at lift-off, kg
    end 
end

% figure;
% imagesc(M, H/1000, v);
% cb = colorbar;
% title(cb, "Take-over Velocity [m/s]");
% axis xy;
% xlabel("Take-over Mach #");
% ylabel("Take-over Altitude [km]");

figure;
imagesc(M, H/1000, dV_ideal);
cb = colorbar;
title(cb, "Ideal \DeltaV [m/s]");
axis xy;
xlabel("Take-over Mach #");
ylabel("Take-over Altitude [km]");

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


dV_design(100,100)