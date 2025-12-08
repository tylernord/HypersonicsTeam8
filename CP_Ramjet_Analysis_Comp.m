% CP Ramjet Analysis Companion Script

%% Param Definition

thrust = 5000; % From traj sim, N
M_cruise = 3.75; % Cruise mach number
alt = 25000; % Cruise altitude, m
S = 5/(3.2808^2); % Capture area, ft^2
gamma = 1.4;
gamma_hot = 1.3;
cp_hot = 1100;

[T,a,P,rho,~,~] = atmosisa(alt,extended="on"); % Atmospheric properties at cruise

T03 = T*(1+(gamma-1)/2 * M_cruise^2); % Combustor inlet stagnation temp, K
V_cruise = M_cruise*a; % Cruise velocity, m/s
mdot_air = V_cruise*rho*S;
T04 = 1800;
P_exit = P;

M4_range = 0.5:0.01:1.0;

P3_vec = []; P03_vec = []; M3_vec = []; AR_vec = []; P04_vec = [];

for i = 1:length(M4_range)
    M4 = M4_range(i);
    [M3,P3,P03,P04,M_exit,AR_exit,comb_AR,total_AR] = CP_Ramjet_Analysis(mdot_air,thrust,P_exit,T03,T04,M4,gamma_hot,cp_hot,V_cruise);
    
    M3_vec(end+1) = M3;
    P3_vec(end+1) = P3;
    P03_vec(end+1) = P03;
    AR_vec(end+1) = total_AR;
    P04_vec(end+1) = P04;
end

%% Plotting

figure(1);
plot(M4_range,M3_vec);
xlabel("Fixed M4");
ylabel("Required M3");
title("M3 vs M4");
grid on;

figure(2);
plot(M4_range,P3_vec./1000);
xlabel("Fixed M4");
ylabel("Required P3 [kPa]");
title("P3 vs M4");
grid on;

figure(3);
plot(M4_range,P03_vec./1000);
xlabel("Fixed M4");
ylabel("Required P_{03} [kPa]");
title("P_{03} vs M4");
grid on;

figure(4);
plot(M4_range,AR_vec);
xlabel("Fixed M4");
ylabel("Area Ratio change across combustor to nozzle");
title("Area Ratio across combustor to nozzle vs M4");
grid on;