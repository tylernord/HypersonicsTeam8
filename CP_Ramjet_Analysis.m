function [M3,P3,P03,P04,M_exit,AR_exit,comb_AR,total_AR] = CP_Ramjet_Analysis(mdot_air,thrust,P_exit,T03,T04,M4,gamma_hot,cp_hot,V_cruise)

% This function computes the upstream conditions for a constant pressure
% ramjet combustor given various performance and downstream conditon inputs
% using gas dynamic equations from class and CEA analysis.

% Inputs:
% mdot_air = Air flow rate through engine, kg/s
% thrust = Total thrust produced by engine, N
% P_exit = Exit pressure of nozzle, Pa
% T03 = Stagnation temperature at combustor inlet, K
% T04 = Fixed stagnation pressure at combustor exit, K
% M4 = Fixed mach number at combustor exit
% gamma_hot = Specific heat ratio for expanding hot gas
% cp_hot = Specific heat capacity of exiting hot gas, J/kg-K
% V_cruise = Cruise velocity of ramjet, m/s

% Outputs:
% M3 = Inlet combustor mach number
% P3 = Inlet combustor static pressure, Pa
% P03 = Inlet combustor stagnation pressure, Pa
% M_exit = Nozzle exit mach number
% AR_exit = Area ratio of nozzle to produce thrust
% comb_AR = Area ratio of combustor exit to inlet

% Assumptions:
% 1. Ideal equation for thrust applies (no fuel fraction considered -->
%    more conservative in jet thrust cause minor impact)
% 2. Constant static pressure combustor
% 3. Fixed stagnation pressure rise across combustor

% Combustor inlet properties
gamma = 1.4; % Specific heat ratio at comb. inlet
R = 287; % Specific gas law constant at comb. inlet, J/kg-K
cp = R/(1 - 1/gamma); % Specific constant pressure specific heat at inlet, J/kg-K

% Flow functions file path
addpath('C:\Users\rugve\OneDrive - purdue.edu\MATLAB\AAE Coursework Scripts\Fall 2025\AAE 537\HW2');

% Combustor Exit Condition Analysis

if M4 ~= 1
    AR_comb_nozz = D_func(1,gamma_hot)/D_func(M4,gamma_hot);
else
    AR_comb_nozz = 1;
end

n_n = 0.98;
V_exit = (thrust + mdot_air*V_cruise)/mdot_air; % Ideal ramjet thrust analysis
P04 = fzero(@(P04) V_exit - (2 * cp_hot * T04 * n_n * (1 - (P_exit/P04)^((gamma_hot-1)/gamma_hot)))^0.5, P_exit*100);
M_exit = (P04/P_exit)^((gamma_hot - 1)/gamma_hot); % Engine exit mach number
AR_exit = 1/M_exit * (((gamma_hot+1)/2)/(1+(gamma_hot-1)/2 * M_exit^2))^((gamma_hot+1)/(2 - 2*gamma_hot)); % Area ratio for nozzle exit

% Combustor Analysis
R_hot = cp_hot*(1 - 1/gamma_hot); % Specific gas law constant for combustor gases, J/kg-K

P4 = P04/((1+(gamma_hot-1)/2 * M4^2)^(gamma_hot/(gamma_hot-1))); % Static pressure, combustor exit, P4
P3 = P4;
M3 = fzero(@(M3) M3/(1+(gamma-1)/2 * M3^2)^0.5 - (M4/(1+(gamma_hot-1)/2 * M4^2)^0.5) * sqrt((gamma_hot*R_hot)/(gamma*R)) * sqrt(T04/T03), M4/2);
comb_AR = (sqrt(gamma_hot/R_hot) * sqrt(R/gamma) * sqrt(T03/T04) * (1+(gamma_hot-1)/2 * M4^2)^0.5 / (M3*(1+(gamma-1)/2 * M3^2)^0.5))^-1;
P03 = P3*((1+(gamma-1)/2 * M3^2)^(gamma/(gamma-1)));

total_AR = comb_AR*(AR_comb_nozz^-1)*AR_exit; % Total area ratio change from combustor inlet to engine exit

end

