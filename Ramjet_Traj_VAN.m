% AAE 537 Project Trajectory Code

% Author(s): Rugved Dikay, ...
% Date: 11/26/2025

% Description: This script computes the trajectory of the ramjet based on
% flight properties of the ramjet (mass, wetted area, fuel fraction...) and
% project mission requirements. As such, this code can be used as both a
% design tool (informing required thrust and inlet conditions), as well as
% a verification of what the trajectory of a certain design will be.

clc;
clear;
close all;

% Flow functions file path
addpath('C:\Users\rugve\OneDrive - purdue.edu\MATLAB\AAE Coursework Scripts\Fall 2025\AAE 537\HW2');
% CEA file path
addpath("C:\Users\rugve\OneDrive - purdue.edu\MATLAB\Useful Scripts\Purdue CEA")

%% Assumptions

% As of 12/06/2025

% 1. Cruise is variable altitude, constant mach (unaccelerated flight in cruise direction)
%      - This means lift of the vehicle changes as a function of time, as
%        vehicle AoA is not changed. The net result is a slight increase in
%        altitude for the vehicle over time, which can be controlled via
%        changing the inert mass of the ramjet.
% 2. Fuel fraction is iteratively determined for each time step to reflect
%    thrust needed for unaccelerated cruise, and prescribed acceleration
%    during the boost phase.
% 3. Drag and lift correlations have been taken from HW1 and a NASA report
%    on all body hypersonic vehicle analysis (assuming wave-riderish
%    geometry)
% 4. Inlet area assumed to be a certain fraction of wetted area
% 5. Perfectly expanded, VARIABLE AREA nozzle (done to determine how
%    significant variation in nozzle exit must be, to see if its a design
%    option)
% 6. Desired combustor inlet pressure computed --> Inform inlet design
%    instead of other way around
% 7. Realistic combustor and nozzle efficiencies considered
% 8. Trajectory analysis is split into an acceleration phase post
%    separation from rocket and a cruise phase.
% 9. For simplicity in iteration/calcs, Cp and gamma at nozzle exit
%    are conisdered as reasonable constants (based on prior
%    estimates/studies)
% 10.This assumes engine is oriented in-line with body and so sees flow
%    with same incidence as alpha compared to flight theta.
% 11. Combustion analysis is done by considering mach number at combustor
%     exit with diff. area ratios of nozzle throat to combustor area. This
%     method of calculation is more scalable for scramjet type device.
% 12. Nozzle throat is assumed choked --> May have to rework this for a
%     DMRJ configuration.
      

%% Param Definition

M0_range = 4.5;
GLOW_range = 200:25:600; % Range of ramjet weights to evaluate effect of lift and inlet crossflow, lbs
stoich_JETA = 15; % Stoichometric air-fuel mass ratio for Jet-A and air

for i = 1:1 % Performs computations for range of Mach numbers
    for j = 17:17 % Performs analysis for different ramjet weights

        % Flight conditions --> Update these !!!
        M0 = M0_range(i); % Design flight mach number
        alt = 25000; % Cruise altitude, m
        release_M = 3.75; % Separation mach number
        release_alt = 25000; % Separation altitude, m
        g = 9.81; % Gravity in m/s^2
        gamma = 1.4; % Specific heat ratio for air before combustor
        beta = @(M) (M^2-1)^0.5; % Air compressibility factor function of mach #

        % Vehicle properties --> Update these !!!
        GLOW = GLOW_range(j); % Initial weight of ramjet, lbs
        S = 10; % Wetted surface area for lift/drag, ft^2
        A_cap_frac = 0.50; % Inlet area/Wetted area
        A_cap = S*A_cap_frac; % Capture area for engine
        AspectRatio = 1.5; % Aerodynamic aspect ratio of vehicle design
        engineOffset = 0.0; % Aeroydnamic offset of engine, degress

        % Trajectory properties --> Update these!!!
        sep_type = "constant alt"; % Trajectory type between sep. and cruise start
        range = 250; % Cruise range, miles

        % Engine properties --> Update these!!!

        % Ideally we have a couple cases:
        % 1. Fuel flow rate is constant for cruise, meaning expansion system and
        %    inlet have to adjust
        % 2. Variable fuel flow rate: Need CEA, nozzle exit conditions, ISP
        % Maybe some type of convergence study
        % Currently using f value and inlet mdot_air assumptions

        % P0_recov = 0.80; % Conservative value for inlet pressure recovery, %
        % T04_max = 1800; % Maximum combustor stagnation temperature, K --> Update this based on lit review!!!
        T04_design = 1700; % Design combustor exit stagnation temperature, K during accel. (Ramjet mode?)
        AR_45 = 1.25; % Area ratio between combustor and nozzle throat --> Edit this parameter!
        n_b = 0.95; % Combustion efficiency
        n_n = 0.98; % Nozzle efficiency
        cp_hot = 1100; % Conservative value for specific heat based on CEA iteration, J/kg-K
        gamma_hot = 1.30; % Conservative value for specific heat ratio based on CEA iteration

        %% Computation

        % Conversions
        range = range*5280/3.2808; % Range in km
        tot_mass = GLOW/2.205; % Vehicle mass in kg
        S = S/(3.2808^2); % Wetted drag area in m^2
        A_cap = A_cap/(3.28082^2); % Inlet area, m^2

        % Flight properties at release
        curr_m = tot_mass;
        init_prop = air_prop(release_alt);
        a = init_prop(2);
        V_release = release_M*a;

        % Flight properties at cruise
        props = air_prop(alt);
        a = props(2);
        V_cruise = M0*a;

        % Data vectors

        % Initialized
        mass_vec = curr_m; % Vehicle mass, kg
        az_vec = 0; % Vehicle vertical acceleration, m/s^2
        Vz_vec = 0; % Vehicle vertical velocity, m/s
        alt_vec = release_alt; % Vehicle altitude, m
        Vx_vec = V_release; % Vehicle horizontal velocity, m/s
        range_vec = 0; % Vehicle range, m
        t_vec = 0; % Time vector, s.
        incidence_vec = [0];

        % Non-initialized
        P00_vec = []; % Freestream stagnation pressure, Pa
        thrust_vec = []; % Vehicle thrust, N
        f_vec = []; % Engine fuel fraction
        T04_vec = []; % Engine stagnation temperature, K
        phi_vec = []; % Equivalence ratio
        SFC_vec = []; % Engine SFC
        LD_vec = []; % Lift drag ratio
        rho_vec = []; % Density, kg/m^3
        drag_vec = []; % Drag force, N
        alpha_vec = []; % Angle of attack of control surfaces, deg.
        theta_vec = []; % Velocity vector angle of vehicle, deg.
        %incidence_vec = []; % Incidence vector for vehicle inlet, deg.
        AR_vec = []; % Area ratio for nozzle
        P03_vec = []; % Stagnation pressure at combustor inlet, Pa
        P04_vec = []; % Stagnation pressure at combustor outlet, Pa
        M3_vec = []; % Mach number at combustor inlet
        M4_vec = []; % Mach number at combustor outlet
        M_exit_vec = []; % Mach number at engine exit

        %% Simulating Ramjet Release to Cruise (R2C)

        fprintf("\nNow starting Release to Cruise Sim:\n");

        dt_sep = 1; % Time step for R2C sim
        t_accel = 20; % Total acceleration time to cruise, s
        dVx = (V_cruise - V_release)/t_accel; % Acceleration of ramjet, m/s^2
        %V_profile = accel_profile(V_sep,V_cruise,t_accel,dt_sep);
        %V_profile = V_release:(V_cruise-V_release)/(t_accel/dt_sep):V_cruise;
        
        while Vx_vec(end) < V_cruise % While ramjet is accelerating

            % Current properties
            Vx = Vx_vec(end);
            Vz = Vz_vec(end);
            t = t_vec(end);
            alt = alt_vec(end);

            fprintf("\nCurrently assessing time step %0.f", t);

            props = air_prop(alt);
            T0 = props(1);
            a = props(2);
            P0 = props(3);
            rho = props(4);

            % Flight parameters
            q_R2C = 0.5*rho*Vx^2; % Dynamic pressure at separation
            M_R2C = Vx/a; % Mach number
            theta_vec(end+1) = atand(Vz/Vx); theta = theta_vec(end); % Vehicle flight path angle

            % Acceleration
            Vx_vec(end+1) = Vx + dVx*dt_sep; % Velocity of ramjet at next time step

            % Lift, Drag, and Thrust Analysis
            Cd1 = Cd_mach(M_R2C); % Based on class given correlations
            if Vz == 0
                Cl_R2C = 0.08; % Lift coefficient per NASA documents
            elseif Vz < 0 && incidence_vec(end) < 0
                Cl_R2C = Cl_R2C + 0.015;
            elseif Vz > 0 && incidence_vec(end) > 0
                Cl_R2C = Cl_R2C - 0.015;
            elseif Vz < 0 && incidence_vec(end) > 0
                Cl_R2C = Cl_R2C + 0.005;
            elseif Vz > 0 && incidence_vec(end) < 0
                Cl_R2C = Cl_R2C - 0.005;    
            end
            [Cd2,alpha,~] = NASA_coeff(M_R2C,beta(M_R2C),Cl_R2C,AspectRatio,0);
            if Cd1 > Cd2 % Choose larger drag coefficient for R2C phase
                Cd = Cd1;
            else
                Cd = Cd2;
            end
            lift = q_R2C*S*Cl_R2C; % Lift force generated, N (should counter weight)
            drag = q_R2C*S*Cd; % Drag force generated, N
            thrust = (curr_m*dVx + drag)/cosd(theta + engineOffset); % Thrust force needed to accelerate, N
            incidence_vec(end+1) = theta + engineOffset; % Incidence angle for engine, deg.

            % Station Properties/Cycle Analysis
            mdot_air = rho*A_cap*Vx_vec(end); % Engine air mdot, kg/s

            P00 = P0*(1+(gamma-1)/2 * M_R2C^2)^(gamma/(gamma-1)); % Inlet stag. pressure
            T00 = T0*(1+(gamma-1)/2 * M_R2C^2); % Inlet stagnation temp.

            T03 = T00; % Stag. temperature at combustor inlet, K
            %T3 = T03/(1+(gamma-1)/2 * M3^2); % Static temperature, combustor inlet, K

            P_exit = P0; % Perfect expansion

            V_exit = (thrust + mdot_air*V_cruise)/mdot_air; % Ideal ramjet thrust analysis
            P04 = fzero(@(P04) V_exit - (2 * cp_hot * 1800 * n_n * (1 - (P_exit/P04)^((gamma_hot-1)/gamma_hot)))^0.5, P00*0.5);
            M_exit = (P04/P_exit)^((gamma_hot - 1)/gamma_hot); % Engine exit mach number
            AR_vec(end+1) = 1/M_exit * (((gamma_hot+1)/2)/(1+(gamma_hot-1)/2 * M_exit^2))^((gamma_hot+1)/(2 - 2*gamma_hot)); % Area ratio for nozzle exit
            
            D5 = D_func(1,gamma_hot); D_exit = D_func(M_exit,gamma_hot); AR = D5/D_exit; % Area ratio of nozzle
            M4 = fzero(@(M) D_func(1,gamma_hot) - AR_45*D_func(M,gamma_hot),1/M0); % Mach number at combustor outlet, subsonic
            M3 = fzero(@(M) sqrt(T04_design/T03) - N_func(M4,gamma_hot)/N_func(M,gamma_hot), 1/M0); % Mach number at combustor inlet, subsonic
            P03 = fzero(@(P03) P04*G_func(M4,gamma_hot) - P03*G_func(M3,gamma_hot), P00*0.7); % Stagnation ressure at combustor inlet
            T3 = T03/(1+(gamma-1)/2 * M3^2); % Static temperature at combustor inlet

            if size(f_vec,2) < 1
                f_guess = 0;
            else
                f_guess = f_vec(end)*n_b;
            end
            [f,~,~] = engineIterator(P03,T03,T3,T04_design,f_guess);

            % Property + Vector updates
            f = f/n_b; % Accounting for combustion efficiency
            mdot_fuel = mdot_air*f; %Consumed fuel mass, kg

            % Vector updates
            f_vec(end+1) = f;
            thrust_vec(end+1) = thrust;
            range_vec(end+1) = Vx*dt_sep + range_vec(end);
            T04_vec(end+1) = T04_design;
            phi_vec(end+1) = stoich_JETA/(1/f);
            SFC_vec(end+1) = (mdot_fuel*3600)/thrust;
            az_vec(end+1) = (lift+thrust*sind(theta + engineOffset)-curr_m*g)/curr_m; % Changing AoA to keep flight level
            Vz_vec(end+1) = Vz_vec(end) + az_vec(end)*dt_sep;
            alt_vec(end+1) = alt + Vz_vec(end)*dt_sep;
            LD_vec(end+1) = lift/drag;
            rho_vec(end+1) = rho;
            drag_vec(end+1) = drag;
            alpha_vec(end+1) = alpha;
            P00_vec(end+1) = P00;
            P03_vec(end+1) = P03;
            P04_vec(end+1) = P04;
            M3_vec(end+1) = M3;
            M4_vec(end+1) = M4;
            M_exit_vec(end+1) = M_exit;

            curr_m = curr_m - dt_sep*mdot_fuel; % Update vehicle mass, kg
            mass_vec(end+1) = curr_m;
            t_vec(end+1) = t_vec(end) + dt_sep;
            disp(alt);

        end

        %%
        save('C:\Users\rugve\OneDrive - purdue.edu\MATLAB\AAE Coursework Scripts\Fall 2025\AAE 537\Project\Trajectory Stuff\Release_to_Cruise_Sim_VAN.mat','t_vec','thrust_vec','T04_vec','f_vec','phi_vec','SFC_vec','alpha_vec','alt_vec','drag_vec','LD_vec','rho_vec','mass_vec','M3_vec','AR_vec','P03_vec','theta_vec','alpha_vec','incidence_vec');
        input("enter");
       
        disp("\nNow starting cruise sim:\n");

        % Cruise simulation
        t_cruise_end = range/V_cruise;
        dt = 1; % Time step for cruise sim, s.
        Cl_cruise = Cl_R2C;
        % dt = round(t_cruise_end/100); % Useful dt for long cruise sim.

        % % Cruise aerodynamics
        % lift_cruise = curr_m*g; % Lift force at start of cruise = weight for SLUF
        % Cd1 = Cd_mach(M0); % Drag ceofficient in cruise based on class given correlations
        % Cl_design = (lift_cruise)/(q_R2C*S); % Design lift coefficient throughout cruise
        % [Cd2,alpha,~] = NASA_coeff(M0,beta(M0),Cl_design,AR,0);
        % if Cd1 > Cd2 % Choose larger drag coefficient for cruise
        %     Cd_cruise = Cd1;
        % else
        %     Cd_cruise = Cd2;
        % end

        while range_vec(end) < range % Run for entire cruise range
            % Current state values
            fprintf("\nCurrent sim time is: %0.f seconds", t);
            alt = alt_vec(end);
            Vz = Vz_vec(end);

            % Atmospheric properties
            props = air_prop(alt);
            T0 = props(1);
            a = props(2);
            P0 = props(3);
            rho = props(4);
            V_cruise = M0*a;
            q_cruise = 0.5*rho*V_cruise^2; % Cruise dynamic pressure
            theta_vec(end+1) = atand(Vz/V_cruise); theta = theta_vec(end); % Vehicle flight path angle
            
            % Lift, Drag, and Thrust Analysis
            Cd1 = Cd_mach(M_R2C); % Based on class given correlations
            if Vz == 0
                Cl_cruise = 0.08; % Lift coefficient per NASA documents
            elseif Vz < 0 && incidence_vec(end) < 0
                Cl_cruise = Cl_cruise + 0.015;
            elseif Vz > 0 && incidence_vec(end) > 0
                Cl_cruise = Cl_cruise - 0.015;
            elseif Vz < 0 && incidence_vec(end) > 0
                Cl_cruise = Cl_cruise + 0.005;
            elseif Vz > 0 && incidence_vec(end) < 0
                Cl_cruise = Cl_cruise - 0.005;
            end
            [Cd2,alpha,~] = NASA_coeff(M0,beta(M0),Cl_cruise,AspectRatio,0);
            if Cd1 > Cd2 % Choose larger drag coefficient for R2C phase
                Cd = Cd1;
            else
                Cd = Cd2;
            end
            lift = q_cruise*S*Cl_cruise;
            drag = q_cruise*S*Cd;
            thrust = drag/cosd(theta + engineOffset);
            incidence_vec(end+1) = theta + engineOffset; % Incidence angle for engine, deg.
            disp(incidence_vec(end));

            % Station Properties/Cycle Analysis
            mdot_air = rho*A_cap*Vx_vec(end); % Engine air mdot, kg/s

            P00 = P0*(1+(gamma-1)/2 * M_R2C^2)^(gamma/(gamma-1)); % Inlet stag. pressure
            T00 = T0*(1+(gamma-1)/2 * M_R2C^2); % Inlet stagnation temp.

            T03 = T00; % Stag. temperature at combustor inlet, K
            %T3 = T03/(1+(gamma-1)/2 * M3^2); % Static temperature, combustor inlet, K

            P_exit = P0; % Perfect expansion

            V_exit = (thrust + mdot_air*V_cruise)/mdot_air; % Ideal ramjet thrust analysis
            P04 = fzero(@(P04) V_exit - (2 * cp_hot * 1800 * n_n * (1 - (P_exit/P04)^((gamma_hot-1)/gamma_hot)))^0.5, P00*0.5);
            M_exit = (P04/P_exit)^((gamma_hot - 1)/gamma_hot); % Engine exit mach number
            AR_vec(end+1) = 1/M_exit * (((gamma_hot+1)/2)/(1+(gamma_hot-1)/2 * M_exit^2))^((gamma_hot+1)/(2 - 2*gamma_hot)); % Area ratio for nozzle exit
            
            D5 = D_func(1,gamma_hot); D_exit = D_func(M_exit,gamma_hot); AR = D5/D_exit; % Area ratio of nozzle
            M4 = fzero(@(M) D_func(1,gamma_hot) - AR_45*D_func(M,gamma_hot),1/M0); % Mach number at combustor outlet, subsonic
            M3 = fzero(@(M) sqrt(T04_design/T03) - N_func(M4,gamma_hot)/N_func(M,gamma_hot), 1/M0); % Mach number at combustor inlet, subsonic
            P03 = fzero(@(P03) P04*G_func(M4,gamma_hot) - P03*G_func(M3,gamma_hot), P00*0.7); % Stagnation ressure at combustor inlet
            T3 = T03/(1+(gamma-1)/2 * M3^2); % Static temperature at combustor inlet

            f_guess = f_vec(end)*n_b;
            [f,~,~] = engineIterator(P03,T03,T3,T04_design,f_guess);

            % Property + Vector updates
            f = f/n_b; % Accounting for combustion efficiency
            t = t + dt; % Cruise time update
            mdot_fuel = mdot_air*f; %Consumed fuel mass, kg

            % Vector updates
            f_vec(end+1) = f;
            thrust_vec(end+1) = thrust;
            range_vec(end+1) = V_cruise*dt + range_vec(end);
            T04_vec(end+1) = T04_design;
            phi_vec(end+1) = stoich_JETA/(1/f);
            SFC_vec(end+1) = (mdot_fuel*3600)/thrust;
            az_vec(end+1) = (lift+thrust*sind(theta + engineOffset)-curr_m*g)/curr_m; 
            Vz_vec(end+1) = Vz_vec(end) + az_vec(end)*dt;
            alt_vec(end+1) = alt + Vz_vec(end)*dt;
            Vx_vec(end+1) = V_cruise;
            LD_vec(end+1) = lift/drag;
            rho_vec(end+1) = rho;
            drag_vec(end+1) = drag;
            alpha_vec(end+1) = alpha;
            P00_vec(end+1) = P00;
            P03_vec(end+1) = P03;
            P04_vec(end+1) = P04;
            M3_vec(end+1) = M3;
            M4_vec(end+1) = M4;
            M_exit_vec(end+1) = M_exit;

            curr_m = curr_m - dt*mdot_fuel; % Update vehicle mass, kg
            mass_vec(end+1) = curr_m;
            t_vec(end+1) = t_vec(end) + dt;

        end

        save('C:\Users\rugve\OneDrive - purdue.edu\MATLAB\AAE Coursework Scripts\Fall 2025\AAE 537\Project\Trajectory Stuff\Full_Trajectory_Sim_VAN.mat','t_vec','thrust_vec','T04_vec','f_vec','phi_vec','SFC_vec','alpha_vec','alt_vec','drag_vec','LD_vec','rho_vec','mass_vec','M3_vec','AR_vec','P03_vec','theta_vec','alpha_vec','incidence_vec');

        %% Plotting

        fprintf("\n");
        figure(10*j-9);
        plot(t_vec./60, mass_vec.*2.205);
        xlabel("Time [min]");
        ylabel("Ramjet Mass [lbs]");
        title("Ramjet Mass vs. Time in Cruise");
        grid on;

        figure(10*j-8);
        plot(t_vec(1:end-1)./60, thrust_vec);
        xlabel("Time [min]");
        ylabel("Ramjet Thrust [N]");
        title("Ramjet Thrust vs. Time in Cruise");
        grid on;

        figure(10*j-7);
        plot(range_vec.*3.2808/5280,alt_vec./1000);
        xlabel("Range [miles]");
        ylabel("Altitude [km]");
        title("Vehicle Altitude vs. Range in Cruise");
        ylim([0,max(alt_vec./1000)*1.5]);
        grid on;

        figure(10*j-6);
        plot(t_vec(1:end-1)./60,f_vec);
        xlabel("Time [min]");
        ylabel("Fuel Fraction");
        title("Vehicle Fuel Fraction vs. Time in Cruise");
        grid on;

        figure(10*j-5);
        plot(t_vec(2:end)./60,LD_vec);
        xlabel("Time [min]");
        ylabel("Lift-Drag Ratio");
        title("Vehicle L/D vs. Time in Cruise");
        ylim([0,5]);
        grid on;

        figure(10*j-4);
        plot(t_vec(2:end)./60,thrust_vec./(mdot_air*f_vec*g));
        xlabel("Time [min]");
        ylabel("Ramjet Specific Impulse [s]");
        title("Ramjet ISP vs Time");
        grid on;

        figure(10*j-3);
        plot(t_vec(2:end)./60,rho_vec);
        xlabel("Time [min]");
        ylabel("Freestream Density [kg/m^3]");
        title("\rho_0 vs. Time");
        grid on;

        figure(10*j-2);
        plot(t_vec(2:end)./60,drag_vec);
        xlabel("Time [min]");
        ylabel("Drag Force [N]");
        title("Drag vs. Time");
        grid on;

        figure(10*j-1);
        plot(t_vec(2:end)./60,phi_vec);
        xlabel("Time [min]");
        ylabel("Equivalence Ratio");
        title("Phi vs. Time");
        grid on;

        figure(1000);
        plot(t_vec(2:end)./60,AR_vec);
        xlabel("Time [min]");
        ylabel("Area Ratio of Nozzle Exit / Throat");
        title("Nozzle Exit Area Ratio vs Time");
        grid on;

        figure(1001);
        plot(t_vec(2:end)./60,P03_vec./P00_vec);
        xlabel("Time [min]");
        ylabel("Inlet Compressive Efficiency");
        title("Required \eta_c vs. Time");
        grid on;

        figure(1002);
        plot(t_vec(2:end)./60,alpha_vec);
        xlabel("Time [min]");
        ylabel("Vehicle Angle of Attack [deg]");
        title("Vehicle Angle of Attack vs Time");
        grid on;

        figure(1003);
        plot(t_vec./60,incidence_vec);
        xlabel("Time [min]");
        ylabel("Angle of Flow Incidence to Engine Inlet [deg]");
        title("Engine Flow Incidence Angle vs. Time");
        grid on;

        figure(1004);
        plot(t_vec(2:end)./60,theta_vec);
        xlabel("Time [min]");
        ylabel("Vehicle Flight Path Angle [deg]");
        title("Vehicle Flight Path Angle vs. Time");
        grid on;

        figure(1005);
        plot(t_vec(2:end)./60,M3_vec);
        xlabel("Time [min]");
        ylabel("Combustor Inlet Mach Number");
        title("Combustor Inlet Mach Number vs. Time");
        grid on;

        figure(1006);
        plot(t_vec(2:end)./60,M3_vec);
        xlabel("Time [min]");
        ylabel("Combustor Inlet Mach Number");
        title("Combustor Inlet Mach Number vs. Time");
        grid on;

        figure(1007);
        plot(t_vec(2:end)./60,M4_vec);
        xlabel("Time [min]");
        ylabel("Combustor Exit Mach Number");
        title("Combustor Exit Mach Number vs. Time");
        grid on;


    end
end

% % Final figure formatting
% figure(1000);
% hold off;
% cb = colorbar;
% xlabel("Flight Mach Number");
% ylabel("Altitude [km]");
% title("Maximum Thrust in Cruise");
% title(cb,"Maximum Thrust [N]");
% grid on;
% 
% figure(1001);
% hold off;
% cb = colorbar;
% xlabel("Flight Mach Number");
% ylabel("Altitude [km]");
% title("ISP at Max Thrust");
% title(cb,"Specific Impulse [s]");
% grid on;

%% Helper Functions

function forces = LD_compute(Cd, V, rho, S)
% This function computes the magnitude of lift and drag force on the
% vehicle based on the vehicle dynamic pressure, Cd, and surface area.

forces = zeros(2,1);

drag = 0.5*S*Cd*rho*V^2; % Drag force, N
lift = drag*(3*exp(-0.06)); % Lift force, N

forces(1) = lift;
forces(2) = drag;
end

function Cd = Cd_mach(M)
% This function computes the drag coefficient (Cd) of the vehicle as
% a function of mach number, based on provided correlations

if M <= 0.4
    Cd = 0.03528;
elseif M > 0.4 && M <= 1.5
    Cd = 0.04556*M^2 - 0.06443*M + 0.05376;
else
    Cd = 0.1077*M^(-1.484);
end
end

function X0 = stag_prop(static_prop, prop_type, mach)
% This function computes the stagnation temperature or pressure seen by
% the vehicle as a function of mach number. Air is assumed here as the
% ideal gas media.

gamma = 1.4; % Ratio of specific heats for air

if prop_type == "P"
    X0 = static_prop * (1 + (gamma-1)/2 * mach^2) ^ (gamma/(gamma-1));
else
    X0 = static_prop * (1 + (gamma-1)/2 * mach^2);
end
end

function properties = air_prop(alt)
% This function computes the properties of the freestream air seen by
% the vehicle given an altitude at which it is flying.

properties = zeros(4,1);

[T,a,P,rho] = atmosisa(alt,extended="on");
properties(1) = T; % Temperature, K
properties(2) = a; % Speed of sound, m/s
properties(3) = P; % Pressure, Pa
properties(4) = rho; % Density, kg/m^3
end

function output_data = air_prop_adv(P,T)
% This function uses NASA CEA to determine certain air properties
% relevant to heat transfer analysis. Function is defined as
% air_prop_adv to reflect advanced properties being provided like
% viscosity, Prandtl number, ect

% Inputs:
% P = Fluid pressure, Pa
% T = Fluid temperature, K

% Outputs:
% output_data: Struct with all output information from simulation

CEA_RUN = true;
CEA_SAVE_FILE = 'cea.mat';

inp = containers.Map;
inp('type') = 'eq';              % Sets the type of CEA calculation
inp('p') = P/6894.76;                % Chamber pressure
inp('p_unit') = 'psi';              % Chamber pressure units
inp('sup') = 1;               % Supersonic area ratios
inp('ox') = 'Air';              % Ox name from thermo.inp
inp('ox_t') = T;                  % Ox inlet temperature
inp('file_name') = 'Scramjet_AirProps1.inp';    % Input/output file name

if CEA_RUN
    data = cea_rocket_run(inp);     % Call the CEA MATLAB code
    save(CEA_SAVE_FILE, 'data');
    output_data = data('eq');
else
    load(CEA_SAVE_FILE);
end

end

function output_data = engineCEA(f,P03,T3)
% This function performs a CEA analysis of a airbreathing combustor to
% determine the output temperature and other relevant properties
% assuming adiabatic combustion and frozen equilibrium

% Inputs:
% f: Fuel/Ox mass flow ratio for combustion
% P03: Combustor stagnation pressure
% T3: Combustor inlet air temperature

% Outputs:
% output_data: Struct with all output information from simulation

CEA_RUN = true;
CEA_SAVE_FILE = 'cea.mat';

inp = containers.Map;
inp('type') = 'eq fr';              % Sets the type of CEA calculation
inp('p') = P03/6894.76;                % Chamber pressure
inp('p_unit') = 'psi';              % Chamber pressure units
inp('o/f') = 1/f;               % Mixture ratio
inp('sup') = 1;               % Supersonic area ratios
inp('fuel') = {'Jet-A(L)'};             % Fuel name from thermo.inp
inp('fuel_t') = 298;                % Fuel inlet temperature, K
inp('ox') = 'Air';              % Ox name from thermo.inp
inp('ox_t') = T3;                  % Ox inlet temperature
inp('file_name') = 'Project_Ramjets_YesBoss.inp';    % Input/output file name

if CEA_RUN
    data = cea_rocket_run(inp);     % Call the CEA MATLAB code
    save(CEA_SAVE_FILE, 'data');
    output_data = data('eq');
else
    load(CEA_SAVE_FILE);
end

end

function [f,cp_hot,gamma_hot] = engineIterator(P03,T03,T3,T04,f_fir_guess)

% This function uses the bisection method to determine the fuel-air
% fraction needed in combustion to produce a desired rise in stagnation
% pressure

T04_calc = T03; % Initialize temperature guess as compressor exit temp, K
f_high = 0.15; % Large guess for fuel fraction
f_low = 0.0001; % Small guess for fuel fraction

counter = 0; % Variable to help expedite iteration for f_guess

while abs(T04_calc - T04)/T04 > 0.01

    if counter == 0 && f_fir_guess ~= 0
        f_guess = f_fir_guess;
    else
        f_guess = (f_high + f_low)/2;
    end

    data = engineCEA(f_guess,P03,T3);
    temps = squeeze(data('t'));
    T04_calc = temps(1); % Extracts combustor chamber temp,

    if T04_calc > T04
        f_high = f_guess;
    else
        f_low = f_guess;
    end

    % fprintf("\n    f = %f\n", f_guess);
    % fprintf("\n    f_high = %f\n", f_guess_high);
    % fprintf("\n    f_low = %f\n", f_guess_low);
    % fprintf("\n  T04_calc = %f\n", T04_calc);

    counter = counter + 1;

end

cp_vals = squeeze(data('cp')); cp_hot = cp_vals(1);
gammas = squeeze(data('gammas')); gamma_hot = gammas(1);
f = f_guess;

end

function [Cd,alpha,C_lift] = NASA_coeff(M0,beta,C_lift,AR,alpha)

% This function computes the vehicle drag coefficient and
% angle of attack based on the correlations provided in the NASA
% Estimated Aerodynamics of All-Body Configurations. This function works in
% 2 different ways, computing alpha and Cd given a Cl, and Cd and Cl given an
% alpha.

% Computing coefficients used in lift coeff. correlation
if beta < 4/AR
    c1 = pi*AR/2 - 0.153*beta*AR^2;
    c2 = interp1([0,4/AR],[0,exp(0.955-(4.35/M0))],beta);
else
    c1 = 4.17/beta - 0.13;
    c2 = exp(0.955 - (4.35/M0));
end

% Computes either alpha or lift coefficient based on given info
CL_0 = 0.02; % Lift coefficient at zero angle of attack
if C_lift ~= 0
    alpha = rad2deg(fzero(@(alpha) C_lift - c1*sin(alpha) - c2*(sin(alpha))^2 - CL_0, deg2rad(5)));
else
    C_lift = c1*sind(alpha) + c2*(sind(alpha))^2 + CL_0;
end

% Induced drag
K_M = 1.0;
C_di = K_M*C_lift*tand(alpha);

% Base drag of body
C_db = 1/(0.91*M0^2 - 0.2*M0 + 1.51);

% Total drag coefficient
Cd = C_db + C_di;

end

%% Archive

% % Seperation to Cruise Simulation
% curr_mass = GLOW; % Mass of ramjet at separation
% T04_vec_sep = [];
% X_vec = [0]; Z_vec = [release_alt];
% f_vec_sep = [];
% if release_alt > release_alt
%     trans_alt = release_alt; % Initialize altitude for transition sim
%     trans_M = release_M; % Initialize mach number for transition sim
%     psi = 40; % Initial flight path angle, degrees
%     psi_vec = psi;
%
%     % Atmospheric properties at separation
%     curr_prop = air_prop(trans_alt);
%     a = curr_prop(2);
%     P = curr_prop(3);
%     rho = curr_prop(4);
%     V_tot = trans_M*a;
%     q = 0.5*rho*V_tot^2; % Dynamic pressure to be maintained
%
%     if sep_type == "constant q"
%         while Z_vec(end) < alt % Before desired altitude reached
%
%             % Air properties
%             curr_prop = air_prop(Z_vec(end));
%             T0 = curr_prop(1);
%             a = curr_prop(2);
%             P = curr_prop(3);
%             rho = curr_prop(4);
%
%             psi = psi_vec(end);
%
%             % Body frame velocities
%             V_tot(end+1) = sqrt(2*q/rho);
%             curr_M = V_tot(end)/a;
%
%             % Inertial frame velocites
%             Vx = V_tot(end)*cosd(psi);
%             Vz = V_tot(end)*sind(psi);
%
%             % Lift and drag on vehicle
%             beta = (curr_M^2-1)^0.5;
%             [Cd,alpha,C_lift] = NASA_coeff(curr_M,beta,AR,alpha+psi);
%             L = 0; D = q*S*Cd; % Assume no lift generated in this phase
%
%             % Total vehicle acceleration
%             dVtot = (V_tot(end) - V_tot(end-1))/dt;
%
%             % Required thrust
%             T = curr_mass*dVtot + D + curr_mass*g*sind(psi_vec(end)); % Required thrust
%
%             % Inertial accelerations
%             dVx = (T*cosd(alpha+psi) - D*cosd(alpha) - L*sind(alpha))/curr_mass;
%             dVz = (T*sind(alpha+psi) + L*cosd(alpha) - D*sind(alpha) - curr_mass*g)/curr_mass;
%
%             Vx = Vx + dVx*dt; Vz = Vz + dVz*dt;
%             X_vec(end+1) = X_vec(end) + Vx*dt; Z_vec(end+1) = Z_vec(end) + Vz*dt;
%             psi_vec(end+1) = atand(Vz/Vx);
%
%             % Station Properties - Engine Analysis
%             mdot_air = rho*S*A_cap_frac*V_tot; % Engine air mdot
%
%             P00 = P*(1+(gamma-1)/2 * curr_M^2)^(gamma/(gamma-1)); % Inlet stag. pressure
%             T00 = T0*(1+(gamma-1)/2 * curr_M^2); % Inlet stagnation temp.
%
%             P03 = P00*P0_recov; % Stag. pressure at combustor inlet
%             T03 = T00; % Stag. temperature at combustor inlet, K
%             T3 = T03/(1+(gamma-1)/2 * M3^2); % Static temperature, combustor inlet, K
%             P04 = P03*P0_loss; % Stag. pressure at combustor oulet
%
%             P_exit = P;
%             thrust_calc = 0;
%             counter = 0;
%             T04_low = T03;
%             T04_high = T03*3;
%             while abs(T - thrust_calc)/T > 0.05
%
%                 if counter == 0 & size(T04_vec_sep,2) ~= 0
%                     T04_guess = T04_vec_sep(end);
%                 else
%                     T04_guess = (T04_low + T04_high)/2;
%                 end
%
%                 if size(f_vec,2) ~= 0
%                     f_fir_guess = f_vec_sep(end)/n_b;
%                 else
%                     f_fir_guess = 0;
%                 end
%
%                 [f,cp_hot,gamma_hot] = engineIterator(P03,T03,T3,T04_guess,f_fir_guess);
%
%                 % Engine exit velocity
%                 V_exit = (2 * cp_hot * T04_guess * n_n * (1 - (P_exit/P04)^((gamma_hot-1)/gamma_hot)))^0.5;
%
%                 thrust_calc = mdot_air*(1+f)*V_exit - mdot_air*V_tot(end);
%
%                 if thrust_calc > T
%                     T04_high = T04_guess;
%                 else
%                     T04_low = T04_guess;
%                 end
%
%                 %fprintf("\nthrust_delta = %f\n", abs(thrust_calc-thrust));
%                 counter = counter + 1;
%             end
%
%             % Property + Vector updates
%             f = f/n_b; % Accounting for combustion efficiency
%             f_vec_sep(end+1) = f;
%             T04_vec_sep(end+1) = T04_guess;
%
%             curr_mass = curr_mass - mdot_air * f;
%         end
%     end
% end

% thrust_calc = 0;
%             T04_low = T03;
%             T04_high = T04_max;
%             fprintf("\nCurrent sim time is: %0.f seconds", t);
%             counter = 0; % Controls guessing in T04 binary search
% 
%             while abs(thrust - thrust_calc)/thrust > 0.05
% 
%                 if counter == 0 & size(T04_vec,2) ~= 0
%                     T04_guess = T04_vec(end);
%                 else
%                     T04_guess = (T04_low + T04_high)/2;
%                 end
% 
%                 % Use iterator to find combustor properties
%                 [f,cp_hot,gamma_hot] = engineIterator(P03,T03,T3,T04_guess);
% 
%                 % Engine exit velocity & thrust --> Converging/Diverging nozzle
%                 % analysis
%                 V_exit = (2 * cp_hot * T04_guess * n_n * (1 - (P_exit/P04)^((gamma_hot-1)/gamma_hot)))^0.5;
%                 thrust_calc = mdot_air*(1+f)*V_exit - mdot_air*V_cruise;
% 
%                 % Update T04 guess
%                 if thrust_calc > thrust
%                     T04_high = T04_guess;
%                 else
%                     T04_low = T04_guess;
%                 end
% 
%                 counter = counter + 1;
%             end
            % 
            % % Iterative solver for combustor exit conditions
            % thrust_calc = 0;
            % T04_low = T03;
            % T04_high = T04_max;
            % fprintf("\nCurrent sim time is: %0.f seconds", t);
            % 
            % while abs(thrust - thrust_calc)/thrust > 0.05
            % 
            %     T04_guess = (T04_high + T04_low)/2;
            % 
            %     [f,cp_hot,gamma_hot] = engineIterator(P03,T03,T3,T04_guess,0);
            % 
            %     % Engine exit velocity & thrust --> Converging/Diverging nozzle analysis
            %     V_exit = (2 * cp_hot * T04_guess * n_n * (1 - (P_exit/P04)^((gamma_hot-1)/gamma_hot)))^0.5;
            %     thrust_calc = mdot_air*(1+f)*V_exit - mdot_air*V_cruise;
            %     % disp(thrust_calc - thrust);
            % 
            %     % Update T04 guess
            %     if thrust_calc > thrust
            %         T04_high = T04_guess;
            %     else
            %         T04_low = T04_guess;
            %     end
            % 
            % end