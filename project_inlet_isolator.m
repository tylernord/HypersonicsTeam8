% AAE 537 Hypersonic Propulsion
% Final Project
% constant width two turn inlet with isolator 
% Tyler Nord
clear; clc;

%% Variable Instantiations
    % Inputs
H = 25*1000;                                                               %cruise altitude, m
M0 = 4.5;                                                                  %cruise Mach number
% target_Po3 = 250000;                                                       %isolator exit stag. pressure, Pa
% A0 = ;                                                                   %capture area, m^2
Isolator_L_margin = 1.25;                                                  %ratio of isolator length to shock train length at cruise
    % Assumed constants
gamma = 1.4;
theta_H = 0.04; %boundary layer momentum thickness to duct height ratio 
Re = 20000; 
    % Initial calcs
[T0, a0, P0, rho0, ~, ~] = atmosisa(H, 'extended', true);                  %cruise static conditions, [K, m/s, Pa, kg/m^3, ~, ~]
Po0 = Po(P0,M0,gamma);                                                     %freestream stagnation pressure, Pa
% target = target_Po3/Po0;                                                   %stagnation pressure ratio
targets = [0.35, 0.4, 0.45, 0.5];                                                              % [[ Placeholder ]]
    % Values used by solver
N_theta = 500;
N_iso = 500;
thetas_deg = linspace(10,40,N_theta);                     % wedge/turning angle sweep [deg]
thetas = thetas_deg*pi/180;                               % in radians
M1 = M0;
    % Preallocate (NaNs to mask invalid cases)
beta1         = NaN(size(thetas));
beta2         = NaN(size(thetas));
M2            = NaN(size(thetas));
M3            = NaN(size(thetas));
T3_T0         = NaN(size(thetas));
Po3_Po0       = NaN(size(thetas));
P3            = NaN(size(thetas));
Po4_Po3_req   = NaN(size(thetas));
Po4_Po3_min   = NaN(size(thetas));
Po4_Po3_max   = NaN(size(thetas));
CompEff       = NaN(size(thetas));
P4_P3_max     = NaN(size(thetas));
P4_P3         = NaN(length(thetas), N_iso);
M4            = NaN(length(thetas), N_iso);
P4            = NaN(length(thetas), N_iso);
Po3           = NaN(length(thetas), N_iso);
Po4           = NaN(length(thetas), N_iso);
Po4_Po3       = NaN(length(thetas), N_iso);
Po4_Po0       = NaN(length(thetas), N_iso);
k             = NaN(length(targets),length(thetas));
L_H           = NaN(length(targets),length(thetas));
M4_opt        = NaN(length(targets),length(thetas));
P4_opt        = NaN(length(targets),length(thetas));
Po3_opt       = NaN(length(targets),length(thetas));
Po4_opt       = NaN(length(targets),length(thetas));
Po4_Po3_opt   = NaN(length(targets),length(thetas));
Po4_Po0_opt   = NaN(length(targets),length(thetas));


%% Solver
for i = 1:length(targets)
    target = targets(i);
    for j = 1:length(thetas)
        theta = thetas(j);
    
    % ---------- Shock 1 (oblique) ----------
        b1 = ShockAngle(theta, M1, gamma);
        if ~isfinite(b1), continue; end
        beta1(j) = b1;
    
        M1n = M1*sin(b1);                              % upstream normal Mach to shock 1
        [T2_T1, p2_p1] = ShockRatios(M1n, gamma);      % canonical normal-shock T & p ratios
    
        % Reconstruct full M2 from downstream normal component
        M2n   = ShockM2n(M1n, gamma);                  % downstream normal Mach
        M2(j) = M2n / sin(b1 - theta);                 % full downstream Mach after shock 1
        if ~isfinite(M2(j)) || M2(j) <= 1, continue; end
    
        % Stagnation pressure ratio across shock 1 (use FULL M1 & M2)
        Po2_Po1 = PoRatioFromFullM(p2_p1, M1, M2(j), gamma);
    
    % ---------- Shock 2 (oblique) ----------
        b2 = ShockAngle(theta, M2(j), gamma);
        if ~isfinite(b2), continue; end
        beta2(j) = b2;
    
        M2n_p = M2(j)*sin(b2);                         % upstream normal Mach to shock 2
        [T3_T2, p3_p2] = ShockRatios(M2n_p, gamma);
    
        % Static temperature ratio relative to freestream static:
        T3_T0(j) = T2_T1 * T3_T2;
    
        % Reconstruct full M3 and stagnation ratio across shock 2
        M3n   = ShockM2n(M2n_p, gamma);
        M3(j) = M3n / sin(b2 - theta);
        if ~isfinite(M3(j)), continue; end
        if M3(j) < 1
            fprintf('\nInlet brings flow subsonic at or above Theta = %.0f degrees\n', thetas_deg(j))
            j_break = j;
            break;
        end
    
        Po3_Po2 = PoRatioFromFullM(p3_p2, M2(j), M3(j), gamma);
        % Accumulate total-loss across both shocks
        Po3_Po0(j) = Po2_Po1 * Po3_Po2;
        P3(j) = p3_p2*p2_p1*P0;
        CompEff(j) = (T3_T0(j)-(1/Po3_Po0(j))^((gamma-1)/gamma))/(T3_T0(j)-1);
      
    % ------------- Isolator -------------
        % required stagnation pressure ratio across isolator
        Po4_Po3_req(i,j) = target / Po3_Po0(j);
        
        % isolator Po limits
        Po4_Po3_min(j) = PoShock(M3(j), gamma);   % normal-shock case (most loss)
        Po4_Po3_max(j) = 1.0;                     % isentropic (no loss) upper bound
        
        % static pressure ratio limit for normal shock (same as before)
        P4_P3_max(j)   = PShock(M3(j), gamma);
        
        % feasibility check: can the isolator actually hit this target?
        if Po4_Po3_req(i,j) < Po4_Po3_min(j) || Po4_Po3_req(i,j) > Po4_Po3_max(j)
            fprintf(['\nTarget Po4/Po0 = %.3f is infeasible at Theta = %.2f deg ', ...
                     '(required Po4/Po3 = %.3f, allowed [%.3f, %.3f])\n'], ...
                     target, thetas_deg(j), Po4_Po3_req(i,j), Po4_Po3_min(j), Po4_Po3_max(j));
        
            % mark this theta as infeasible and skip isolator sizing
            k(i,j)           = NaN;
            M4_opt(i,j)      = NaN;
            P4_opt(i,j)      = NaN;
            Po4_opt(i,j)     = NaN;
            Po4_Po3_opt(i,j) = NaN;
            Po4_Po0_opt(i,j) = NaN;
            L_H(i,j)         = NaN;
            continue;
        end
    
        %possible pressure ratios from a shock train (this is my variable for the design space)
        P4_P3(j,:) = linspace(1.01, P4_P3_max(j), N_iso);
        %compute the mach number at the exit
        gm1 = gamma-1;
            %use formula from Lec 11 slide 23
        f = @(M4,P4_P3) (1./P4_P3).*(M3(j)./M4).*sqrt((1+gm1/2*M3(j).^2)./(1+gm1/2*M4.^2)) ...
                        - (((1+gamma*M3(j).^2)./P4_P3)-1)./(gamma*M4.^2);
        M4(j,:) = arrayfun(@(p) fsolve(@(M4) f(M4,p), 0.2, ...
                         optimoptions('fsolve','Display','off')), P4_P3(j,:));
        %compute stagnation pressure loss across isolator
        P4(j,:) = P3(j).*P4_P3(j,:);
        Po3(j,:) = Po(P3(j),M3(j),gamma);
        Po4(j,:) = Po(P4(j,:),M4(j,:),gamma);
        Po4_Po3(j,:) = Po4(j,:)./Po3(j,:);
        Po4_Po0(j,:) = Po4_Po3(j,:).*Po3_Po0(j);
        %find index k such that Po4_Po3(j,k) = Po4_Po3_req(j)
        [~, k(i,j)] = min(abs(Po4_Po3(j,:) - Po4_Po3_req(i,j)));   % closest match
        %compute shock train length via correlation in Waltrup & Billing
        L_H(i,j) = Isolator_L_margin.*sqrt(theta_H).*(50.*(P4_P3(j,k(i,j))-1)+170.* ...
            (P4_P3(j,k(i,j))-1).^2)./((M3(j)^2-1).*Re^0.25);
        %collect post isolator values for optimal Po ratio case at each theta
        M4_opt(i,j) = M4(j,k(i,j));
        P4_opt(i,j) = P4(j,k(i,j));
        Po3_opt(i,j) = Po3(j,k(i,j));
        Po4_opt(i,j) = Po4(j,k(i,j));
        Po4_Po3_opt(i,j) = Po4_Po3(j,k(i,j));
        Po4_Po0_opt(i,j) = Po4_Po0(j,k(i,j));
    end
end

%% Plotting 
close all;

figure;
for n = 1:length(targets)
    plot(thetas*180/pi,Po4_Po0_opt(n,:),'-b','DisplayName','selected Po ratio'); hold on; 
    yline(targets(n),':k','DisplayName','target Po ratio'); 
end
grid on;
xlabel('Ramp Angle, Theta [deg]');
ylabel('Stag. Pressure Ratio Across Inlet+Isolator, Po4/Po0');
legend('location','best')

figure;
% --- Choose 4 colors that are clearly NOT red/black ---
clr = [0      0.4470 0.7410;  ... % blue
       0.4660 0.6740 0.1880;  ... % green
       0.3010 0.7450 0.9330;  ... % cyan
       0.4940 0.1840 0.5560]; ... % purple
plot(thetas*180/pi,Po4_Po3_min,'--r','DisplayName','min isolator Po recovery (normal shock)'); hold on;
yline(1,'-.r','DisplayName','max isolator Po recovery (isentropic)')
ax = gca;
ax.ColorOrder = clr;    % custom 4 colors
ax.ColorOrderIndex = 1; % restart the sequence before the loop
for n = 1:length(targets)
    plot(thetas*180/pi,Po4_Po3_req(n,:), ':k', 'HandleVisibility', 'off');
    plot(thetas*180/pi,Po4_Po3_opt(n,:), '-', 'DisplayName', compose("Po3/Po0 = %.3f", targets(n))); grid on;
end
title(compose("Isolator Design Possibilities\n2-turn inlet @ M_{0} = %.2f", M1))
xlabel('Ramp Angle, Theta [deg]');
ylabel('Stag. Pressure Ratio Across Isolator, Po3/Po2');
legend('location','best')
xlim([thetas_deg(1), thetas_deg(j_break)])

chosen_theta = 18; %deg
target_idx = 4; %Po3/Po0 selected as an index of the vector 'targets'
idx = find(abs(thetas_deg - chosen_theta) == min(abs(thetas_deg - chosen_theta)));
PlotInletGeometry(thetas(idx), beta1(idx), beta2(idx), M0, L_H(target_idx,idx))
TabulateResults(M1, M2(idx), M3(idx), M4_opt(target_idx,idx), Po3_Po0(idx), Po4_Po0_opt(target_idx,idx), Po4_opt(target_idx,idx)/1000, CompEff(idx), thetas_deg(idx), beta1(idx)*180/pi, beta2(idx)*180/pi)

%% Functions

function Po = Po(P,M,gamma)
    g=gamma;gm1=g-1;gp1=g+1;
    Po = P.*(1+gm1/2.*M.^2).^(g/gm1);
end

function [Po1_P1] = PoIsentropicRel(M, gamma)
    g = gamma; gm1 = g-1;
    Po1_P1 = (1+gm1*M^2/2)^(g/gm1);
end

function [T2_T1, p2_p1] = ShockRatios(M1n, gamma)
% Canonical normal-shock relations for a perfect gas
    g   = gamma; gp1 = g + 1; gm1 = g - 1;

    % Static pressure and density ratios
    p2_p1   = 1 + 2*g/gp1 * (M1n^2 - 1);
    rho2_r1 = (gp1*M1n^2) / (gm1*M1n^2 + 2);

    % Static temperature ratio
    T2_T1   = p2_p1 / rho2_r1;
end

function M2n = ShockM2n(M1n, gamma)
% Downstream **normal** Mach number after a normal shock.
% For oblique shocks, reconstruct full M2 via:  M2 = M2n / sin(beta - theta)
    g = gamma; gm1 = g - 1;
    den = g*M1n^2 - 0.5*gm1;
    if den <= 0, M2n = NaN; return, end
    M2n_sq = (1 + 0.5*gm1*M1n^2) / den;
    if M2n_sq < 0, M2n = NaN; return, end
    M2n = sqrt(M2n_sq);
end

function Po2_Po1 = PoRatioFromFullM(p2_p1, M1, M2, gamma)
% Stagnation-pressure ratio using FULL upstream/downstream Machs:
% Po2/Po1 = (p2/p1)*[(1 + (γ-1)/2 M2^2)^(γ/(γ-1))]/[(1 + (γ-1)/2 M1^2)^(γ/(γ-1))]
    g = gamma; gm1 = g - 1;
    Po2_Po1 = p2_p1 * ( (1 + 0.5*gm1*M2^2)^(g/gm1) / (1 + 0.5*gm1*M1^2)^(g/gm1) );
end

function [P2_P1] = PShock(M1, gamma)
    g = gamma; gm1 = g-1; gp1 = g+1; 
    Po2_Po1 = (gp1*M1^2/2/(1 + gm1*M1^2/2))^(g/gm1) * (1 + 2*g*(M1^2 - 1)/gp1)^(-1/gm1);
    den = g*M1^2 - 0.5*gm1;
    M2n = sqrt((1 + 0.5*gm1*M1^2) / den);
    Po1_P1 = PoIsentropicRel(M1, g);
    Po2_P2 = PoIsentropicRel(M2n, g);
    P2_P1 = Po2_Po1*Po1_P1/Po2_P2;
end

function [Po2_Po1] = PoShock(M1, gamma)
    g = gamma; gm1 = g-1; gp1 = g+1; 
    Po2_Po1 = (gp1*M1^2/2/(1 + gm1*M1^2/2))^(g/gm1) * (1 + 2*g*(M1^2 - 1)/gp1)^(-1/gm1);
end

function beta = ShockAngle(theta, M1, gamma)
% Weak-shock angle beta [rad] for given deflection theta [rad], upstream M1.
% Returns NaN if M1<=1 or theta exceeds theta_max (tiny tolerance).
% Uses residual minimization over (μ, π/2) with a relaxed tolerance.

    if M1 <= 1
        beta = NaN; return
    end

    [theta_max, ~] = ThetaMaxAttached(M1, gamma);
    if theta > theta_max + 1e-6
        beta = NaN; return
    end

    g = gamma;
    f = @(b) tan(theta) - 2*cot(b).*(M1^2.*sin(b).^2 - 1)./(M1^2.*(g + cos(2*b)) + 2);

    mu  = asin(1/M1);
    blo = mu + 1e-6;                 % avoid singularities at endpoints
    bhi = pi/2 - 1e-6;

    % Minimize residual over attached range (relaxed tolerance)
    obj = @(b) abs(f(b));
    [b_opt, res] = fminbnd(obj, blo, bhi);

    if ~isfinite(res) || res > 1e-3
        beta = NaN; return
    end

    beta = b_opt;
end

function [theta_max, beta_at_max] = ThetaMaxAttached(M, gamma)
% Max deflection for an attached oblique shock (weak branch), via
% maximizing theta(beta) on (mu, pi/2).
    g = gamma;
    mu = asin(1/M);
    eps = 1e-6;
    theta_of_beta = @(b) atan( 2*cot(b).*(M.^2.*sin(b).^2 - 1)./(M.^2.*(g + cos(2*b)) + 2) );
    [b_opt, neg_th] = fminbnd(@(b) -theta_of_beta(b), mu + 10*eps, pi/2 - 10*eps);
    theta_max   = -neg_th;
    beta_at_max = b_opt;
end

function [] = TabulateResults(M1, M2, M3, M4, Po3_Po0, Po4_Po0, Po4, CompEff, theta, beta1, beta2)
    fprintf('\nInlet)\n\tM0 = %.3f, M2 = %.3f, M3 = %.3f\n\t', M1, M2, M3)
    fprintf('theta = %.2f deg, beta1 = %.2f deg, beta2 = %.2f deg\n\t', ...
    theta, beta1, beta2)
    fprintf('Po3/Po1 = %.6f\n\t', Po3_Po0)
    fprintf('Inlet Compression Efficiency = %.5f\n', CompEff)
    fprintf('Isolator)\n\tM4 = %.3f\n\t', M4)
    fprintf('Po4/Po1 = %.6f, Po4 = %.3f kPa\n', Po4_Po0, Po4)
end

function [] = PlotInletGeometry(theta, beta1, beta2, M0, L_H)
    % Arbitrary first-shock length
    Lshock1 = 1/sin(beta1);
    
    % Starting point for first shock (top inlet)
    x0 = 0; 
    y0 = 1;
    
    % --- First shock (downward) ---
    x1 = x0 + Lshock1 * cos(beta1);
    y1 = y0 - Lshock1 * sin(beta1);
    
    num = 1 - x1*tan(theta) - y1;
    den = sin(beta2 - theta) + cos(beta2 - theta)*tan(theta);
    Lshock2 = num / den;
    
    % Coordinates of second shock end
    x2 = x1 + Lshock2 * cos(beta2 - theta);
    y2 = y1 + Lshock2 * sin(beta2 - theta);

    % Isolator geometry
    H = y2-y1;
    L = L_H*H;
    
    % --- Plot ---
    figure; hold on; axis equal; grid on
    xlabel('x'); ylabel('y')
    
    % Top inlet line
    x_top = [0, x2 + L];
    y_top = [1, 1];
    plot(x_top, y_top, 'k', 'LineWidth', 2)
    
    % First shock (down)
    plot([x0, x1], [y0, y1], '--b', 'LineWidth', 1.5)
    
    % Second shock (up)
    plot([x1, x2], [y1, y2], '--b', 'LineWidth', 1.5)
    
    % Cowl lip (from inlet top downward)
    x_cowl = linspace(0, x2, 100);
    y_cowl = 1 - x_cowl * tan(theta);
    plot(x_cowl, y_cowl, 'k', 'LineWidth', 2)
    
    % Inlet walls
    plot([x2, x2 + L], [y2, y2], 'k', 'LineWidth', 2)
    plot([x1, x2 + L], [y1, y1], 'k', 'LineWidth', 2)
    
    title(sprintf('Two-Turn Inlet Geometry with Isolator\n M_0 = %.2f, Theta = %.2f deg', M0, theta*180/pi))
    % ylim([-0.5, 1.5])
end
