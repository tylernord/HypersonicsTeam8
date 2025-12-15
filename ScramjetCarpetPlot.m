%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%           Project            %%%%%%%%%%%%
%%%%%%%%%%%%     Combustor/Expansion      %%%%%%%%%%%%
%%%%%%%%%%%%        Gage Bachmann         %%%%%%%%%%%%
%%%%%%%%%%%%           AAE 537            %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all; fclose all;
%% Combustor/Expansion Parametric Curves
%Givens/Constants
gamma = 1.4;
R = 287; 
%f = 0.005; %Assume very low
% A3_A2 = 4; %Expanison Area Ratio. Constant [can change input]
% 
% AltRange = 20000:1000:35000; %Specify alt range [m]
% M1range = 0.5:0.5:2; %Combustor Inlet Mach [-]
% A2_A1 = 1:1:5; %Combustor Area Ratio [-]
% FlightMach = 4;
% 
% for i = 1:length(M1range) %Get Plots Per M1
%     for j = 1:length(A2_A1)
%         for k = 1:length(AltRange)
%             [T1(k),~,P3(k),~,~,~] = atmosisa(AltRange(k), extended=true); %Get Static Pressure and Temp in Pa and K
%             M2(i, j) = sqrt((M1range(i).^2)./A2_A1(j));
%             %M3(i, j) = NewtonsMethodDFindingM2(M2(i, j), 1.2, gamma, 1, A3_A2); 
%             P1(j, k, i) = P3(k).*(((1 + (gamma - 1)./2.*M3(i, j).^2)./...
%                 (1 + (gamma - 1)./2.*M2(i, j).^2)).^(gamma./(gamma - 1)));
%             T02_T01(j, k, i) = ((M1range(i)./M2(i, j)).^2).*((1 + (gamma - 1)./2.*M2(i, j).^2)./...
%                 (1 + (gamma - 1)./2.*M1range(i).^2));
% 
%             T01(j, k, i) = T1(k).*(1 + (gamma - 1)./2.*FlightMach.^2);
%             T02(j, k, i) = T01(j, k, i).*T02_T01(j, k, i); %constant through inlet and isolator
%             T3(j, k, i) = T02(j, k, i)./(1 + (gamma - 1)./2.*M3(i, j).^2);
%             u1(j, k, i) = M1range(i).*sqrt(gamma.*R.*T1(k));
%             u3(j, k, i) = M3(i, j).*sqrt(gamma.*R.*T3(j, k, i));
% 
%             ST(j, k, i) = (1 + f).*u3(j, k, i) - u1(j, k, i);
%         end
%     end
%     figure(i)
%     subplot(3, 1, 1)
%     plot(AltRange(:), P1(:, :, i))
%     %xlabel("Altitude [m]")
%     ylabel("P1 [Pa]")
%     title(['P1 For M1 = ' num2str(M1range(i))]);
%     legend("A2/A1 = 1", "A2/A1 = 2", "A2/A1 = 3", "A2/A1 = 4", "A2/A1 = 5")
%     subplot(3, 1, 2)
%     plot(AltRange(:), T02_T01(:, :, i))
%     xlabel("Altitude [m]")
%     ylabel("T02/T01 [-]")
%     title(['T02/T01 For M1 = ' num2str(M1range(i))]);
%     subplot(3, 1, 3)
%     plot(AltRange(:), ST(:, :, i))
%     xlabel("Altitude [m]")
%     ylabel("Specific Thrust [N/(kg/s)]")
%     title(['Specific Thrust For M1 = ' num2str(M1range(i))]);
%     hold off
% 
% end

%% Carpet Plot for Scramjet
M0 = 4:1:6;
alt = 35000; %Altitude at 25km
A1 = 0.064516; %M^2
g = 9.8; 

%Find the corresponding altitudes and mach numbers for what you need
T04 = [1600, 1800, 2000]; %Combustor exit temperatures [K]
PiC = [0.3, 0.4, 0.5]; %Compressor Ratio Values [-]
gamma1 = 7./5; 
MW = 28.97./1000; %Molecular Weight
R = 8.3145; %Universal Gas Constant
Cp1 = (gamma1.*(R./MW))./(gamma1 - 1);
nc = 0.85; %Compressor efficiency [-]
nb = 0.98; %Combustor Efficiency [-]
%nt = 0.91; %Turbine Efficiency [-]
nn = 0.96; %Nozzle Efficiency [-]
ni = 0.1;

L_D = [9, 8.5, 8];
%n0 = 65; %Overall Efficiency
massratio = 6; 
LHV = 43e6; %Lower heating value of JP10


%Isolator:
Re_theta = 20000; %Reynolds Number
theta_H = 0.04; %Another parameter

%Start Loop
for i = 1:length(M0)
    ST = zeros(length(PiC), length(T04));
    SFC = zeros(length(PiC), length(T04));
    OF = zeros(length(PiC), length(T04));
    for j = 1:length(PiC)
        for k = 1:length(T04)
            %Station 1: Before Inlet
            nd = 1 - (0.00689).*M0(i).^2.516;
            [T0,a0,P0,rho0,~,~] = atmosisa(alt, extended=true); %Temp [K] and Pressure [Pa]
            P01(j, k) = P0.*(1 + nd.*(gamma1 - 1)./2.*M0(i).^2).^(gamma1./(gamma1 - 1));
            T01(j, k) = T0.*(1 + (gamma1 - 1)./2.*M0(i).^2);
            u0 = M0(i).*a0;
    
            %Station 2: After Inlet
            P02(j, k) = PiC(j).*P01(j, k);
            T02(j, k) = T01(j, k).*(1 + (1./nc).*(((PiC(j)).^((gamma1 - 1)./gamma1)) - 1));
            %Get M2
            [~, ~, G0, ~, ~, ~] = MachNumberFunctions(M0(i), gamma1);
            G2 = G0./PiC(j);
            M2guess = linspace(0.1, 5, 100);
            for m = 1:length(M2guess)
                [~, ~, G2guess, ~, ~, ~] = MachNumberFunctions(M2guess(m), gamma1);
                error = abs(G2guess - G2);
                if abs(G2guess - G2) < 1e-2
                    M2(j, k) = M2guess(m);
                    break
                end
            end

            %Station 3
            %Assuming a 75% static pressure rise
            P3_P2(j, k) = ni.*(1 + (2.*gamma1)./(gamma1 + 1).*(M2(j, k).^2 - 1));
            LoverD(j, k) = sqrt(theta_H)./((M2(j, k).^2 - 1).*(Re_theta.^0.25))...
                .*(50.*(P3_P2(j, k) - 1) + 170.*(P3_P2(j, k) - 1).^2);
            C1(j, k) = (1./(gamma1.*M2(j, k))).*((1 + gamma1.*M2(j, k).^2) - P3_P2(j, k)).*(1 + (gamma1 - 1)./2.*M2(j, k).^2).^(-0.5);
            M3(j, k) = C1(j, k)./(sqrt(1 - (gamma1 - 1)./2.*C1(j, k).^2));
            T03(j, k) = T02(j, k); %Stagnation temperature is constant accross isolator
            P03(j, k) = P02(j, k).*P3_P2(j, k).*((1 + (gamma1 - 1)./2.*M3(j, k).^2)./(1 + (gamma1 - 1)./2.*M2(j, k).^2)).^(gamma1./(gamma1 - 1));
            

            %Station 4
            if T03(j, k) > T04(k)
                ST(j, k) = NaN; %Don't plot these
                SFC(j, k) = NaN; %Don't plot these
                OF(j, k) = NaN;
            else
    
                % % Guess O/F Ratio
                
                OF(j, k, i) = 20; %Esitmated first OF Ratio Guess
                Te = 0; %Initialize the exit temp
                error = inf;
                iter = 1;
                
                while abs(error) > 5
                   
                    if iter == 1
                        %fprintf("Starting Iteration For CEA\n");
                    elseif mod(iter, 1) == 0
                        %fprintf("CEA iteration %i] Estimated Error = %0.4f \n", iter, error);
                    end
    
                    [gamma2, Cp2, Te, MW] = HW3CEA(P03(j, k), T03(j, k), OF(j, k, i));
                     
               
                     error = T04(k) - Te;
                
                     OF(j, k, i) = OF(j, k, i) - error*0.03;
                     % if iter > 40
                     %     break
                     % end
                     iter = iter + 1;
                     fclose all;
                    delete("Project1.inp")
                    delete("Project1.out")
                end
                R2(j, k) = 8314./MW;
                OF(j, k, i) = OF(j, k, i) + error*0.03;
                %fprintf("O/F Ratio for T04 %0.f and PiC %0.f is %0.4f \n", T04(j), PiC(k), OF{i}(k, j));
                f_cea(j, k, i) = 1./OF(j, k, i);
                f_real(j, k, i) = f_cea(j, k, i)./nb; %Get the real fuel to air ratio
                M4(j, k) = sqrt((M3(j, k).^2)./((T04(k)./T03(j, k)) + ((gamma2 - 1)./2.*M3(j, k).^2).*((T04(k)./T03(j, k)) - 1)));
                A4_A3(j, k) = (M3(j, k).^2)./(M4(j, k).^2);
                P04(j, k) = P03(j, k).*((1 + (gamma2 - 1)./2.*M4(j, k).^2)./(1 + (gamma - 1)./2.*M3(j, k).^2)).^(gamma2./(gamma2 - 1)); %Assuming constant pressure combustor
                P4(j, k) = P04(j, k)./((1 + (gamma2 - 1)./2.*M4(j, k).^2).^(gamma2./(gamma2 - 1)));
                T4(j, k) = T04(k).*(1 - nn.*(1 - (P0./P04(j, k)).^((gamma2 - 1)./gamma2)));
                a9(j, k) = sqrt(gamma2.*R2(j, k).*T4(j, k));
                u4(j, k) = a9(j, k).*M4(j, k);
              
    
                %Exit
                u9(j, k) = sqrt(2.*Cp2.*T04(k).*nn.*(1 - ((P0./P04(j, k)).^((gamma2 - 1)./gamma2))));
                ST(j, k, i) = (1 + f_real(j, k, i)).*u4(j, k) - u0; %Specific Thrust [Ns/kg]
                SFC(j, k, i) = (f_real(j, k, i)./ST(j, k, i)).*3600; %Specific Fuel Consumption [kg/Nh]
                if SFC(j, k, i) < 0
                    SFC(j, k, i) = NaN;
                end
                if ST(j, k, i) < 0
                    ST(j, k, i) = NaN;
                end
            
                %Get Thrust
                T(j, k, i) = ST(j, k, i).*(rho0.*u0.*A1);
                mdotf(j, k, i) = (rho0.*u0.*A1).*f_real(j, k, i);
                Isp(j, k, i) = T(j, k, i)./(mdotf(j, k, i).*g);
                %Find range equaiton
                
             
                v(j, k, i) = (u0./u4(j, k));
                np(j, k, i) = (2.*v(j, k, i).*(f_real(j, k, i) + 1 - v(j, k, i)))./...
                    (2.*v(j, k, i).*(f_real(j, k, i) + 1 - v(j, k, i)) + (f_real(j, k, i) + 1).*(1 - v(j, k, i)).^2);
                nt(j, k, i) = (f_real(j, k, i).*(0.5 + (v(j, k, i).^2)./2) + 0.5 - (v(j, k, i).^2)./2)./...
                    (f_real(j, k, i).*v(j, k, i).^2.*(LHV./(u0.^2)));
                %n0(j, k, i) = np(j, k, i).*(nt(j, k, i));
                n0(j, k, i) = u0./(SFC(j, k, i).*LHV);
                Range(j, k, i) = n0(j, k, i).*LHV.*L_D(i).*(1./g).*log(massratio);
            end
        end
        
    end

    figure(i)
    hold on;
    for n = 1:size(ST, 1)
        plot(ST(n, :, i), SFC(n, :, i), 'b-o');  % Blue lines across rows PIc
        hold on
    end
    for m = 1:size(ST, 2)
        plot(ST(:, m, i), SFC(:, m, i), 'r--o');  % Red lines down columns TO4
        hold on
    end
    xlabel('Specific Thrust [Ns/kg]')
    ylabel('Specific Fuel Consumption [kg/hN]')
    title(["M0 = ", M0(i), " and Alt = ", alt, "km"])
end