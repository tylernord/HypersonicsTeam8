%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%           Project            %%%%%%%%%%%%
%%%%%%%%%%%%        Inlet/Isolator        %%%%%%%%%%%%
%%%%%%%%%%%%        Gage Bachmann         %%%%%%%%%%%%
%%%%%%%%%%%%           AAE 537            %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all; fclose all;
%% Starting and Unstart Mach Number
gamma = 1.4; %Assumed
IsoHeights = linspace(1, 10, 100); %Array of isolator heights
depth = 10; %inches
H1 = linspace(5, 15, 5); %Inlet is 10 inches
A1 = H1.*depth;
A2 = IsoHeights.*depth;
M1_init = 0.5;
M2_init = 1.5;
figure (7)
legendEntries = {};   % initialize as empty cell array
for j = 1:length(H1)
    for i = 1:length(A2)
        if i == 99
            fdsfdfs= 2;
        end
        [M1(j, i)] = NewtonsMethodDSolvingM1(1, M1_init, gamma, A1(j), A2(i));
        if isnan(M1)
            Mcrit(j, i) = NaN;
        else
            [Mcrit(j, i)] = NewtonsMethodForFindM2(M1(j, i), M2_init, gamma);
        end
       
    end
     h(j) = plot(Mcrit(j, :), IsoHeights(:));
     ylabel("Isolator Heights [in]")
     xlabel("Critical Mach Number [-]")
     hold on
     legendEntries{end+1} = sprintf("Inlet Height %0.2f", H1(j));
end
legend(h, legendEntries);







%% Inlet Parametric Curves
%Specifiy Ranges
%ArrayNum = 10;
M0 = linspace(3.5, 4.5, 100); %linspace(4, 5.5, ArrayNum); %Specify Inlet Mach Array
%theta = linspace(7, 15, ArrayNum); %Specify theta [can change]
AltRange = linspace(20000, 35000, 100); %Specify alt range [m]
theta = zeros(1, length(M0));



%Specify Target Values

%T_ratio = 6; %Specify as constant [can change later]

%Constants
gamma = 1.4; %Assumed
multiplier = 10;
TotalLength = 3.6*multiplier;
NormalizedHeight = 1*multiplier;
hd = 0.21*multiplier; % [can change this value]


%Start For Loop 
%for i = 1:length(theta)
%j = 1; %Initialize j

error = inf;
for i = 1:length(M0)
    thetaguess = 7;
    while abs(error) > 0.2
        theta(i) = thetaguess;
        %[T0,~,P0,~,~,~] = atmosisa(AltRange, extended=true); %Get Static Pressure and Temp in Pa and K
        Beta1(i) = FindBeta(thetaguess, gamma, M0(i));
        width(i) = NormalizedHeight./tand(Beta1(i));
        Cx(i) = (NormalizedHeight - hd)./tand(theta(i)) - width(i);
        Cx_norm(i) = Cx(i)./(Cx(i) + width(i));
        width_norm(i) = width(i)./(Cx(i) + width(i));
        l1(i) = sqrt(NormalizedHeight.^2 + width(i).^2);
        l2(i) = sqrt(hd.^2 + Cx(i).^2);
        
        %Get Beta2 and M2 from front
        M1N(i) = M0(i).*sind(Beta1(i));
        M2N(i) = NewtonsMethodForFindM2(M1N(i), 0.5, gamma);
        M2(i) = M2N(i)./sind(Beta1(i) - thetaguess);
        %Beta2_forward(i) = FindBeta(thetaguess, gamma, M2(i));
        Beta2_forward(i) = SolveBetaWeak(thetaguess, M2(i));

        %Find Beta from constant duct
        real(i) = asind(sind(Beta1(i) - thetaguess).*(l1(i)./l2(i)));
        if ~isreal(real(i))
            Beta2(i) = NaN;
        else
            Beta2(i) = real(i);
        end
        error = Beta2_forward(i) - Beta2(i);
        thetaguess = thetaguess + error.*0.01;
    end
    theta(i) = thetaguess - error.*0.01;
    M2(i) = FindMach(theta(i), Beta2(i));
    M2N1(i) = M2(i).*sind(Beta2(i));
    M3N(i) = NewtonsMethodForFindM2(M2N1(i), 0.5, gamma);
    M3(i) = M3N(i)./sind(Beta2(i) - theta(i));
    Tratio(i) = (1 + (gamma-1)./2.*M0(i).^2)./(1 + (gamma - 1)./2.*M3(i).^2);

    %Get Pressure Recovery
    [~, ~, G1(i), ~, ~, ~] = MachNumberFunctions(M0(i), gamma);
    [~, ~, G3(i), ~, ~, ~] = MachNumberFunctions(M3(i), gamma);
    [~, ~, G1N(i), ~, ~, ~] = MachNumberFunctions(M1N(i), gamma);
    [~, ~, G3N(i), ~, ~, ~] = MachNumberFunctions(M3N(i), gamma);
    Pratio(i) = (G1(i)./G3(i)).*((1 + (gamma - 1)./2.*M0(i).^2)./(1 + (gamma - 1)./2.*M3(i).^2)).^(gamma./(gamma - 1));
    PressureRecovery(i) = (G1N(i)./G3N(i));
    [P2_P1(i)] = ObliqueShockPRatio(M0(i), Beta1(i), gamma);
    [P3_P2(i)] = ObliqueShockPRatio(M2(i), Beta2(i), gamma);
    P3_P1(i) = P2_P1(i).*P3_P2(i);

    %Isolator:
    Re_theta = 20000; %Reynolds Number
    theta_H = 0.04; %Another parameter
    P2_P1max = 1 + (2.*gamma)./(gamma + 1).*(M3.^2 - 1); %Max Pressure Rise

    IsoL(i) = TotalLength - (NormalizedHeight - hd)./tand(theta(i));
    error = inf;
    P2_P1guess = 1;
    while abs(error) > 1e-4
        LoverD = sqrt(theta_H)./((M3(i).^2 - 1).*(Re_theta.^0.25))...
            .*(50.*((P2_P1guess) - 1) + 170.*((P2_P1guess)- 1).^2);
        error = LoverD - (IsoL(i)./hd);
        P2_P1guess = P2_P1guess - error.*0.01;
    end
    P2_P1_new(i) = P2_P1guess + error.*0.01;
    P2_P1_percent_new(i) = P2_P1_new(i)./P2_P1max(i);
    C1_new(i) = (1./(gamma.*M3(i))).*((1 + gamma.*M3(i).^2) - P2_P1_new(i)).*(1 + (gamma - 1)./2.*M3(i).^2).^(-0.5);
    IsoExitMach_new(i) = C1_new(i)./(sqrt(1 - (gamma - 1)./2.*C1_new(i).^2));
    P0ratio(i) = P2_P1_new(i).*((1 + (gamma - 1)./2.*IsoExitMach_new(i).^2)./(1 + (gamma - 1)./2.*M3(i).^2)).^(gamma./(gamma - 1));
    TotalPressRec(i) = PressureRecovery(i).*P0ratio(i);
    CombustorInletM(i) = NewtonsMethodForFindM2(IsoExitMach_new(i), 0.5, gamma);
    [~, ~, G4(i), ~, ~, ~] = MachNumberFunctions(IsoExitMach_new(i), gamma);
    [~, ~, G5(i), ~, ~, ~] = MachNumberFunctions(CombustorInletM(i), gamma);
    NSRec(i) = (G4(i)./G5(i));
    TotalPressRec(i) = PressureRecovery(i).*P0ratio(i).*NSRec(i);



    error = inf;
   
    TipSurfaceX = [0, (NormalizedHeight - hd)./tand(theta(i))];
    TipSurfaceY = [NormalizedHeight, hd];
    CxX = [width(i), TotalLength];
    CxY = [0, 0];
    InletX = [(NormalizedHeight - hd)./tand(theta(i)), TotalLength];
    InletY = [hd, hd];
    TopX = [0, TotalLength];
    TopY = [NormalizedHeight, NormalizedHeight];
    FirstShockX = [0, width(i)];
    FirstShockY = [NormalizedHeight, 0];
    SecondShockX = [width(i), (NormalizedHeight - hd)./tand(theta(i))];
    SecondShockY = [0, hd];
    adjust = 0;

    figure(1)
    plot(TipSurfaceX+adjust, TipSurfaceY+adjust, 'k')
    hold on 
    plot(CxX+adjust, CxY+adjust, 'k');
    hold on
    plot(InletX+adjust, InletY+adjust, 'g')
    hold on
    plot(TopX+adjust, TopY+adjust, 'k')
    hold on
    plot(FirstShockX+adjust, FirstShockY+adjust, '--r')
    hold on
    plot(SecondShockX+adjust, SecondShockY+adjust, '--r')
    


    xlim([0, TotalLength + 2])
    ylim([0, 1.5])
    xlabel('Inlet and Isolator Length [in]')
    ylabel('Inlet and Isolator Height [in]')
    title(['Two-Turn Inlet Geometry For M0 = ' num2str(M0(i)), ' where Theta = ', num2str(theta(i))]);
    axis equal
    %pause(0.5) %Pause for 5 seconds
    hold off
end

figure(2)
plot(M0(:), Tratio(:))
xlabel("Flight Mach Number [-]")
ylabel("Static Temperature Ratio T3/T1 [-]")

figure(3)
plot(M0(:), Pratio(:))
hold on
plot(M0(:), P3_P1(:))
hold on
plot(M0(:), PressureRecovery(:))
xlabel("Flight Mach Number [-]")
ylabel("Ratios [-]")
legend("Static Press Ratio", "P0 Ratio NASA Eqn", "P0 Ratio GN/GN")

figure(4)
plot(M0, CombustorInletM)
xlabel("Flight Mach Number")
ylabel("Combustor Inlet Mach Number [Assuming NS]")

figure(5)
plot(M0, P2_P1_percent_new)
xlabel("Flight Mach Number")
ylabel("% of Max Isolator P2/P1")

figure(6)
plot(M0, TotalPressRec)
xlabel("Flight Mach Number")
ylabel("Total Pressure Recovery of Inlet and Isolator")

%% MIL-E-5007D Inlet Performance Range

%Uncomment if you want to see the plot

% M0 = linspace(0, 6, 100);
% for i = 1:length(M0)
%     if M0(i) < 1
%         P2_P1(i) = 1;
%     elseif M0(i) >=1 && M0(i) < 5
%         P2_P1(i) = 1 - 0.075.*(M0(i) - 1).^1.35;
%     else
%         P2_P1(i) = 800./(M0(i).^4 + 925);
%     end
% end
% figure(5)
% plot(M0, P2_P1)


