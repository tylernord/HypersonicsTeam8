%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%        Homework 3 P2         %%%%%%%%%%%%
%%%%%%%%%%%%        Gage Bachmann         %%%%%%%%%%%%
%%%%%%%%%%%%           AAE 537            %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [RealGamma, Cp, Te, MW] = HW3CEA(P03, T03, OF)
    %fclose all;
    addpath('C:\Users\Gage Bachmann\Downloads\cea_rocket_run-main\cea_rocket_run-main\');
    
    inp = containers.Map;
    inp('type') = 'eq';              % Sets the type of CEA calculation
    inp('p') = P03./100000;                      % Chamber pressure Pa to Bar
    inp('p_unit') = 'bar';              % Chamber pressure units
    inp('o/f') = OF;                % Mixture ratio (Check if this is correct)
    %inp('pip') = 1;                    % Supersonic area ratios (Check if this is correct)
    inp('fuel') = 'JP-10(L)';
    inp('ox') = 'Air';
    inp('fuel_wt%') = 100;
    inp('ox_wt%') = 100;
    inp('fuel_t') = 298.15; % Assuming room temp fuel inlet temp
    inp('ox_t') = T03;
    inp('fuel_t_unit') = 'K';
    inp('ox_t_unit') = 'K';
    inp('file_name') = 'Project.inp';    % Input/output file name
    inp('verbose') = 0;
    
    data = cea_rocket_run(inp);     % Call the CEA MATLAB code
    data_eq = data('eq');
    
    %Call Values
    gammas = squeeze(data_eq('gammas')); %Get Gamma
    L = squeeze(data_eq('(dlv/dlp)t')); %Get whatever the fuck this is
    RealGamma = -gammas.*L; %Get the actual gamma value
    Cp = squeeze(data_eq('cp')); %Get Gamma
    Te = squeeze(data_eq('t'));
    MW = squeeze(data_eq('m'));





    % cstar_p2 = squeeze(data_eq('cstar')); %Get C*
    % cstar_p2 = cstar_p2(2); %Just get one
    % cf_p2 = data_eq('cf'); %Get Cf
    % cf_p2 = cf_p2(2); %Get the one at the throat
    % Pressures = squeeze(data_eq('p')); %Get pressures
    % Pc_p2 = Pressures(1); %Get chamber Pressure
    % Pe_p2 = Pressures(2);
    % Isp_p2 = squeeze(data_eq('isp'));
    % Isp_p2 = Isp_p2(2)./9.81; %Get Isp in Seconds
    % gamma_p2 = squeeze(data_eq('gammas'));
    % gammae_p2 = gamma_p2(2);
    % gamma_p2 = (gamma_p2(1) + gamma_p2(2))./2; %Take an average of the gamma values
    % Me_p2 = squeeze(data_eq('mach'));
    % Me_p2 = Me_p2(2); %Get exit Mach Number
    % Te_p2 = squeeze(data_eq('t'));
    % Te_p2 = Te_p2(2); %Get exit Temperature
    % MW = squeeze(data_eq('mw'));
    % MWe = MW(2)./1000; %Get Molecular Weight at the Exit in kmol
    % ivac = squeeze(data_eq('ivac'));
    % ivac = ivac(2);
end