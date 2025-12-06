function [D, dDOverdM, G, dGOverdM, N, dNOverdM] = MachNumberFunctions(Mach, gamma)

%First find the Mach number function D
D = Mach./((1+((gamma - 1)./2).*Mach.^2).^((gamma+1)./(2*(gamma-1))));
dDOverdM = (D./Mach).*(1-Mach.^2)./(1+((gamma-1)./2)*Mach.^2);

%Then find the mach number function G and its derivative
G = (1 + gamma * (Mach.^2))./((1+((gamma - 1)./2)*(Mach.^2)).^(gamma./(gamma-1))); %Find G
dGOverdM = gamma.*G.*Mach.*((2./(1 + gamma.*(Mach.^2))) - (1./(1 + ((gamma - 1)./2).*(Mach.^2)))); %dG/dM

%Then find the mach number function N and its derivative
N = D./G; %This is why we solved for D at the beginning
dNOverdM = N.*((1./Mach) + (((gamma - 1)./2).*(Mach./(1 + ((gamma - 1)./2).*(Mach.^2)))) - ((2.*gamma.*Mach)./(1 + gamma.*(Mach.^2)))); %Find dN/dM