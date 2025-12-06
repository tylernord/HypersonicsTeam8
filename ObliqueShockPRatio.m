function [Pratio] = ObliqueShockPRatio(M1, beta, gamma)

Pratio = (((gamma + 1).*M1.^2.*sind(beta).^2)./((gamma - 1).*M1.^2.*sind(beta).^2 + 2)).^(gamma./(gamma - 1))...
    .*((gamma + 1)./(2.*gamma.*M1.^2.*sind(beta).^2 - (gamma - 1))).^(1./(gamma - 1));