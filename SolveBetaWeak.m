function beta_weak_deg = SolveBetaWeak(theta_deg, M)
% SolveBetaWeak  Computes the weak-shock solution for β (degrees)
% from the θ–β–M oblique shock relation.
%
% INPUTS:
%   theta_deg : flow deflection angle θ (deg)
%   M         : upstream Mach number
%
% OUTPUT:
%   beta_weak_deg : weak shock angle β (deg), or NaN if no solution.

    gamma = 1.4;
    theta = deg2rad(theta_deg);

    % θ–β–M residual function
    f = @(beta) tan(theta) - ...
        ( 2.*cot(beta) .* (M.^2 .* sin(beta).^2 - 1) ./ ...
        (M.^2 .* (gamma + cos(2.*beta)) + 2) );

    % Beta must lie between Mach angle and 90°
    beta_min = asin(1/M) + 1e-6;   
    beta_max = pi/2 - 1e-6;        

    % Scan to find weak-shock sign change (first root)
    Nscan = 2000;
    bvec = linspace(beta_min, beta_max, Nscan);
    fvec = f(bvec);

    idx = find(fvec(1:end-1).*fvec(2:end) < 0, 1, 'first');

    if isempty(idx)
        beta_weak_deg = NaN;  % No solution → detached shock region
        return;
    end

    % Solve for β inside this bracket
    beta_weak = fzero(f, [bvec(idx), bvec(idx+1)]);

    beta_weak_deg = rad2deg(beta_weak);

end