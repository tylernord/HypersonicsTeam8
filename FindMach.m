function M = FindMach(theta_deg, beta_deg)
% Safe Mach solver. Returns NaN if no physical root exists.

    theta = deg2rad(theta_deg);
    beta  = deg2rad(beta_deg);
    gamma = 1.4;

    f = @(M) tan(theta) - ...
        (2*cot(beta) * (M.^2 .* sin(beta).^2 - 1) ./ ...
        (M.^2 * (gamma + cos(2*beta)) + 2));

    % Scan Mach space for sign changes:
    M_min = 1.0001;
    M_max = 20;
    Nscan = 1000;

    Mvec = linspace(M_min, M_max, Nscan);
    Fvec = f(Mvec);

    idx = find(Fvec(1:end-1).*Fvec(2:end) < 0, 1);

    if isempty(idx)
        % No real solution → physically impossible (detached region)
        M = NaN;
        return;
    end

    bracket = [Mvec(idx), Mvec(idx+1)];
    M = fzero(f, bracket);

end