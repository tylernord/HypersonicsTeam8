function Beta = FindBeta(theta, gamma, M1)

    beta = linspace(0, deg2rad(90), 1000);
    theta_target = deg2rad(theta);
    
    % Convert theta from degrees to radians for MATLAB
    
    for i = 1:length(beta)
        theta_new(i) = atan(2 .* cot(beta(i)) .* ((M1.^2 .* sin(beta(i)).^2 - 1) ./...
            (M1.^2 .* (gamma + cos(2 .* beta(i))) + 2)));
        if abs(theta_new(i) - theta_target) < 10e-3
            index = i;
            break
        end
    end
    
    Beta = rad2deg(beta(i));

end