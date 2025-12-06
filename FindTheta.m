function theta = FindTheta(Beta, M1, gamma)

    theta = abs(atand(2 .* cotd(Beta) .* ((M1.^2 .* sind(Beta).^2 - 1) ./...
         (M1.^2 .* (gamma + cosd(2 .* Beta)) + 2))));
end