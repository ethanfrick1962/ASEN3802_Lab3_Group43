function [XB, YB, YC] = Airfoil_Generator(MAX_CAMBER, LOC_CAMBER, MAX_THICK, C, RES)
% Name:     AIRFOIL_GENERATOR
% Summary:  Generates 4-digit NACA airfoils.

    %% Input and Output Arguments

    arguments (Input)
        MAX_CAMBER
        LOC_CAMBER
        MAX_THICK
        C
        RES
    end

    arguments (Output)
        XB
        YB
        YC
    end

    %% Defining Key Parameters

    % Linspace of x and x/c
    X = linspace(0, C, RES);
    X_OVER_C = X / C;

    % Defining m, t, and p
    M = MAX_CAMBER / 100;
    P = LOC_CAMBER / 10;
    T = MAX_THICK / 100;

    %% Calculating y-positions

    % Calculating thickness y-position
    YT = (T / 0.2) * C * ( ...
        0.2969 * sqrt(X_OVER_C) ...
      - 0.1260 * (X_OVER_C) ...
      - 0.3516 * (X_OVER_C).^2 ...
      + 0.2843 * (X_OVER_C).^3 ...
      - 0.1036 * (X_OVER_C).^4 );

    % Kutta condition
    YT(end) = 0;
    
    % Setting camber y-position
    YC = zeros(size(X_OVER_C));

    %  Yc with x < p*c and x >= p*c
    for I = 1:length(X)
        if X_OVER_C(I) < P
            YC(I) = ((M / (P^2)) * X(I)) * (2*P - X_OVER_C(I));
        else
            YC(I) = ((M / ((1 - P)^2)) * (C - X(I))) .* (1 + X_OVER_C(I) - 2*P);
        end
    end

    %% Calculating upper and lower x and y

    % Ksi
    ANGLE = atan(diff(YC ./ X));
    ANGLE(isnan(ANGLE)) = 0;
    ANGLE(end + 1) = 0;

    % Upper and lower x-positions
    XU = X - (YT .* sin(ANGLE));
    XL = X + (YT .* sin(ANGLE));

    % Upper and lower y-positions
    YU = YC + (YT .* cos(ANGLE));
    YL = YC - (YT .* cos(ANGLE));

    XB = [XU fliplr(XL)];
    YB = [YU fliplr(YL)];

    XB = XB([RES+1:RES*2, 2:RES]);
    YB = YB([RES+1:RES*2, 2:RES]);
    
end