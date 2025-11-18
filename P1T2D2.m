%% Clear Workspace

clear; clc; close all;

%% Plot NACA Airfoils w/ AoA Variation

% Define the range of AoA and airfoil parameters
DAOA = 0.5;
AOA = -15:DAOA:15;
AIRFOIL_THICKNESS = [6; 12; 18];
CHORD = 1;
RES_OPTIMAL = 131;

% Put all into one plot
figure();
hold on
grid on

% Loop through airfoils and create AoA plot
for I = 1:length(AIRFOIL_THICKNESS)

    % Generate airfoil based on thickness
    [XB, YB, YC] = Airfoil_Generator(0, 0, AIRFOIL_THICKNESS(I), CHORD, RES_OPTIMAL);

    for ANGLE = 1:length(AOA)

        % Calculate lift and drag coefficients for the current angle of attack
        CL(ANGLE) = Vortex_Panel(XB, YB, AOA(ANGLE));

    end

    % Plot coefficient of lift for this angle of attack
    plot(AOA, CL);

    % Truncate to four decimal places
    N = 4;

    % Calculate zero-lift angle via data
    ZERO_LIFT_ANGLE_EXP(I) = fix(CL(AOA == 0) * 10^N) / (10^N);

    % Calculating lift slope via data
    LIFT_SLOPE_EXP(I) = mean(diff(CL(~AOA == 0) / DAOA));

end

% Plot metadata
title('NACA airfoils (0006, 0012, 0018) coefficient of lift compared over angle of attack');
xlabel('Angle of attack (°)');
ylabel('Coefficient of lift (dimensionless)');
legend('NACA 0006', 'NACA 0012', 'NACA 0018', 'Location', 'best');

% Save figure
print(gcf, 'NACA_Airfoil_Uncambered.png', '-dpng', '-r300'); 