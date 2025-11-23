%% Clear Workspace

clear; clc; close all;

%% Part 1 - Test Case

% See the solution struct for the results if you edit this test case.
% All units are either dimensionless or in [m] and [rad].
b = 10;
a0_t = 6.0;
a0_r = 6.0;
c_t = 0.75;
c_r = 1.25;
aero_t = -0.035;
aero_r = -0.035;
geo_t = 0.0873;
geo_r = 0.0873;
N = 20;

% Calculate test case
[solution.e, solution.c_L, solution.c_Di] = PLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t, geo_r, N);

% Display test case results
disp('Test Case Results:');
disp(' ');
disp(['e = ' num2str(solution.e)]);
disp(['c_L = ' num2str(solution.c_L)]);
disp(['c_Di = ' num2str(solution.c_Di)]);

%% Part 2 - Figure 5.20 Replication

% Creating AR and taper cases
AR = [4 6 8 10];
TR = linspace(0, 1, 50);

% Preallocate the results matrix
this_e = zeros(length(AR), length(TR)); 

% Calculate span efficiency factor for every case
for AR_i = 1:length(AR)

    for TR_i = 1:length(TR)

        % Calculate taper and wingspan based on AR and ratio
        c_r = 1;
        c_t = TR(TR_i) * c_r;
        b = (AR(AR_i) * (c_t + c_r)) / 2;

        % Calculate span efficiency factor
        [this_e(AR_i, TR_i), ~, ~] = PLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t, geo_r, N);

    end

end

% Converting span efficiency factor to delta
this_delta = (1 ./ this_e) - 1;

% Plotting delta with span efficiency factor
figure();
hold on
grid on

% Metadata generation
title('Delta vs taper ratio');
xlabel('Taper ratio, $\frac{c_t}{c_r}$', 'Interpreter', 'latex');
ylabel('$\delta$', 'Interpreter', 'latex', Rotation=0);
legend('show', 'Location', 'best');

% Iteratively plotting by aspect ratio
for AR_i = 1:length(AR)

    plot(TR, this_delta(AR_i, :), 'DisplayName', ['AR = ' num2str(AR(AR_i))]);
    
end

% Create figure (for team, will be commented in submission)
% print(gcf, 'Figure_5-20_Replication.png', '-dpng', '-r300');

%% PLLT Function Definition

% Name:         PLLT 
% Description:  Calculates lift, induced drag, and efficiency factor through PLLT
% Inputs:       Geometric properties of wing
% Outputs:      Lift, induced drag, and efficiency factor for wing
function [e, c_L, c_Di] = PLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t, geo_r, N)

    arguments (Input)
        b       % Wingspan
        a0_t    % Sectional lift slope at tip (1/rad)
        a0_r    % Sectional lift slope at root (1/rad)
        c_t     % Chord length at tip (m)
        c_r     % Chord length at root (m)
        aero_t  % Zero-lift AoA at tip (rad)
        aero_r  % Zero-lift AoA at root (rad)
        geo_t   % Geometric AoA at root (Geo twist + alpha, rad)
        geo_r   % Geometric AoA at root (Geo twist + alpha, rad)
        N       % Number of odd terms in circulation expansion
    end
    
    arguments (Output)
        e       % Span efficiency factor
        c_L     % Coefficient of lift
        c_Di    % Coefficient of induced drag
    end
    
    % Defining properties of wing
    %
    % Finding the proper way to convert theta -> y for the sake of this linear
    % interpolation was lengthy... but makes sense since it's not just some
    % direct conversion via y = -b/2 * cos(theta), but rather the factor of
    % linearization.
    AR = (2*b)/(c_t+c_r);
    theta = ((1:N)' * pi) / (2 * N);
    linfactor = cos(theta);
    
    % Linear interpolation of properties
    myChord =       c_r         + (c_t - c_r) * linfactor;
    myZeroLift =    aero_r      + (aero_t - aero_r) * linfactor;
    myAoA =         geo_r       + (geo_t - geo_r) * linfactor;
    myLiftCurve =   a0_r        + (a0_t - a0_r) * linfactor;
    
    % Building solution matrix
    M = zeros(N,N);
    for j = 1:N
        n = 2*j - 1;
        M(:,j) = ((4*b .* sin(n .* theta)) ./ (myLiftCurve .* myChord)) + n .* (sin(n .* theta) ./ sin(theta));
    end
    
    % Solving for coefficients
    A = M \ (myAoA - myZeroLift);
    
    % Calculating delta
    myDelta = 0;
    for j = 2:N
        n = 2*j - 1;
        myDelta = myDelta + n * (A(j)/A(1))^2;
    end
    
    % Final calculations
    c_L = A(1) * pi * AR;
    c_Di = (c_L^2 * (1+myDelta))/(pi * AR);
    e = 1/(1+myDelta);

end