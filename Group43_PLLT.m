%% Clear Workspace

clear; clc; close all;

%% Test Case

% Recall order: b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t, geo_r, N
[e, c_L, c_Di] = PLLT(10, 6.0, 6.0, 0.75, 1.25, -0.035, -0.035, 0.0873, 0.0873, 5);

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
    AR = b/((c_t+c_r)/2);
    theta = (1:N)' * pi / (2 * N);
    linfactor = cos(theta);
    
    % Linear interpolation of properties
    myChord =       c_r         + (c_t - c_r) * abs(linfactor);
    myZeroLift =    aero_r      + (aero_t - aero_r) * abs(linfactor);
    myAoA =         geo_r       + (geo_t - geo_r) * abs(linfactor);
    myLiftCurve =   a0_r        + (a0_t - a0_r) * abs(linfactor);
    
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