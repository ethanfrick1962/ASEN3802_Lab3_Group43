%% Clear Workspace

clear; clc; close all;

%% Part 1 - Test Case

% See the solution struct for the results if you edit this test case.
% All units are either dimensionless or in [m] and [rad].

% Angle of attack (5 degrees)
aoa = 5 * pi/180;

% Wing parameters
b = 10.9728;
a0_t = 0.1191*180/pi;
a0_r = 0.1192*180/pi;
c_t = 1.092;
c_r = 1.627632;
aero_t = 0;
aero_r = deg2rad(-2.1402);
geo_t = 0;
geo_r = 0;
N = 50;

% Calculate test case
for i = 1:1:401
geo_t = pi/180 * (-5 + (i * 0.05));
geo_r = pi/180 * (-3 + (i * 0.05));
[solution.e(i), solution.c_L(i), solution.c_Di(i)] = PLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t, geo_r, N);
end

% Display test case results
disp('Test Case Results:');
disp(' ');
disp(['e = ' num2str(solution.e)]);
disp(['c_L = ' num2str(solution.c_L)]);
disp(['c_Di = ' num2str(solution.c_Di)]);

%% Part 2 - Figure 5.20 Replication

% % Creating AR and taper cases
% AR = [4 6 8 10];
% TR = linspace(0, 1, 50);
% 
% % Preallocate the results matrix
% this_e = zeros(length(AR), length(TR)); 
% 
% % Calculate span efficiency factor for every case
% for AR_i = 1:length(AR)
% 
%     for TR_i = 1:length(TR)
% 
%         % Calculate taper and wingspan based on AR and ratio
%         c_r = 1;
%         c_t = TR(TR_i) * c_r;
%         b = (AR(AR_i) * (c_t + c_r)) / 2;
% 
%         % Calculate span efficiency factor
%         [this_e(AR_i, TR_i), ~, ~] = PLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t, geo_r, N);
% 
%     end
% 
% end
% 
% % Converting span efficiency factor to delta
% this_delta = (1 ./ this_e) - 1;
% 
% % Plotting delta with span efficiency factor
% % figure();
% hold on
% grid on
% 
% % Metadata generation
% title('Delta vs taper ratio');
% xlabel('Taper ratio, $\frac{c_t}{c_r}$', 'Interpreter', 'latex');
% ylabel('$\delta$', 'Interpreter', 'latex', Rotation=0);
% legend('show', 'Location', 'best');
% 
% % Iteratively plotting by aspect ratio
% for AR_i = 1:length(AR)
% 
%     % plot(TR, this_delta(AR_i, :), 'DisplayName', ['AR = ' num2str(AR(AR_i))]);
% 
% end
% 
% % Create figure (for team, will be commented in submission)
% % print(gcf, 'Figure_5-20_Replication.png', '-dpng', '-r300');

%% Code
data_0012 = load("0012 cdvscl.csv");
data_2412 = load("2412 cdvscl.csv");

X0012_cl = data_0012(1:end, 1);
X0012_cd = data_0012(1:end, 2);
X2412_cl = data_2412(1:end, 1);
X2412_cd = data_2412(1:end, 2);

cd_fit_0012 = polyfit(X0012_cl,X0012_cd,2);
cd_fit_2412 = polyfit(X2412_cl,X2412_cd,2);

AoA = -5:0.05:15;
S = 23112;

y = -18*12:1:18*12;
c = linspace(64,43,round(length(y)/2));
chalf = c(1,2:end);
c = cat(2,flip(chalf),64,chalf);

a_0012 = 0.1060;
a_2412 = 0.1041;

a = linspace(a_2412,a_0012,round(length(y)/2));
ahalf = a(1,2:end);
a = cat(2,flip(ahalf),a_2412,ahalf);

aL0_0012 = 0.0004;
aL0_2412 = -2.1402;

aL0 = linspace(aL0_2412,aL0_0012,round(length(y)/2));
aL0half = aL0(1,2:end);
aL0 = cat(2,flip(aL0half),aL0_2412,aL0half);

twist = linspace(0,2,round(length(y)/2));
twisthalf = twist(1,2:end);
twist = cat(2,flip(twisthalf),0,twisthalf);

local_cl = zeros(size(y));
local_cd = zeros(size(y));
weight = linspace(0,1,round(length(y)/2));
weighthalf = weight(1,2:end);
weight = cat(2,flip(weighthalf),0,weighthalf);
total_cd = zeros(size(AoA));

for i = 1:length(AoA)
    local_AoA = twist + AoA(i);
    for j = 1:length(y)
        local_cl(j) = (local_AoA(j) - aL0(j)) * a(j);
        local_cd0012 = cd_fit_0012(1) * local_cl(j) ^ 2 + cd_fit_0012(2) * local_cl(j) + cd_fit_0012(3);
        local_cd2412 = cd_fit_2412(1) * local_cl(j) ^ 2 + cd_fit_2412(2) * local_cl(j) + cd_fit_2412(3);
        local_cd(j) = weight(j) * local_cd0012 + (1 - weight(j)) * local_cd2412;
    end
    cd_integral = trapz(local_cd .* c);
    total_cd(i) = cd_integral / S;
end

figure()
plot(AoA,total_cd, LineWidth=0.8);
hold on
grid on 
title("Coefficient of Drag Analysis")
plot(AoA, solution.c_Di, LineWidth=0.8)
plot(AoA, solution.c_Di + total_cd, LineWidth=0.8)
legend("C_{Do}", "C_{Di}", "C_D", Location="northwest")
xlabel("Angle of Attack(Deg)")
ylabel("C_{Dx}")
print(gcf, 'CD_Cessna.png', '-dpng', '-r300');

figure()
plot(AoA, solution.c_L, LineWidth=0.8)
hold on
grid on
title("Coefficient of Lift")
xlabel("Angle of Attack(Deg)")
ylabel("C_L")
print(gcf, 'CL_Cessna.png', '-dpng', '-r300');

W_cessna = 2500;
rho = 0.001756;
S_cessna = 0.5*b*(c_t + c_r)*(3.28084)^2;
V_cessna = sqrt((W_cessna .* 2)./(rho .* solution.c_L .* S_cessna));
V_cessna = V_cessna(54:end);
q = 0.5*rho.*V_cessna.^2;
T_cessna = (solution.c_Di(54:end) + total_cd(54:end)) .* q .* S_cessna;
figure()
plot(AoA(54:end), T_cessna)
hold on 
plot(AoA(54:end), V_cessna)
legend("Thrust", "Velocity")

figure()
plot(V_cessna, T_cessna)


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