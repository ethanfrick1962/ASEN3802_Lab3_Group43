%% Clear Workspace
clear; clc; close all;

%% Part 1 - Test Case

% Angle of attack (radians)
aoa = 5 * pi/180;

% Wing parameters
b = 10.9728;          % Span [m]
a0_t = 0.1191*180/pi; % Tip lift slope
a0_r = 0.1192*180/pi; % Root lift slope
c_t = 1.092;           % Tip chord [m]
c_r = 1.627632;        % Root chord [m]
aero_t = 0;            % Tip zero-lift AoA [rad]
aero_r = deg2rad(-2.1402); % Root zero-lift AoA [rad]
geo_t = 0;             % Tip geometric twist [rad]
geo_r = 0;             % Root geometric twist [rad]
N = 50;                % PLLT terms

% Solve PLLT across AoA sweep
AoA = -5:0.05:15; % degrees
solution.e = zeros(size(AoA));
solution.c_L = zeros(size(AoA));
solution.c_Di = zeros(size(AoA));

for i = 1:length(AoA)
    geo_t_i = geo_t + deg2rad(AoA(i));
    geo_r_i = geo_r + deg2rad(AoA(i));
    [solution.e(i), solution.c_L(i), solution.c_Di(i)] = PLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t_i, geo_r_i, N);
end

%% Load airfoil CD vs CL data
data_0012 = load("0012 cdvscl.csv");
data_2412 = load("2412 cdvscl.csv");

X0012_cl = data_0012(:,1);
X0012_cd = data_0012(:,2);
X2412_cl = data_2412(:,1);
X2412_cd = data_2412(:,2);

cd_fit_0012 = polyfit(X0012_cl,X0012_cd,2);
cd_fit_2412 = polyfit(X2412_cl,X2412_cd,2);

%% Wing discretization
y = linspace(-b/2, b/2, 101);  % spanwise stations
c = linspace(c_r, c_t, length(y)); % chord distribution
a = linspace(0.1041, 0.1060, length(y)); % lift slope
aL0 = linspace(deg2rad(-2.1402), 0.0004, length(y)); % zero-lift AoA
twist = linspace(0, deg2rad(2), length(y)); % geometric twist
weight = linspace(0,1,length(y)); % blending weight for airfoils

%% Compute profile drag for each AoA
total_cd = zeros(size(AoA));

for i = 1:length(AoA)
    local_AoA = twist + deg2rad(AoA(i)); % radians
    local_cl = (local_AoA - aL0) .* a;
    
    % Clip negative local_cl to 0 to avoid NaNs (if AoA below zero-lift)
    local_cl(local_cl < 0) = 0;
    
    local_cd0012 = polyval(cd_fit_0012, local_cl);
    local_cd2412 = polyval(cd_fit_2412, local_cl);
    local_cd = weight .* local_cd0012 + (1-weight) .* local_cd2412;
    
    % Integrate along span
    cd_integral = trapz(y, local_cd .* c);
    S = trapz(y, c); % total wing area
    total_cd(i) = cd_integral / S;
end

%% Plot Drag vs AoA
figure();
plot(AoA, total_cd, 'LineWidth', 1.2); hold on; grid on;
plot(AoA, solution.c_Di, 'LineWidth', 1.2);
plot(AoA, solution.c_Di + total_cd, 'LineWidth', 1.2);
xlabel('Angle of Attack (deg)'); ylabel('C_D');
title('Drag Coefficient Analysis');
legend('C_{Do}','C_{Di}','C_D','Location','northwest');

%% Plot Lift vs AoA
figure();
plot(AoA, solution.c_L, 'LineWidth',1.2); grid on;
xlabel('Angle of Attack (deg)'); ylabel('C_L');
title('Lift Coefficient');

%% Compute Minimum Airspeed
W_cessna = 2500; % lb
rho = 0.001756; % slugs/ft^3
S_cessna = trapz(y, c)*(3.28084)^2; % ft^2

[C_L_max, idx_max] = max(solution.c_L); % maximum lift
V_min = sqrt((2*W_cessna)/(rho*S_cessna*C_L_max)); % ft/s
AoA_min = AoA(idx_max); % degrees
C_D_total = total_cd(idx_max) + solution.c_Di(idx_max);
q = 0.5*rho*V_min^2;
T_min = C_D_total * q * S_cessna; % lb

disp(['Minimum Airspeed: ', num2str(V_min), ' ft/s']);
disp(['AoA at Minimum Airspeed: ', num2str(AoA_min), ' deg']);
disp(['Thrust required: ', num2str(T_min), ' lb']);

%% Proper Thrust vs Velocity Curve
V_range = linspace(50, 250, 200); % ft/s
T_vs_V = zeros(size(V_range));

% Interpolate total drag coefficient vs CL
CL_array = solution.c_L;
CD_array = solution.c_Di + total_cd;

for i = 1:length(V_range)
    V = V_range(i);
    
    % Required CL to support weight at this V
    CL_req = (2*W_cessna) / (rho * S_cessna * V^2);
    
    % Interpolate corresponding drag coefficient
    if CL_req <= max(CL_array) && CL_req >= min(CL_array)
        CD_total = interp1(CL_array, CD_array, CL_req, 'linear');
    else
        CD_total = NaN; % Out of range
    end
    
    % Compute thrust required
    q = 0.5 * rho * V^2;
    T_vs_V(i) = CD_total * q * S_cessna;
end

%% Plot

figure();
plot(V_range, T_vs_V, 'LineWidth',1.5);
grid on;
xlabel('Velocity (ft/s)');
ylabel('Thrust Required (lb)');
title('Thrust Required vs Velocity (Level Flight)');


%% PLLT Function
function [e, c_L, c_Di] = PLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t, geo_r, N)
    AR = (2*b)/(c_t + c_r);
    theta = ((1:N)' * pi) / (2*N);
    linfactor = cos(theta);
    
    myChord = c_r + (c_t - c_r) * linfactor;
    myZeroLift = aero_r + (aero_t - aero_r) * linfactor;
    myAoA = geo_r + (geo_t - geo_r) * linfactor;
    myLiftCurve = a0_r + (a0_t - a0_r) * linfactor;
    
    M = zeros(N,N);
    for j = 1:N
        n = 2*j-1;
        M(:,j) = ((4*b .* sin(n .* theta)) ./ (myLiftCurve .* myChord)) + n * (sin(n .* theta) ./ sin(theta));
    end
    
    A = M \ (myAoA - myZeroLift);
    
    myDelta = 0;
    for j = 2:N
        n = 2*j-1;
        myDelta = myDelta + n*(A(j)/A(1))^2;
    end
    
    c_L = A(1)*pi*AR;
    c_Di = (c_L^2 * (1 + myDelta))/(pi*AR);
    e = 1/(1+myDelta);
end
