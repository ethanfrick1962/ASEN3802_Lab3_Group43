%% Clear Workspace

clear; clc; close all;

%% Generate airfoil output

my_res = 100;
my_chord = 1;

[x_0018, y_0018, yc_0018] = Airfoil_Generator(0, 0, 18, my_chord, my_res);
[x_2418, y_2418, yc_2418] = Airfoil_Generator(2, 4, 18, my_chord, my_res);
[x_0012, y_0012, yc_0012] = Airfoil_Generator(0, 0, 12, my_chord, my_res);


%% Plot

% Plotting NACA 0018
figure();
hold on
grid on

title("NACA 0018 Airfoil")
plot(x_0018, y_0018);
% plot(x_0018(1:my_res), yc_0018);
xlim([-0.2 my_chord + 0.2]);
ylim([-1 1]);

% Plotting NACA 2418
figure();
hold on
grid on

title("NACA 2418 Airfoil")
plot(x_2418, y_2418);
plot(x_2418(1:my_res), yc_2418);
xlim([-0.2 my_chord + 0.2]);
ylim([-1 1]);


%% Vortex Panel Output

% CL = Vortex_Panel(x_0012, y_0012, 5);

%% Task II

% CL Determination
true_Res = 50;
my_res2 = floor(linspace(10, 1000, true_Res));
CL2 = zeros(1, true_Res);
[x_0012, y_0012, yc_0012] = Airfoil_Generator(0, 0, 12, my_chord, 2500);
CLTrue = Vortex_Panel(x_0012, y_0012, 5);

for i = 1:true_Res
    [x_0012, y_0012, yc_0012] = Airfoil_Generator(0, 0, 12, my_chord, my_res2(i));
    CL2(i) = Vortex_Panel(x_0012, y_0012, 5);

end
%%
my_res2 = [my_res2, 2500];
CL2 = [CL2, CLTrue];

%% Plotting NACA 0012 Cl
figure();
hold on
grid on

title("CL by Panel Number")
xlabel("Number of Panels")
ylabel("Coefficient of Lift")
plot(my_res2, CL2);

% Find approximate resolution
approxTarget = 0.99*CL2(51);
% Determine the index of the closest resolution to the target
[~, targetIndex] = min(abs(CL2 - approxTarget));
panelRes = my_res2(targetIndex);