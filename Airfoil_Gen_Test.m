%% Clear Workspace

clear; clc; close all;

%% Adding functions

addpath('func');

%% Generate airfoil output

my_res = 100;
my_chord = 1;

[x_0018, y_0018, yc_0018] = Airfoil_Generator(0, 0, 18, my_chord, my_res);
[x_2418, y_2418, yc_2418] = Airfoil_Generator(2, 4, 18, my_chord, my_res);

%% Plot

% Plotting NACA 0018
figure();
hold on
grid on

plot(x_0018, y_0018);
plot(x_0018(1:my_res), yc_0018);
xlim([-0.2 my_chord + 0.2]);
ylim([-1 1]);

% Plotting NACA 2418
figure();
hold on
grid on

plot(x_2418, y_2418);
plot(x_2418(1:my_res), yc_2418);
xlim([-0.2 my_chord + 0.2]);
ylim([-1 1]);
