clear
clc
close all

% NACA 0012
[XB1,YB1,YC1] = Airfoil_Generator(0,0,12,1,131);
% NACA 2412
[XB2,YB2,YC2] = Airfoil_Generator(2,4,12,1,131);
% NACA 4412
[XB3,YB3,YC3] = Airfoil_Generator(4,4,12,1,131);

AoA = -8:0.1:8;
cL_0012 = zeros(size(AoA));
cL_2412 = zeros(size(AoA));
cL_4412 = zeros(size(AoA));
for i = 1:length(AoA)
    cL_0012(i) = Vortex_Panel(XB1,YB1,AoA(i));
end
for i = 1:length(AoA)
    cL_2412(i) = Vortex_Panel(XB2,YB2,AoA(i));
end
for i = 1:length(AoA)
    cL_4412(i) = Vortex_Panel(XB3,YB3,AoA(i));
end

plot(AoA,cL_0012);
hold on
grid on
plot(AoA,cL_2412);
plot(AoA,cL_4412);
legend('NACA 0012', 'NACA 2412', 'NACA 4412', location='southeast');
ylabel('Coefficient of Lift');
xlabel('Angle of Attack (degrees)');
title('Coefficient of Lift vs. AoA for NACA 0012, 2412, 4412')
xlim('padded')
ylim('padded')
%print('Part 3', '-dpng', '-r300');

slope_0012 = (cL_0012(150) - cL_0012(5)) / (AoA(150) - (AoA(5)));
slope_2412 = (cL_2412(150) - cL_2412(5)) / (AoA(150) - (AoA(5)));
slope_4412 = (cL_4412(150) - cL_4412(5)) / (AoA(150) - (AoA(5)));