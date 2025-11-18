clear
clc
close all

data_0012 = load("0012 cdvscl.csv");
data_2412 = load("2412 cdvscl.csv");

X0012_cl = data_0012(1:end, 1);
X0012_cd = data_0012(1:end, 2);
X2412_cl = data_2412(1:end, 1);
X2412_cd = data_2412(1:end, 2);

cd_fit_0012 = polyfit(X0012_cl,X0012_cd,2);
cd_fit_2412 = polyfit(X2412_cl,X2412_cd,2);

AoA = -5:0.01:8;
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


plot(AoA,total_cd);

hold on
