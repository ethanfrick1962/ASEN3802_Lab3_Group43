clear
clc
close all

% load in the cd vs cl curve data
data_0012 = load("0012 cdvscl.csv");
data_2412 = load("2412 cdvscl.csv");

% extract the cl and cd values for each airfoil
X0012_cl = data_0012(1:end, 1);
X0012_cd = data_0012(1:end, 2);
X2412_cl = data_2412(1:end, 1);
X2412_cd = data_2412(1:end, 2);

% fit the drag curves for each airfoil
cd_fit_0012 = polyfit(X0012_cl,X0012_cd,2);
cd_fit_2412 = polyfit(X2412_cl,X2412_cd,2);

% create an AoA vector to sweep through
AoA = -5:0.01:8;

% area in square inches
S = 23112;

% create wingspan vector in inches
y = -18*12:1:18*12;

% create chord vector in inches
c = linspace(64,43,round(length(y)/2));
% symmetric wing
chalf = c(1,2:end);
% join halves and insert middle value
c = cat(2,flip(chalf),64,chalf);

% lift slope from previous part
% taken fron the tables in the previous report
a_0012 = 0.1060;
a_2412 = 0.1041;

% create lift slope vector
% assume lift slope varies linearly between that of 0012 and 2412 from tip
% to root
a = linspace(a_2412,a_0012,round(length(y)/2));
% symmetric wing
ahalf = a(1,2:end);
% join halves and insert middle value
a = cat(2,flip(ahalf),a_2412,ahalf);

% zero lift AoA from previous part
% taken from the tables in the previous report
aL0_0012 = 0.0004;
aL0_2412 = -2.1402;

% create zero lift AoA vector
% assume zero lift AoA varies linearly between that of 0012 and 2412 from
% tip to root
aL0 = linspace(aL0_2412,aL0_0012,round(length(y)/2));
% symmetric wing
aL0half = aL0(1,2:end);
% join halves and insert middle value
aL0 = cat(2,flip(aL0half),aL0_2412,aL0half);

% create twist vector
twist = linspace(0,2,round(length(y)/2));
% symmetric wing
twisthalf = twist(1,2:end);
% join halves and insert middle value
twist = cat(2,flip(twisthalf),0,twisthalf);

% preallocate variables
local_cl = zeros(size(y));
local_cd = zeros(size(y));

% create cd weighting vector to average the airfoil values (explanation in
% loop)
weight = linspace(0,1,round(length(y)/2));
% symmetric wing
weighthalf = weight(1,2:end);
% join halves and insert middle value
weight = cat(2,flip(weighthalf),0,weighthalf);

% preallocate variable
total_cd = zeros(size(AoA));

% work through the AoA vector
for i = 1:length(AoA)
    % get a vector of local AoA based on twist at each point + wing AoA
    local_AoA = twist + AoA(i);

    % work across the wingspan
    for j = 1:length(y)

        % get the cL at each point along the span based on the local AoA,
        % local lift slope, and local zero lift AoA
        local_cl(j) = (local_AoA(j) - aL0(j)) * a(j);

        % get the cd for 0012 at this cL
        local_cd0012 = cd_fit_0012(1) * local_cl(j) ^ 2 + cd_fit_0012(2) * local_cl(j) + cd_fit_0012(3);
        % get the cd for 2412 at this cL
        local_cd2412 = cd_fit_2412(1) * local_cl(j) ^ 2 + cd_fit_2412(2) * local_cl(j) + cd_fit_2412(3);

        % take a weighted average of the two cd values to get an estimate
        % of the actual cd
        % the weighting assumes the airfoil is 100% 0012 at the tip, 100%
        % 2412 at the root, and a linear mix of the two airfoils in between
        local_cd(j) = weight(j) * local_cd0012 + (1 - weight(j)) * local_cd2412;
    end
    % integrate the local cd across the wingspan
    cd_integral = trapz(local_cd .* c);
    % divide by the wing area to get profile drag coefficient at this AoA
    total_cd(i) = cd_integral / S;
end

% plot the profile drag coefficient as a function of AoA
plot(AoA,total_cd);
grid on
hold on