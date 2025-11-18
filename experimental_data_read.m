clear
clc
close all

data = load("0006.csv");

AoA = data(:,1);
cL = data(:,2);

p = polyfit(AoA,cL,1);

