function loaddata

%this program loads the experimental TGA data from the .xls worksheet

clc
close all
clear all

data = xlsread('MC carpel H 5-min'); % the H 5-min stands for the heating rate of 5 K/min
time = data(:,1);
Temp  = data(:,2);
mexp  = data(:,3);
dmexp = data(:,4);
xexp = (10.32687*ones(1,length(mexp))'- mexp)/(10.32687-3.45199); % (mi-mt)/(mi-mf)
dxdt = dmexp/(10.32687-3.45199);

C=[time mexp dmexp xexp dxdt Temp];
save data_carpel_5.txt C -ASCII 
