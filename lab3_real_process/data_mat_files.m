addpath('C:\Users\Gasper\OneDrive\FAKS\MAGISTERIJ\letnik_1\SEMESTER_2\MODUL_A_IDENTIFIKACIJA\lab_vaje\vaja_3\data')
%% clearing:
clear all;
close all;
%% visualize data from .mat files:
simOut = load('simOut.mat'); %% ta je slaba
simOutEna = load('simOutEna.mat'); %% to mamo v nonparametric

figure;
plot(simOut.simOut.t, simOut.simOut.Y);
title('simOut.mat');

figure;
plot(simOutEna.simOut_ena.t, simOutEna.simOut_ena.Y);
title('simOutEna.mat');



