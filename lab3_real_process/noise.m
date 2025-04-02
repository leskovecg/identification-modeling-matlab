addpath('C:\vaja_3\simulink');

%% clearing
clear all 
close all

%% nois:
modelName = 'nois';
load_system(modelName);

constant = ones(1501,1) * 3;
Ts = 0.1; 
t = 0:Ts:150;
constant_input = [t', constant];

n = zeros(length(t), 1);

%% plotanje

figure;
plot(out.t, constant); 
hold on;
plot(out.t, out.n); 
xlabel('t [s]');
ylabel('u(t), y(t)');
title(sprintf('noise'));
grid on;

%% fft:

N = fft(out.n);

figure;
plot(out.t, N); 
xlabel('t [s]');
ylabel('N(w), y(t)');
title(sprintf('noise'));
grid on;

