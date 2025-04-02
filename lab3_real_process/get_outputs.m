% addpath('C:\vaja_3\simulink');

%% clearing
clear all 
close all
%% input signals:
% Number of frequencies to test
num_freqs = 11;

% Define the frequency range
frequencies = logspace(log10(0.05), log10(0.2), num_freqs); % from 0.05 Hz to 0.2 Hz

% Time vector and sampling time
Ts = 0.1; 
t = 0:Ts:150;

% Preallocate arrays for input and output
x = zeros(length(t), num_freqs); 
x1 = zeros(length(t), num_freqs); 
y = zeros(length(t), num_freqs);

% Generate sinusoidal input signals
for i = 1:num_freqs
    x(:, i) = cos(2 * pi * frequencies(i) * t) + 3;
    x1(:, i) = sin(2 * pi * frequencies(i) * t) + 3;
end
%% getting output one by one and saving them
Ts = 0.1; 
modelName = 'proces';
load_system(modelName);

i = 11; %%% nastavljaj ročno !!!

x = x(:, i); % vhodni signal je sinus
x1 = x1(:, i); % Za izracun bomo rabili se kosinus

% to gre v simulink:
x_inp = [t', x]; 
x1_inp = [t', x1];
%% plot output and input
figure;
plot(out.t, x(:, 1)); 
hold on;
plot(out.t, out.y_out); 
xlabel('t [s]');
ylabel('u(t), y(t)');
title(sprintf('Frequency = %.3f Hz', frequencies(i)));
grid on;









