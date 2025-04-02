
addpath('C:\Users\Uporabnik\OneDrive\FAKS\MAGISTERIJ\letnik_1\SEMESTER_2\MODUL_A_IDENTIFIKACIJA\lab_vaje\vaja_3\results\')
%% clearing
clear all;
close all;
%% input signals:
num_freqs = 11;
frequencies = logspace(log10(0.05), log10(0.2), num_freqs); % from 0.05 Hz to 0.2 Hz

Ts = 0.1; 
t = 0:Ts:150;
t = t';

x = zeros(length(t), num_freqs); 
x1 = zeros(length(t), num_freqs); 

for i = 1:num_freqs
    x(:, i) = cos(2 * pi * frequencies(i) * t);
    x1(:, i) = sin(2 * pi * frequencies(i) * t);
end
%% plot all signals:
figure;
for i = 1:11
    i_s = num2str(i);
    result = ['out_', i_s, '.mat'];
    load(result);
   
    subplot(3, 4, i);
    plot(out.t, x(:, i));
    hold on;
    plot(out.t, x1(:, i));
    hold on;
    plot(out.t, out.y_out);
    hold on;
    xlabel('t [s]');
    ylabel('x(t), x1(t), y(t)');
    title(sprintf('Frequency = %.3f Hz', frequencies(i)));
    grid on;
end
%% cut and plot signals 
figure;
for i = 1:11
    
    i_s = num2str(i);
    result = ['out_', i_s, '.mat'];
    load(result);
    
    whole_time = 150;
    cutted_time = 30;
    time_left = whole_time - cutted_time;

    T = 1 / frequencies(i);
    
    periods = floor(time_left / T);

    time_of_periods = periods * T;
    samples = time_of_periods / Ts;
    
    start = floor(length(out.t) - samples);

    t = out.t(start:end, 1);
    neki = out.y_out(start:end, 1);
    vrenost = mean(neki);


    subplot(3, 4, i);
    plot(t, x(start:end, i)); 
    hold on;
    plot(t, x1(start:end, i)); 
    hold on;
    plot(t, out.y_out(start:end, 1));
    hold on;
    plot(t, out.y_out(start:end, 1) - vrenost);
    xlabel('t [s]');
    ylabel('x(t), x1(t), y(t), y(t) - mean');
    title(sprintf('Frequency = %.3f Hz', frequencies(i)));
    legend('x(t)', 'x1(t)', 'y(t)', 'y(t) - mean');
    grid on;

    lege = legend;
    lege.Position = [0.6, 0.1, 0.3, 0.2];
end


