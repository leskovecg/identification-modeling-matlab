
%% clearing:
clear all;
close all;
%% signals:
num_freqs = 11;
frequencies = logspace(log10(0.05), log10(0.2), num_freqs); % from 0.05 Hz to 0.2 Hz

Ts = 0.1; 
t = 0:Ts:150;

x = zeros(length(t), num_freqs); 
x1 = zeros(length(t), num_freqs);
x_moved= zeros(length(t), num_freqs); 

for i = 1:num_freqs
    x(:, i) = cos(2 * pi * frequencies(i) * t);
    x1(:, i) = sin(2 * pi * frequencies(i) * t);
    x_moved(:, i) = cos(2 * pi * frequencies(i) * t) + 3; 
end
%% plot signals:
figure;

for i = 1:11
    i_s = num2str(i);
    result = ['out_', i_s, '.mat'];
    load(result);
    
    subplot(3, 4, i);
    plot(out.t, x(:, i));
    hold on;
    plot(out.t, out.y_out);
    hold on;
    plot(out.t, x1(:, i));
    hold on;
    plot(out.t, x_moved(:, i));
    
    xlabel('t [s]');
    ylabel('x(t), y(t), x1(t), x(t) + 3');
    title(sprintf('Frequency = %.3f Hz', frequencies(i)));
    legend('x(t)', 'y(t)', 'x1(t)', 'x(t) + 3');
    grid on;

    leg = legend;
    leg.Position = [0.6, 0.1, 0.3, 0.2];
end



