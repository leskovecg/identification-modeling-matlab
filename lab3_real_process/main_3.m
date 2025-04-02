
% addpath('C:\Users\Uporabnik\OneDrive\FAKS\MAGISTERIJ\letnik_1\SEMESTER_2\MODUL_A_IDENTIFIKACIJA\lab_vaje\vaja_3\results\')
% addpath('C:\Users\Uporabnik\OneDrive\FAKS\MAGISTERIJ\letnik_1\SEMESTER_2\MODUL_A_IDENTIFIKACIJA\lab_vaje\vaja_3\functions')
% addpath('C:\Users\Uporabnik\OneDrive\FAKS\MAGISTERIJ\letnik_1\SEMESTER_2\MODUL_A_IDENTIFIKACIJA\lab_vaje\vaja_3\data')
% addpath('C:\Users\Uporabnik\OneDrive\FAKS\MAGISTERIJ\letnik_1\SEMESTER_2\MODUL_A_IDENTIFIKACIJA\lab_vaje\vaja_3\simulink')
%% clearing
clear all;
close all;
%% 1. part 
%% static characteristic Ts=0.1
U_static = -5:1:5;
Y_static = [-5; -4; -2.7; -1.4; -0.2; 0; 0.7; 1.8; 2.9; 3.9; 5.1];

figure;
plot(U_static, Y_static, 'o-');
xlabel('U');
ylabel('Y');
title('static characteristic');
grid on;
%% static characteristic Ts=0.05
U_static = 2:1:5;
Y_static = [1.6; 2.7; 3.8; 5.6];

figure;
plot(U_static, Y_static, 'o-');
xlabel('U');
ylabel('Y');
title('static characteristic');
grid on;
%% dealing with prbs signal (preprocessing signals):
Ts_prbs = 0.1;

% input signal prbs:
start = 2.7*ones(200,1);
middle = start(1) + prbs(10);
ending = start(1) + ones(200,1); 
u_prbs = [start; middle; ending];

% input signal prbs + time:
t_prbs = (0:length(u_prbs)-1)'*Ts_prbs;
u_prbs_matrix = [t_prbs, u_prbs];

% load the data:
simOutEna = load('simOutEna.mat');
y_prbs = simOutEna.simOut_ena.Y(1:end-1); 

% preprocessing signals (input, output, time):
% a) subtracting operating point
% b) deleting transient
cutting = 200;
u_prbs_preprocessed = u_prbs - mean(u_prbs(1:200));
u_prbs_preprocessed = u_prbs_preprocessed(cutting:end);
y_prbs_preprocessed = y_prbs - mean(y_prbs(20:200));
y_prbs_preprocessed = y_prbs_preprocessed(cutting:end);
t_prbs_preprocessed = t_prbs(cutting:end);

% visualize input prbs / output signal (before and after preprocessing):
figure;
subplot(2,2,1);
plot(t_prbs, u_prbs);
title('input signal');
xlabel('t [s]');
ylabel('u_p_r_b_s(t)');
grid on;
subplot(2,2,3);
plot(t_prbs, y_prbs);
title('output signal');
xlabel('t [s]');
ylabel('y_p_r_b_s(t)');
grid on;

subplot(2,2,2);
plot(t_prbs_preprocessed, u_prbs_preprocessed);
title('input preprocessed signal');
xlabel('t [s]');
ylabel('u_p_r_b_s_ _p_r_e_p_r_o_c_e_s_s_e_d(t)');
grid on;
subplot(2,2,4);
plot(t_prbs_preprocessed, y_prbs_preprocessed);
title('output preprocessed signal');
xlabel('t [s]');
ylabel('y_p_r_b_s_ _p_r_e_p_r_o_c_e_s_s_e_d(t)');
grid on;
%% calculate/plot the frequency response of the system:
N_prbs = length(u_prbs_preprocessed);
T_prbs = t_prbs_preprocessed(2) - t_prbs_preprocessed(1);
U_prbs = fft(u_prbs_preprocessed)/T_prbs;
Y_prbs = fft(y_prbs_preprocessed)/T_prbs;
G_prbs = Y_prbs./U_prbs;
Gdb_prbs = 20*log10(abs(G_prbs));

% phase in degrees:
phase_prbs = unwrap(angle(G_prbs))*180/pi;

% frequencies:
f_prbs = (0:N_prbs-1)/(N_prbs*Ts_prbs);

% plot the frequency response of the system:
figure;
subplot(2,1,1);
semilogx(f_prbs(1:613), Gdb_prbs(1:613));
xlabel('f [Hz]');
ylabel('|G(j2πf)| [dB]');
title('magnitude response');
grid on;

subplot(2,1,2);
semilogx(f_prbs(1:613), phase_prbs(1:613));
xlabel('f [Hz]');
ylabel('∠G(j2πf) [deg]');
title('phase response');
grid on;
%% dealing with noise signal:
load('noise.mat');
u_noise = ones(1201 + 23,1) * 3; % constant input
t_noise = out.t(301 - 23:end, 1); % time
y_noise = out.n(301 -23:end, 1); % output
y_noise = y_noise - mean(y_noise);

N_noise = length(t_noise);
T_noise = t_noise(2) - t_noise(1);
freq_vector = (0:N_noise-1)/(N_noise*T_noise);
freq_vector_halfed = freq_vector(1:N_noise/2);

% noise curve parameters:
N0 = 60; % low-frequency amplitude
Ninf = 5; % high-frequency amplitude
w = linspace(0, (N_noise-1)*T_noise, N_noise);
wg = 1;

% model of the colour noise:
Nw = Ninf + (N0 - Ninf)./sqrt(1 + (w/wg).^2);

% fft of the noise:
Y_noise = fft(y_noise)/T_noise;
Y_noise = abs(Y_noise);

% visualization of the noise and its specter:
figure;
subplot(3,1,1);
plot(t_noise, u_noise); 
xlabel('t [s]');
ylabel('c(t)');
title('input signal (constant value)');
grid on;

subplot(3,1,2);
plot(t_noise, y_noise); 
xlabel('t [s]');
ylabel('n(t)');
title('output signal (noise)');
grid on;

subplot(3,1,3); 
plot(freq_vector_halfed, Y_noise(1:N_noise/2));
hold on;
plot(freq_vector_halfed, Nw(1:N_noise/2), 'r');
xlabel('w [Hz]');
ylabel('Y_n_o_i_s_e(w), N(w)');
title('comparison between fft of the noise signal and model of the noise');
legend('noise FFT', 'noise model');
grid on;
%% dealing with standard deviation:
N = N_noise;
% the measured one:
stddev_of_error_absG_measured = Y_noise ./ abs(Y_prbs);

% the theoretical one:
stddev_of_error_absG_theoretical = Nw' ./ abs(Y_prbs);

% a primer:
figure;
plot(freq_vector_halfed, stddev_of_error_absG_measured(1:N/2), 'b');
hold on;
plot(freq_vector_halfed, stddev_of_error_absG_theoretical(1:N/2), 'r');
xlabel('w [Hz]');
ylabel('standard deviation of the absolute error');
legend('measured', 'theoretical (model)');

% b primer:
middle = G_prbs;
middle_db = 20*log10(abs(middle));

stddev_of_error_absG_theoretical = 20*log10(abs(stddev_of_error_absG_theoretical));

lower = middle_db + stddev_of_error_absG_theoretical;
% lower_db = 20*log10(abs(lower));

upper = middle_db - stddev_of_error_absG_theoretical;
% upper_db = 20*log10(abs(upper));

figure;
semilogx(freq_vector_halfed, middle_db(1:N/2), 'r')
hold on
semilogx(freq_vector_halfed, lower(1:N/2), 'b')
hold on
semilogx(freq_vector_halfed, upper(1:N/2), 'g')
xlabel('w (Hz)')
ylabel('standard deviation of the absolute error')
legend('middle', 'lower', 'upper');
%% 2. del --> NONPARAMETRIC 

% orthogonal correlation:
n=1501;
Ts = 0.1;
num_freqs = 11; 
frequencies = logspace(log10(0.05), log10(0.2), num_freqs); % from 0.05 Hz to 0.2 Hz
ampl = nan(1,num_freqs); 
faza = nan(1,num_freqs);

for i = 1:num_freqs
    
    i_s = num2str(i);
    result = ['out_', i_s, '.mat'];
    load(result);

    whole_time = 150;
    cutted_time = 30;
    time_left = whole_time - cutted_time;
    period = 1 / frequencies(i);
    n_periods = floor(time_left / period);
    time_of_periods = n_periods * period;
    samples = time_of_periods / Ts;
    start = floor(length(out.t) - samples);

    t = out.t(start:end, 1);
    x1 = sin(2*pi*frequencies(i)*t);
    x = cos(2*pi*frequencies(i)*t);
    y = out.y_out(start:end, 1);

    re = y'*x*2/(length(y) * 0.5);
    im = -y'*x1*2/(length(y) * 0.5);
    ampl(i) = sqrt(re^2+im^2);
    faza(i) = atan2(im,re);
    
    figure;
    plot(t, x1);
    hold on;
    plot(t, x);
    hold on;
    plot(t, y);
end

figure;
subplot(2, 1, 1);
semilogx(frequencies, ampl, '.k');
hold on;
semilogx(f_prbs(1:613), Gdb_prbs(1:613));
title('amplitude response (orthogonal correlation method)');
xlabel('f [Hz]');
ylabel('|G(j2πf)| [dB]');
legend('identified','true');

subplot(2, 1, 2);
semilogx(frequencies, faza*180/pi, '.k');
hold on;
semilogx(f_prbs(1:613), phase_prbs(1:613));
title('phase response (orthogonal correlation method)');
xlabel('f [Hz]');
ylabel('∠(G(j2πf)) [°]');
legend('identified', 'true');

