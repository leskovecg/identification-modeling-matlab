% IDENTIFICATION 
% GAŠPER LESKOVEC
% Identification and Modelling of Nonlinear Processes using Simulink S-Function: 
% A Comprehensive Parametric and Non-Parametric Approach with LS, RLS, and Instrumental Variables Methods 

% NOTES:
% sfun05 file, prbs file, simulink schemes has to be in THE SAME folder as main file !! 
%% ...
clear all;
close all;
%% static characteristic:
Ts = 0.01;
t = 0:Ts:10;
t = t';

Y = zeros(1,length(t))';

modelName = 'for_static_characteristic';
load_system(modelName);

U_static = zeros(1,length(-5:0.1:5))';
Y_static = zeros(1,length(-5:0.1:5))';
i = 1;

figure;
hold on;
xlabel('Time (s)');
ylabel('Output');
title('Time Series Plot');
grid on;

for U = -5:0.05:5
    U_static(i) = U;
    U = (U * ones(1,length(t)))';
    U = [t, U];
    simOut = sim(modelName);
    t_out = simOut.t.Time;
    y_out = simOut.Y.Data;
    plot(t_out, y_out);

    Y_static(i) = mean(y_out(500:end));
    
    i = i + 1;
end

hold off;

x = 0.5:0.05:1.5;
y = 0.1807*x - 0.0785;  

figure;
plot(U_static, Y_static, 'o-');
hold on;
plot(x, y, 'r--'); 
xlabel('U');
ylabel('Y');
title('static characteristic');
grid on;

% IMPORTANT INFO!!!
% x1 = 0.7   in x2 = 1.3     --> input
% y1 = 0.048 in y2 = 0.1564  --> output
%% dealing with prbs signal (preprocessing signals):
modelName = 'for_all';
load_system(modelName);
Ts_prbs = 0.1;
Ts = Ts_prbs;

% input signal prbs:
start_part = 0.7*ones(200,1);
middle_part = start_part(1) + 0.6*prbs(10);
ending_part = start_part(1) + 0.6*ones(200,1); 
u_prbs = [start_part; middle_part; ending_part];

% input signal prbs + time:
t_prbs = (0:length(u_prbs)-1)'*Ts_prbs;
u_matrix = [t_prbs, u_prbs];

% output signals of prbs (we override u_prbs and t_prbs in case of synchronization):
simOut_prbs = sim(modelName);
y_prbs = simOut_prbs.y.data;
u_prbs = simOut_prbs.u.data;
t_prbs = simOut_prbs.t.data;

% preprocessed input/output + time signals (cutting and subtracting mean value)
cutting = 195; 
u_prbs_preprocessed = u_prbs - mean(u_prbs(1:200));
u_prbs_preprocessed = u_prbs_preprocessed(cutting:end);
y_prbs_preprocessed = y_prbs - mean(y_prbs(50:175));
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

% calculate/plot the frequency response of the system:
T = t_prbs_preprocessed(2) - t_prbs_preprocessed(1);
N_prbs = length(u_prbs_preprocessed);
U_prbs = fft(u_prbs_preprocessed)/T;
Y_prbs = fft(y_prbs_preprocessed)/T;
G_prbs = Y_prbs./U_prbs;
Gdb_prbs = 20*log10(abs(G_prbs));

% phase in degrees:
phase_prbs = unwrap(angle(G_prbs))*180/pi;

% frequencies:
N = length(t_prbs_preprocessed);
freq_vector = (0:N-1)/(N*T);
freq_vector_halfed = freq_vector(1:N/2);

% plot the frequency response of the system:
figure;
subplot(2,1,1);
semilogx(freq_vector_halfed, Gdb_prbs(1:N/2));
xlabel('f [Hz]');
ylabel('|G(j2πf)| [dB]');
title('amplitude response');
grid on;

subplot(2,1,2);
semilogx(freq_vector_halfed, phase_prbs(1:N/2));
xlabel('f [Hz]');
ylabel('∠G(j2πf) [deg]');
title('phase response');
grid on;
%% a parametric identification --> LS 
% using prbs signal as verification signal here  
% validation signal will be chirp signal and it will be built in the next section 

Ts = 0.1;
d = 0;
n_s = 2;
t = t_prbs_preprocessed;
u = u_prbs_preprocessed;
y = y_prbs_preprocessed;
N = length(t);
T = t(2) - t(1);
f = ((0:N-1)/(N*T))';
U_for = u(n_s + 1:end);
Y_for = y(n_s + 1:end);

% construct the input matrix PSI based on the system order n:
PSI_for = zeros(length(Y_for), n_s*2);

for i = 1:n_s
    PSI_for(:,i) = -y(n_s-i+1:end-i);
    PSI_for(:,n_s+i) = u(n_s-i+1:end-i);
end

% add the delay to the input matrix:
if d > 0
    PSI_for = [zeros(d, n_s*2); PSI_for(1:end-d,:)];
elseif d < 0
    PSI_for = [PSI_for(-d+1:end,:); zeros(-d, n_s*2)];
end

% obtain the parameter estimates:
theta_for = PSI_for \ Y_for;
a = theta_for(1:n_s)';
b = theta_for(n_s+1:end)';

G = tf(b, [1 a], Ts);

[mag, phase] = bode(G, 2*pi*f);
mag = mag(:);
phase = phase(:);

mag_db = 20*log10(abs(mag));
phase_dgr = phase;

% comparison between y_sim and y:
y_sim = lsim(G, u, t); 

% frequencies:
freq_vector = (0:N-1)/(N*T); % od 0 do 1
freq_vector_halfed = freq_vector(1:N/2);

figure;
plot(t, y, 'r', 'LineWidth', 1.5); 
hold on;
plot(t, y_sim, 'b--', 'LineWidth', 1.5); 
hold off;
xlabel('t [s]');
ylabel('y(t), y_s_i_m(t)');
legend('measured output (PRBS)', 'simulated output (Parametric ID - LS)');
title('comparison of measured and simulated outputs using PRBS signal (varification)');

figure;
subplot(2,1,1);
semilogx(freq_vector_halfed, Gdb_prbs(1:N/2));
hold on; 
semilogx(freq_vector_halfed, mag_db(1:N/2));
xlabel('f [Hz]');
ylabel('|G(j2πf)| [dB]');
title('amplitude response comparison');
legend('amplitude response from PRBS (FFT)', 'amplitude response from Parametric Identification (LS)', 'Location', 'best');

subplot(2,1,2);
semilogx(freq_vector_halfed, phase_prbs(1:N/2));
hold on;
semilogx(freq_vector_halfed, phase_dgr(1:N/2));
xlabel('f [Hz]');
ylabel('∠G(j2πf) [°]');
title('phase response comparison');
legend('phase response from PRBS (FFT)', 'phase response from Parametric Identification (LS)', 'Location', 'best');

% standard deviation of the parameters:
residuals = Y_for - PSI_for * theta_for;
sigma2 = var(residuals); % calculate the variance of the residuals
cov_theta = inv(PSI_for' * PSI_for) * sigma2; % calculate the covariance matrix of the estimated parameters
var_theta = diag(cov_theta); % extract the variances of the estimated parameters (diagonal elements)
std_theta = sqrt(var_theta); % calculate the standard deviation of the estimated parameters
disp('standard deviation of the parameters (LS):');
disp(std_theta);
%% dealing with chirp signal (validation signal): 
% validation signal with low, middle, high frequencies

% dealing with validation signal (preprocessing signals):
modelName = 'for_all';
load_system(modelName);
Ts_validation = 0.1;
Ts = Ts_validation;

% define the time parameters - chirp signal:
sample_interval = 0.1;
signal_duration = 150;
constant_samples = 200;
constant_value = 0.7;
chirp_start_frequency = 0.01;
chirp_end_frequency = 0.6;
chirp_min_amplitude = 0.7;
chirp_max_amplitude = 1.3;

n_samples = signal_duration / sample_interval; % calculate the number of samples
constant_signal = constant_value * ones(1, constant_samples); % generate the constant part of the signal
t_chirp = linspace(0, signal_duration - (constant_samples * sample_interval), n_samples - constant_samples); % generate the time array for the chirp signal
chirp_signal = chirp(t_chirp, chirp_start_frequency, t_chirp(end), chirp_end_frequency, 'linear', 180); % generate the chirp signal with a slow increase starting from phase 180
chirp_signal = (chirp_max_amplitude - chirp_min_amplitude) * (chirp_signal - min(chirp_signal)) / (max(chirp_signal) - min(chirp_signal)) + chirp_min_amplitude; % scale the chirp signal
combined_signal = [constant_signal, chirp_signal]; % combine the constant and chirp signals

t_validation = [0:sample_interval:(n_samples-1)*sample_interval];

u_matrix = [t_validation', combined_signal'];

% output signals of chirp validation signal:
simOut_sinus = sim(modelName);
y_validation = simOut_sinus.y.data;
u_validation = simOut_sinus.u.data;
t_validation = simOut_sinus.t.data;

% preprocess signals:
cutting = 330; 
u_validation_preprocessed = u_validation - mean(u_validation(1:end));
u_validation_preprocessed = u_validation_preprocessed(cutting:end);
y_validation_preprocessed = y_validation - mean(y_validation(1:end));
y_validation_preprocessed = y_validation_preprocessed(cutting:end);
t_validation_preprocessed = t_validation(cutting:end);

figure;
subplot(2,2,1);
plot(t_validation, u_validation);
title('input signal');
xlabel('t [s]');
ylabel('u_c_h_i_r_p(t)');
grid on;
subplot(2,2,3);
plot(t_validation, y_validation);
title('output signal');
xlabel('t [s]');
ylabel('y_c_h_i_r_p(t)');
grid on;

subplot(2,2,2);
plot(t_validation_preprocessed, u_validation_preprocessed);
title('input preprocessed signal');
xlabel('t [s]');
ylabel('u_c_h_i_r_p_ _p_r_e_p_r_o_c_e_s_s_e_d(t)');
grid on;
subplot(2,2,4);
plot(t_validation_preprocessed, y_validation_preprocessed);
title('output preprocessed signal');
xlabel('t [s]');
ylabel('y_c_h_i_r_p_ _p_r_e_p_r_o_c_e_s_s_e_d(t)');
grid on;

% comparison between y_sim and y (validation):
y_sim_val = lsim(G, u_validation_preprocessed, t_validation_preprocessed); 

% visualization:
figure;
plot(t_validation_preprocessed, y_validation_preprocessed, 'r', 'LineWidth', 1.5); 
hold on;
plot(t_validation_preprocessed, y_sim_val, 'b--', 'LineWidth', 1.5); 
hold off;
xlabel('t [s]');
ylabel('y(t), y_s_i_m(t)');
legend('measured output (CHIRP)', 'simulated output (Parametric ID - LS)', 'Location', 'best');
title('comparison of measured and simulated outputs using CHIRP signal (validation)');

% the error of the obtained model should be estimated:
% compute the errors for training data:
error_train = y - y_sim;
rmse_train = sqrt(mean(error_train.^2));
mae_train = mean(abs(error_train));

% compute the errors for validation data:
error_val = y_validation_preprocessed - y_sim_val;
rmse_val = sqrt(mean(error_val.^2));
mae_val = mean(abs(error_val));

% display the results:
disp('estimated error of the obtained model (LS):');
fprintf('Training Data: RMSE = %.4f, MAE = %.4f\n', rmse_train, mae_train);
fprintf('Validation Data: RMSE = %.4f, MAE = %.4f\n', rmse_val, mae_val);
%% a non-parametric identification:

modelName = 'for_nonparametric_identification';
load_system(modelName);

% orthogonal correlation:
n=1501;
Ts = 0.1;
num_freqs = 11; 
frequencies = logspace(log10(0.01), log10(4.5), num_freqs); % from 0.01 Hz to 4.5 Hz
ampl = zeros(1, num_freqs);
faza = nan(1, num_freqs);
t_nonparametric = (0:1500-1)'*Ts;

for i = 1:num_freqs

    x_nonparametric =  1 + 0.3 *cos(2 * pi * frequencies(i) * t_nonparametric);
    x_nonparametric_cos =  cos(2 * pi * frequencies(i) * t_nonparametric);
    x1_nonparametric = sin(2 * pi * frequencies(i) * t_nonparametric);
    
    u_matrix_cos = [t_nonparametric, x_nonparametric_cos];
    u1_matrix = [t_nonparametric, x1_nonparametric];
    u_matrix = [t_nonparametric, x_nonparametric];

    % output signals:
    simOut_nonparametric = sim(modelName);
    y_nonparametric = simOut_nonparametric.y.data(1:end);
    t_nonparametric = simOut_nonparametric.t.data(1:end);
    x_nonparametric = simOut_nonparametric.u_cos.data(1:end);
    x1_nonparametric = simOut_nonparametric.u1.data(1:end);

    % preprocess:
    whole_time = 150;
    cutted_time = 30;
    time_left = whole_time - cutted_time;
    period = 1 / frequencies(i);
    n_periods = floor(time_left / period);
    time_of_periods = n_periods * period;
    samples = time_of_periods / Ts;
    start = floor(length(t_nonparametric) - samples);

    t_nonparametric_preprocessed = t_nonparametric(start:end); 
    x1_nonparametric_preprocessed = x1_nonparametric(start:end);
    x_nonparametric_preprocessed = x_nonparametric(start:end);
    y_nonparametric_preprocessed = y_nonparametric(start:end);

    value = mean(y_nonparametric_preprocessed);
    y_nonparametric_preprocessed = y_nonparametric_preprocessed - value;

    re = y_nonparametric_preprocessed'*x_nonparametric_preprocessed*2/(length(y_nonparametric_preprocessed)*0.3);
    im = -y_nonparametric_preprocessed'*x1_nonparametric_preprocessed*2/(length(y_nonparametric_preprocessed)*0.3);
    ampl(i) = sqrt(re^2+im^2);
    faza(i) = atan2(im,re);

    title_str = sprintf('frequency = %.3f Hz', frequencies(i)); 
    
    figure;
    plot(t_nonparametric_preprocessed, x_nonparametric_preprocessed);
    hold on;
    plot(t_nonparametric_preprocessed, x1_nonparametric_preprocessed);
    hold on;
    plot(t_nonparametric_preprocessed, y_nonparametric_preprocessed);
    xlabel('t [s]');
    ylabel('x(t), x1(t), y(t)');
    title(title_str);
    legend('x(t) = cos(2πft)', 'x1(t) = sin(2πft)', 'y(t)', 'Location', 'best');
end

ampl = 20*log10(abs(ampl));

% phase unwraping:
faza=(unwrap(faza))*180/pi;

% frequencies:
N = length(t_nonparametric);
T = t_nonparametric(2) - t_nonparametric(1);
freq_vector = (0:N-1)/(N*T); 
freq_vector_halfed = freq_vector(1:N/2);

% visualization: 
figure;
subplot(2, 1, 1);
semilogx(frequencies, ampl, '.k');
hold on;
semilogx(freq_vector_halfed, Gdb_prbs(1:N/2));
title('amplitude response (orthogonal correlation method)');
xlabel('f [Hz]');
ylabel('|G(j2πf)| [dB]');
legend('identified','true');

subplot(2, 1, 2);
semilogx(frequencies, faza, '.k');
hold on;
semilogx(freq_vector_halfed, phase_prbs(1:N/2));
title('phase response (orthogonal correlation method)');
xlabel('f [Hz]');
ylabel('∠(G(j2πf)) [°]');
legend('identified', 'true');
%% noise:
Ts = 0.1;
modelName = 'for_all';
load_system(modelName);

u_noise = 1*ones(length(t_prbs), 1);
u_matrix = [t_prbs, u_noise];

simOut_noise = sim(modelName);
y_noise = simOut_noise.y.data;
t_noise = simOut_noise.t.data;

cut = 195;
t_noise = t_noise(cut:length(t_prbs));
y_noise = y_noise(cut:length(t_prbs));
y_noise = y_noise - mean(y_noise);
u_noise = u_noise(cut:length(t_prbs));

N_noise = length(t_noise);
T_noise = t_noise(2) - t_noise(1);
freq_vector_noise = (0:N_noise-1)/(N_noise*T_noise);
freq_vector_halfed_noise = freq_vector_noise(1:N_noise/2);
Y_noise = fft(y_noise)/T_noise;
Y_noise = abs(Y_noise);

% noise curve parameters:
N0 = 4.5; % low-frequency amplitude
Ninf = 0.3; % high-frequency amplitude
w = linspace(0, (N_noise-1)*T_noise, N_noise);
wg = 2.5;

% model of the colour noise:
Nw = Ninf + (N0 - Ninf)./sqrt(1 + (w/wg).^2);

figure;
subplot(3,1,1);
plot(t_noise, u_noise);
xlabel('t [s]');
ylabel('c(t)');
title('input signal (constant value)');

subplot(3,1,2);
plot(t_noise, y_noise);
xlabel('t [s]');
ylabel('n(t)');
title('output signal (noise)');

subplot(3,1,3);
plot(freq_vector_halfed_noise, Y_noise(1:N_noise/2));
hold on;
plot(freq_vector_halfed_noise, Nw(1:N_noise/2), 'r');
xlabel('f [Hz]');
ylabel('N(w)');
title('comparison between fft of the noise signal and model of the noise');
legend('noise FFT', 'noise model');
%% dealing with standard deviation:

% the measured one:
stddev_of_error_absG_measured = Y_noise ./ abs(Y_prbs);

% the theoretical one:
stddev_of_error_absG_theoretical = Nw' ./ abs(Y_prbs);

% a primer:
figure;
plot(freq_vector_halfed_noise, stddev_of_error_absG_measured(1:N_noise/2), 'b');
hold on;
plot(freq_vector_halfed_noise, stddev_of_error_absG_theoretical(1:N_noise/2), 'r');
xlabel('w (Hz)');
ylabel('standard deviation of the absolute error');
legend('measured', 'theoretical (model)');

% b primer:
middle_G_prbs = G_prbs;
middle_db = 20*log10(abs(middle_G_prbs));

stddev_of_error_absG_theoretical = 20*log10(stddev_of_error_absG_theoretical);

lower_db =  middle_db + stddev_of_error_absG_theoretical;

upper_db = middle_db - stddev_of_error_absG_theoretical;

figure;
semilogx(freq_vector_halfed_noise, middle_db(1:N_noise/2), 'r')
hold on
semilogx(freq_vector_halfed_noise, lower_db(1:N_noise/2), 'b')
hold on
semilogx(freq_vector_halfed_noise, upper_db(1:N_noise/2), 'g')
xlabel('w (Hz)')
ylabel('standard deviation of the absolute error')
legend('middle', 'lower', 'upper');
%% an arbitrary parametric identification method --> RLS:

Ts = 0.1;
d = 0;
n_s = 2;
t = t_prbs_preprocessed;
u = u_prbs_preprocessed;
y = y_prbs_preprocessed;
N = length(t);
T = t(2) - t(1);
f = ((0:N-1)/(N*T))';

U_for = u(n_s + 1:end);
Y_for = y(n_s + 1:end);

% Construct the input matrix PSI based on the system order n
PSI_for = zeros(length(Y_for), n_s*2);

for i = 1:n_s
    PSI_for(:,i) = -y(n_s-i+1:end-i);
    PSI_for(:,n_s+i) = u(n_s-i+1:end-i);
end

% Add the delay to the input matrix
if d > 0
    PSI_for = [zeros(d, n_s*2); PSI_for(1:end-d,:)];
elseif d < 0
    PSI_for = [PSI_for(-d+1:end,:); zeros(-d, n_s*2)];
end

% RLS method
lambda = 0.99; % forgetting factor
P = 1e6 * eye(n_s*2);
theta_for = zeros(n_s*2, 1);

theta_for_all = zeros(n_s*2, size(PSI_for, 1));

% create a figure for live parameter updates:
figure;
subplot(2, 2, 1);
title('a1');
subplot(2, 2, 2);
title('a2');
subplot(2, 2, 3);
title('b1');
subplot(2, 2, 4);
title('b2');

for k = 1:size(PSI_for, 1)
    
    phi = PSI_for(k, :)';
    K = P * phi / (lambda + phi' * P * phi);
    P = (P - K * phi' * P) / lambda;
    e = Y_for(k) - phi' * theta_for;
    theta_for = theta_for + K * e;
    theta_for_all(:, k) = theta_for;

    % update the live parameter plots:
    subplot(2, 2, 1);
    plot(theta_for_all(1, 1:k), 'r');
    xlabel('Iteration');
    ylabel('a1');

    subplot(2, 2, 2);
    plot(theta_for_all(2, 1:k), 'g');
    xlabel('Iteration');
    ylabel('a2');

    subplot(2, 2, 3);
    plot(theta_for_all(3, 1:k), 'b');
    xlabel('Iteration');
    ylabel('b1');

    subplot(2, 2, 4);
    plot(theta_for_all(4, 1:k), 'm');
    xlabel('Iteration');
    ylabel('b2');
 
    % pause(0.001);
end

%%%
std_a = std(theta_for_all(1:n_s, :), 0, 2);
std_b = std(theta_for_all(n_s+1:end, :), 0, 2);

% Display the results
disp('standard deviation of the parameters (RLS):');
fprintf('Standard Deviation of a: %s\n', mat2str(std_a));
fprintf('Standard Deviation of b: %s\n', mat2str(std_b));
%%%

a = theta_for(1:n_s)';
b = theta_for(n_s+1:end)';

G = tf(b, [1 a], Ts);

[mag(:), phase(:)] = bode(G, 2*pi*f);

mag_db = 20*log10(abs(mag));
phase_dgr = phase;

% comparison between y_sim and y - varification:
y_sim = lsim(G, u, t); 

% validation:
y_sim_val = lsim(G, u_validation_preprocessed, t_validation_preprocessed);

figure;
plot(t, y, 'r', 'LineWidth', 1.5); 
hold on;
plot(t, y_sim, 'b--', 'LineWidth', 1.5); 
hold off;
xlabel('t [s]');
ylabel('y(t), y_s_i_m(t)');
legend('measured output (PRBS)', 'simulated output (Parametric ID - RLS)', 'Location', 'best');
title('comparison of measured and simulated outputs using PRBS signal (varification)');

figure;
plot(t_validation_preprocessed, y_validation_preprocessed, 'r', 'LineWidth', 1.5); 
hold on;
plot(t_validation_preprocessed, y_sim_val, 'b--', 'LineWidth', 1.5); 
hold off;
xlabel('t [s]');
ylabel('y(t), y_s_i_m(t)');
legend('measured output (CHIRP)', 'simulated output (Parametric ID - RLS)', 'Location', 'best');
title('comparison of measured and simulated outputs using CHIRP signal (validation)');

% frequencies:
freq_vector = (0:N-1)/(N*T); % od 0 do 1
freq_vector_halfed = freq_vector(1:N/2);

figure;
subplot(2,1,1);
semilogx(freq_vector_halfed, Gdb_prbs(1:N/2));
hold on; 
semilogx(freq_vector_halfed, mag_db(1:N/2));
xlabel('f [Hz]');
ylabel('|G(j2πf)| [dB]');
title('amplitude response comparison');
legend('amplitude response from PRBS (FFT)', 'amplitude response from Parametric Identification (RLS)', 'Location', 'best');

subplot(2,1,2);
semilogx(freq_vector_halfed, phase_prbs(1:N/2));
hold on;
semilogx(freq_vector_halfed, phase_dgr(1:N/2));
xlabel('f [Hz]');
ylabel('∠(G(j2πf)) [°]');
title('phase response comparison');
legend('phase response from PRBS (FFT)', 'phase response from Parametric Identification (RLS)', 'Location', 'best');

% the error of the obtained model should be estimated:
% compute the errors for training data:
error_train = y - y_sim;
rmse_train = sqrt(mean(error_train.^2));
mae_train = mean(abs(error_train));

% compute the errors for validation data:
error_val = y_validation_preprocessed - y_sim_val;
rmse_val = sqrt(mean(error_val.^2));
mae_val = mean(abs(error_val));

% display the results:
disp('estimated error of the obtained model (RLS):');
fprintf('Training Data: RMSE = %.4f, MAE = %.4f\n', rmse_train, mae_train);
fprintf('Validation Data: RMSE = %.4f, MAE = %.4f\n', rmse_val, mae_val);
%% an arbitrary parametric identification method --> instrumental variables:

Ts = 0.1;
d = 0;
n_s = 2;
t = t_prbs_preprocessed;
u = u_prbs_preprocessed;
y = y_prbs_preprocessed;
N = length(t);
T = t(2) - t(1);
f = ((0:N-1)/(N*T))';
U_for = u(n_s + 1:end);
Y_for = y(n_s + 1:end);

% Construct the input matrix PSI based on the system order n
PSI_for = zeros(length(Y_for), n_s*2);

for i = 1:n_s
    PSI_for(:,i) = -y(n_s-i+1:end-i);
    PSI_for(:,n_s+i) = u(n_s-i+1:end-i);
end

% Add the delay to the input matrix
if d > 0
    PSI_for = [zeros(d, n_s*2); PSI_for(1:end-d,:)];
elseif d < 0
    PSI_for = [PSI_for(-d+1:end,:); zeros(-d, n_s*2)];
end

% Obtain the parameter estimates
theta_for = PSI_for \ Y_for;
a = theta_for(1:n_s)';
b = theta_for(n_s+1:end)';

% transfer function
G = tf(b, [1 a], Ts);

% bode diagram
[mag, phase] = bode(G, 2*pi*f);

mag = mag(:);
phase = phase(:);

mag_db = 20*log10(abs(mag));
phase_dgr = phase;

% calculate simulated response y_sim (varification):
y_sim_G = lsim(G, u, t);

% validation:
y_sim_val_G = lsim(G, u_validation_preprocessed, t_validation_preprocessed);

figure;
plot(t, y, 'r', 'LineWidth', 1.5); 
hold on;
plot(t, y_sim_G, 'b--', 'LineWidth', 1.5); 
hold off;
xlabel('t [s]');
ylabel('y(t), y_s_i_m(t)');
legend('measured output (PRBS)', 'simulated output (Parametric ID - IV)');
title('comparison of measured and simulated outputs using PRBS signal (varification), G', 'FontSize', 8);

figure;
plot(t_validation_preprocessed, y_validation_preprocessed, 'r', 'LineWidth', 1.5); 
hold on;
plot(t_validation_preprocessed, y_sim_val_G, 'b--', 'LineWidth', 1.5); 
hold off;
xlabel('t [s]');
ylabel('y(t), y_s_i_m(t)');
legend('measured output (CHIRP)', 'simulated output (Parametric ID - IV)');
title('comparison of measured and simulated outputs using CHIRP signal (validation), G', 'FontSize', 8);

% construct W matrix:
W = zeros(length(y_sim_G(n_s+1:end)), n_s*2);

for i = 1:n_s
    W(:,i) = -y_sim_G(n_s-i+1:end-i);
    W(:,n_s+i) = u(n_s-i+1:end-i);
end

% Instrumental Variables estimation with W matrix
theta_iv_W = (W' * PSI_for) \ (W' * Y_for);
a_W = theta_iv_W(1:n_s)';
b_W = theta_iv_W(n_s+1:end)';

% Transfer function estimation with W matrix
G_W = tf(b_W, [1 a_W], Ts);

[mag_W(:), phase_W(:)] = bode(G_W, 2*pi*f);

mag_db_W = 20*log10(abs(mag_W));
phase_dgr_W = phase_W;

% calculate simulated response y_sim (varification):
y_sim_G_W = lsim(G_W, u, t);

% validation:
y_sim_val_G_W = lsim(G_W, u_validation_preprocessed, t_validation_preprocessed);

figure;
plot(t, y, 'r', 'LineWidth', 1.5); 
hold on;
plot(t, y_sim_G_W, 'b--', 'LineWidth', 1.5); 
hold off;
xlabel('t [s]');
ylabel('y(t), y_s_i_m(t)');
legend('measured output (PRBS)', 'simulated output (Parametric ID - IV)');
title('comparison of measured and simulated outputs using PRBS signal (varification), G_W', 'FontSize', 8);

figure;
plot(t_validation_preprocessed, y_validation_preprocessed, 'r', 'LineWidth', 1.5); 
hold on;
plot(t_validation_preprocessed, y_sim_val_G_W, 'b--', 'LineWidth', 1.5); 
hold off;
xlabel('t [s]');
ylabel('y(t), y_s_i_m(t)');
legend('measured output (CHIRP)', 'simulated output (Parametric ID - IV)');
title('comparison of measured and simulated outputs using CHIRP signal (validation), G_W', 'FontSize', 8);

% frequencies:
freq_vector = (0:N-1)/(N*T); % od 0 do 1
freq_vector_halfed = freq_vector(1:N/2);

figure;
subplot(2,1,1);
semilogx(freq_vector_halfed, Gdb_prbs(1:N/2));
hold on; 
semilogx(freq_vector_halfed, mag_db(1:N/2));
hold on; 
semilogx(freq_vector_halfed, mag_db_W(1:N/2));
xlabel('f [Hz]');
ylabel('|G(j2πf)| [dB]');
title('amplitude response comparison');
legend('amplitude response from PRBS (FFT)', 'amplitude response from Parametric Identification (IV), G', 'amplitude response from Parametric Identification (IV), G_W', 'Location', 'best');

subplot(2,1,2);
semilogx(freq_vector_halfed, phase_prbs(1:N/2));
hold on;
semilogx(freq_vector_halfed, phase_dgr(1:N/2));
hold on; 
semilogx(freq_vector_halfed, phase_dgr_W(1:N/2));
xlabel('f [Hz]');
ylabel('∠(G(j2πf)) [°]');
title('phase response comparison');
legend('phase response from PRBS (FFT)', 'phase response from Parametric Identification (IV), G', 'phase response from Parametric Identification (IV), G_W', 'Location', 'best');

% standard deviation of the parameters:
residuals = Y_for - W * theta_iv_W;
sigma2 = var(residuals); % calculate the variance of the residuals
cov_theta = inv(W' * W) * sigma2; % calculate the covariance matrix of the estimated parameters
var_theta = diag(cov_theta); % extract the variances of the estimated parameters (diagonal elements)
std_theta = sqrt(var_theta); % calculate the standard deviation of the estimated parameters
disp('standard deviation of the parameters (IV):');
disp(std_theta);

% the error of the obtained model should be estimated:
% compute the errors for training data:
error_train = y - y_sim_G;
rmse_train = sqrt(mean(error_train.^2));
mae_train = mean(abs(error_train));

% compute the errors for validation data:
error_val = y_validation_preprocessed - y_sim_val_G;
rmse_val = sqrt(mean(error_val.^2));
mae_val = mean(abs(error_val));

% display the results:
disp('estimated error of the obtained model (IV), G:');
fprintf('Training Data: RMSE = %.4f, MAE = %.4f\n', rmse_train, mae_train);
fprintf('Validation Data: RMSE = %.4f, MAE = %.4f\n', rmse_val, mae_val);

% the error of the obtained model should be estimated:
% compute the errors for training data:
error_train = y - y_sim_G_W;
rmse_train = sqrt(mean(error_train.^2));
mae_train = mean(abs(error_train));

% compute the errors for validation data:
error_val = y_validation_preprocessed - y_sim_val_G_W;
rmse_val = sqrt(mean(error_val.^2));
mae_val = mean(abs(error_val));

% display the results:
disp('estimated error of the obtained model (IV), G_W:');
fprintf('Training Data: RMSE = %.4f, MAE = %.4f\n', rmse_train, mae_train);
fprintf('Validation Data: RMSE = %.4f, MAE = %.4f\n', rmse_val, mae_val);














