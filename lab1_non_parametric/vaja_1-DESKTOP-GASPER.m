clear all
close all
%% PART 1
%% Nal 1
load MEASUREF.DAT
A = MEASUREF;
t = A(:,1);
u = A(:,2);
y = A(:,3);
n = A(:,4);

figure;
subplot(3,1,1);
stem(t, u);
xlabel('t');
ylabel('u');
title('Vhodni signal');

subplot(3,1,2);
stem(t, y);
xlabel('t');
ylabel('y');
title('Izhodni signal');

subplot(3,1,3);
stem(t, n);
xlabel('t');
ylabel('n');
title('Šum');

%% Nal 2
N = length(t);
T = t(2) - t(1);

F = 1/(N*T)
Fmax = 1/2*T

%% Nal 3
U = fft(u)/T; %%% zakaj deljeno s T???
%U1 = U(1);

% F1 = (0:N-1)/(N*T);

Y = fft(y)/T; %%% zakaj deljeno s T???
%Y1 = Y(1);
% 
% figure;
% subplot(2,1,1);
% stem(F1(1:N/2), abs(U(1:N/2)));
% xlabel('Frekvenca [Hz]');
% ylabel('|X(f)|');
% title('Absolutna vrednost Fourierove transformacije vhodnega signala');
% 
% subplot(2,1,2);
% stem(F1(1:N/2), abs(Y(1:N/2)));
% xlabel('Frekvenca [Hz]');
% ylabel('|Y(f)|');
% title('Absolutna vrednost Fourierove transformacije izhodnega signala');

stem(abs(U))

%% Nal 4
G = Y./U;
Gdb = 20*log10(abs(G));
faza = unwrap(angle(G))*180/pi; %Funkcija unwrap se uporablja za odstranjevanje skokov v faznem kotu

% Narišemo amplitudno karakteristiko v linearno-logaritemski skali
figure;
subplot(2,1,1);
semilogx(F1, Gdb);
xlabel('Frekvenca [Hz]');
ylabel('Amplituda [dB]');
title('Amplitudna karakteristika');

% Narišemo fazno karakteristiko v linearno-logaritemski skali
subplot(2,1,2);
semilogx(F1, faza);
xlabel('Frekvenca [Hz]');
ylabel('Fazni kot [°]');
title('Fazna karakteristika');

%% Nal 5
% Št. ničel enako št. polov vsaj 1 dominanten, mogoče 2, kateri je hiter
% Sistem 2 reda in nima konjugiranih kompleksnih polov

%% Nal 6
Ginv = ifft(G);
t = linspace(0, (N-1)*T, N); %vektor časovnih vzorcev

figure;
subplot(2,1,1);
plot(t, real(G));
xlabel('Čas (s)');
ylabel('Realni del impulznega odziva');

subplot(2,1,2);
plot(t, imag(G));
xlabel('Čas (s)');
ylabel('Imaginarni del impulznega odziva');

%% Nal 7 
noise = A(:,4); % noise is frequent
NOISE = fft(noise);
NOISEabs = abs(NOISE);

N0 = 0.25; % low-frequency amplitude
Ninf = 0.02; % high-frequency amplitude
w = linspace(0, (N-1)*T, N);
wg = 10;

Nw = Ninf + (N0 - Ninf)./sqrt(1 + (w/wg).^2);

figure;
plot(F1(1:N/2), NOISEabs(1:N/2), 'b');
hold on;
plot(F1(1:N/2), Nw(1:N/2), 'r');
xlabel('w (Hz)');
ylabel('N(w)');
legend('NOISEabs', 'Nw');

%% Nal 8
% the measured one:
stddev_of_error_absG_measured = NOISEabs ./ abs(U);

% the theoretical one:
stddev_of_error_absG_theoretical = Nw' ./ abs(U);

%% Nal 9
% a primer:
figure;
plot(F1(1:N/2), stddev_of_error_absG_measured(1:N/2), 'b');
hold on;
plot(F1(1:N/2), stddev_of_error_absG_theoretical(1:N/2), 'r');
xlabel('w (Hz)');
ylabel('standard deviation of the absolute error');
legend('measured', 'theoretical (model)');

%% b primer:
middle = G;
middle_db = 20*log10(abs(middle));
lower = middle - stddev_of_error_absG_theoretical;
lower_db = 20*log10(abs(lower));
upper = middle + stddev_of_error_absG_theoretical;
upper_db = 20*log10(abs(upper));

figure;
semilogx(F1(1:N/2), middle_db(1:N/2), 'r')
hold on
semilogx(F1(1:N/2), lower_db(1:N/2), 'b')
hold on
semilogx(F1(1:N/2), upper_db(1:N/2), 'g')
xlabel('w (Hz)')
ylabel('standard deviation of the absolute error')
legend('middle', 'lower', 'upper');

%% PART 2
clear all
close all

load MEASUREF.DAT
A = MEASUREF;
u = A(:,2); % given
y = procesf(u);

%% Nal 1
T = 1;
N = 256;
U = fft(u)/T;
Y = fft(y)/T;
G = Y./U;
Gdb = 20*log10(abs(G));
faza = unwrap(angle(G))*180/pi; %Funkcija unwrap se uporablja za odstranjevanje skokov v faznem kotu
F1 = (0:N-1)/(N*T);
%% Nal 2
% Narišemo amplitudno karakteristiko v linearno-logaritemski skali
figure;
subplot(2,1,1);
semilogx(F1, Gdb);
xlabel('Frekvenca [Hz]');
ylabel('Amplituda [dB]');
title('Amplitudna karakteristika');

% Narišemo fazno karakteristiko v linearno-logaritemski skali
subplot(2,1,2);
semilogx(F1, faza);
xlabel('Frekvenca [Hz]');
ylabel('Fazni kot [°]');
title('Fazna karakteristika');
