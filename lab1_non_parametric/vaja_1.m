clear all
close all
%% PART 1
%% Nal 1
load MEASUREF.DAT
t = MEASUREF(:,1);
u = MEASUREF(:,2);
y = MEASUREF(:,3);
n = MEASUREF(:,4);

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

F = 1/(N*T);
fmax = 1/2*T;

%% Nal 3
U = fft(u)/T; % zakaj delimo s T????
Y = fft(y)/T;
freq_vector = (0:N-1)/(N*T); % od 0 do 1

% upoštevamo samo prvi del spektra (mirrored frequencies)
freq_vector_halfed = freq_vector(1:N/2);
U_halfed = U(1:N/2);
Y_halfed = Y(1:N/2);

figure;
subplot(2,1,1);
stem(freq_vector_halfed, abs(U_halfed));
xlabel('Frekvenca [Hz]');
ylabel('|U(f)|');
title('Absolutna vrednost Fourierove transformacije vhodnega signala');

subplot(2,1,2);
stem(freq_vector_halfed, abs(Y_halfed));
xlabel('Frekvenca [Hz]');
ylabel('|Y(f)|');
title('Absolutna vrednost Fourierove transformacije izhodnega signala');

%% Nal 4
G = Y./U;
Glin = abs(G);
Gdb = 20*log10(abs(G));
Gphase = unwrap(angle(G))*180/pi; %Funkcija unwrap se uporablja za odstranjevanje skokov v faznem kotu

% upoštevamo samo prvi del spektra (mirrored frequencies)
Glin_halfed = Glin(1:N/2);
Gdb_halfed = Gdb(1:N/2);
Gphase_halfed = Gphase(1:N/2);

%%% upoštevana oba dela spektra:
% Narišemo amplitudno karakteristiko v linearno-linearni skali
figure;
subplot(3,1,1);
plot(freq_vector, Glin);
xlabel('Frekvenca [Hz]');
ylabel('Amplituda');
title('Amplitudna karakteristika');

% Narišemo amplitudno karakteristiko v linearno-logaritemski skali
subplot(3,1,2);
semilogx(freq_vector, Gdb);
xlabel('Frekvenca [Hz]');
ylabel('Amplituda [dB]');
title('Amplitudna karakteristika');

% Narišemo fazno karakteristiko v linearno-logaritemski skali
subplot(3,1,3);
semilogx(freq_vector, Gphase);
xlabel('Frekvenca [Hz]');
ylabel('Fazni kot [°]');
title('Fazna karakteristika');

%%% upoštevan prvi del spektra:
% Narišemo amplitudno karakteristiko v linearno-linearni skali
figure;
subplot(3,1,1);
plot(freq_vector_halfed, Glin_halfed);
xlabel('Frekvenca [Hz]');
ylabel('Amplituda');
title('Amplitudna karakteristika');

% Narišemo amplitudno karakteristiko v linearno-logaritemski skali
subplot(3,1,2);
semilogx(freq_vector_halfed, Gdb_halfed);
xlabel('Frekvenca [Hz]');
ylabel('Amplituda [dB]');
title('Amplitudna karakteristika');

% Narišemo fazno karakteristiko v linearno-logaritemski skali
subplot(3,1,3);
semilogx(freq_vector_halfed, Gphase_halfed);
xlabel('Frekvenca [Hz]');
ylabel('Fazni kot [°]');
title('Fazna karakteristika');

%% Nal 5
% Št. ničel enako št. polov vsaj 1 dominanten, mogoče 2, kateri je hiter
% Sistem 2 reda in nima konjugiranih kompleksnih polov

%% Nal 6
ginv = ifft(G);
t = linspace(0, (N-1)*T, N); %vektor časovnih vzorcev

y_check = ginv .* u; %% tabi mogla bit enaka z y_check == y pa ni ???? šum itak ni dodan
nois = y - y_check;  


figure;
subplot(5,1,1);
plot(t, ginv);
xlabel('Čas (s)');
ylabel('Impulse response');

subplot(5,1,2);
plot(t, u);
xlabel('Čas (s)');
ylabel('u');

subplot(5,1,3);
plot(t, y_check);
xlabel('Čas (s)');
ylabel('y_check');

subplot(5,1,4);
plot(t, y);
xlabel('Čas (s)');
ylabel('y');

subplot(5,1,5);
plot(t, nois);
xlabel('Čas (s)');
ylabel('nois');


figure;
subplot(2,1,1);
% plot(t, real(G));
plot(t, real(ginv));
xlabel('Čas (s)');
ylabel('Realni del impulznega odziva');

subplot(2,1,2);
% plot(t, imag(G));
plot(t, imag(ginv));
xlabel('Čas (s)');
ylabel('Imaginarni del impulznega odziva');

%% Nal 7 
NOISE = fft(n); % frequent noise
NOISEabs = abs(NOISE);

N0 = 0.25; % low-frequency amplitude
Ninf = 0.02; % high-frequency amplitude
w = linspace(0, (N-1)*T, N);
wg = 10;

Nw = Ninf + (N0 - Ninf)./sqrt(1 + (w/wg).^2);

figure;
plot(freq_vector_halfed, NOISEabs(1:N/2), 'b');
hold on;
plot(freq_vector_halfed, Nw(1:N/2), 'r');
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
plot(freq_vector_halfed, stddev_of_error_absG_measured(1:N/2), 'b');
hold on;
plot(freq_vector_halfed, stddev_of_error_absG_theoretical(1:N/2), 'r');
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
semilogx(freq_vector_halfed, middle_db(1:N/2), 'r')
hold on
semilogx(freq_vector_halfed, lower_db(1:N/2), 'b')
hold on
semilogx(freq_vector_halfed, upper_db(1:N/2), 'g')
xlabel('w (Hz)')
ylabel('standard deviation of the absolute error')
legend('middle', 'lower', 'upper');