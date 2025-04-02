clear all
close all
%% PART 2
%% Nal 1

u = 4*rand(1, 256) - 2;
y = procesf(u');

figure;
subplot(2,1,1);
plot(u);
xlabel('t');
ylabel('input');

subplot(2,1,2);
plot(y);
xlabel('t');
ylabel('output');

%%
n_s = 1; % choose whatever order you like
d = 1;
Ts = 1;

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

G = tf(b, [1 a], Ts);

%%
load MEASUREF.DAT
UU = MEASUREF(:,2);
YY = MEASUREF(:,3);

Ts = 1;
t = 0:Ts:256;
N = length(t);
T = t(2) - t(1);
f = (0:N-1)/(N*T);

% iz druge vaje:
[mag(:), phase(:)] = bode(G, 2*pi*f);

mag_db = 20*log10(abs(mag));
phase_dgr = phase;

% iz prve vaje:
T = 1;
Unal1 = fft(UU)/T;
Ynal1 = fft(YY)/T;
freq_vector = (0:N-1)/(N*T); % od 0 do 1
freq_vector_halfed = freq_vector(1:N/2);

Gnal1 = Ynal1./Unal1;
Gdb_nal1 = 20*log10(abs(Gnal1));
faza = unwrap(angle(Gnal1))*180/pi; %Funkcija unwrap se uporablja za odstranjevanje skokov v faznem kotu

% Narišemo amplitudno karakteristiko v linearno-logaritemski skali
figure;
subplot(2,1,1);
semilogx(freq_vector_halfed, Gdb_nal1(1:N/2));
xlabel('Frekvenca [Hz]');
ylabel('Amplituda [dB]');
title('Amplitudna karakteristika');
hold on; 
semilogx(freq_vector_halfed, mag_db(1:N/2));
xlabel('Frekvenca [Hz]');
ylabel('Amplituda [dB]');
title('Amplitudna karakteristika bode');
legend('iz 1. vaje', 'iz 2. vaje');

% Bode output
subplot(2,1,2);
semilogx(freq_vector_halfed, faza(1:N/2));
xlabel('Frekvenca [Hz]');
ylabel('Fazni kot [°]');
title('Fazna karakteristika');
hold on;
semilogx(freq_vector_halfed, phase_dgr(1:N/2));
xlabel('Frekvenca [Hz]');
ylabel('Fazni kot [°]');
title('Fazna karakteristika bode');
legend('iz 1. vaje', 'iz 2. vaje');
