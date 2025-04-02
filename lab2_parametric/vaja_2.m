close all
clear all
%% Data:

load MEASUREF.DAT
t = MEASUREF(:,1);
u = MEASUREF(:,2);
y = MEASUREF(:,3);
n = MEASUREF(:,4);

N = length(t);
T = t(2) - t(1);
Ts = 1;
F = 1/(N*T);
f = (0:N-1)/(N*T);
%% Graph:

figure;
plot(t,y);
xlabel("Time");
ylabel("y");
hold on;
stairs(t,u);
xlabel("t");
ylabel("u");
legend("Output", "Input");
%% Naloga 1

% n_s of system (red je 2)  
n_s = 2;

% delay d from graph
d = 0;
%% Naloga 2
% %%%%%%%%%%%%%% Without for loop - try just one order %%%%%%%%%%%%%%%%%%%%%%
% 
% % TRANSFER FUNCTION:
% % Gz = Y(z)/U(z) = ((b1*z^-1 + b2*z^-2)/(1 + a1*z^-1 + a2*z^-2))*z^-d
% 
% % DIFFERENCE-EQUATION FORM:
% % y(k) = - a1*y(k-1) - a2*y(k-2) + b1*u(k-1) + b2*u(k-2)
% 
% % y = PSI * theta
% PSI = [-y(2:end-1) -y(1:end-2) u(2:end-1) u(1:end-2)];
% Y = y(3:end);
% theta = PSI \ Y;
% 
% a1 = theta(1);
% a2 = theta(2);
% b1 = theta(3);
% b2 = theta(4);
% 
% G = tf([b1 b2], [1 a1 a2], Ts);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%% With for loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
n_s = 6; % choose whatever order you like
d = 0;
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
%% Naloga 3

% iz druge vaje:
[mag(:), phase(:)] = bode(G, 2*pi*f);

mag_db = 20*log10(abs(mag));
phase_dgr = phase;

% iz prve vaje:
Unal1 = fft(u)/T;
Ynal1 = fft(y)/T;
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

%% Naloga 4

y_sim = lsim(G, u, t); % simulated output of discrete model
y = MEASUREF(:,3); % measured output

figure;
plot(t, y, 'r', 'LineWidth', 1.5); % Plot measured output in red
hold on;
plot(t, y_sim, 'b--', 'LineWidth', 1.5); % Plot simulated output in blue dashed line
hold off;
xlabel('Time (s)');
ylabel('Output');
legend('Measured Output', 'Simulated Output');
title('Comparison of Measured and Simulated Outputs');
%% covariance matrix and standard deviations of estimated parameters:
% finish this task

%% Naloga 5

sys_c = d2c(G); % Continuous-time transfer function
