% MATLAB kod za generisanje PPG signala i poređenje metoda

clc;
clear;
close all;

% Parametri simulacije
fs = 100;                 % Frekvencija uzorkovanja (Hz)
t = 0:1/fs:10;            % Vremenski opseg (10 sekundi)
f_heart = 1.2;            % Frekvencija pulsa (Hz)
PPG_clean = 0.5 * sin(2 * pi * f_heart * t) + 0.5; % Idealni PPG signal

% Dodavanje šumova
motion_noise = 0.1 * sin(2 * pi * 0.5 * t);  % Pokret (0.5 Hz sinusni šum)
ambient_noise = 0.05 * randn(size(t));       % Ambijentalni šum (beli šum)
PPG_noisy = PPG_clean + motion_noise + ambient_noise; % Ukupan šumni signal

%% Standardna metoda: Butterworth niskopropusni filtar
fc = 2;  % Granična frekvencija (2 Hz)
[b, a] = butter(4, fc / (fs / 2), 'low');  % Butterworth filtar (4. reda)
PPG_standard = filtfilt(b, a, PPG_noisy); % Filtrovani signal standardnom metodom

%% Predložena metoda: Kalmanov filter
% Inicijalizacija Kalmanovog filtra
n = length(t);
X = [0; 0];            % Početno stanje: [DC; AC]
P = eye(2);            % Početna kovarijansa
F = [1 0; 0 1];        % Matrica tranzicije stanja
H = [1 0; 0 1];        % Matrica merenja
Q = 0.01 * eye(2);     % Procesni šum
R = 0.1 * eye(2);      % Šum merenja
PPG_filtered = zeros(1, n); % Memorija za filtrirani signal

% Primena Kalmanovog filtra
for i = 1:n
    % Predikcija
    X_pred = F * X;
    P_pred = F * P * F' + Q;

    % Ažuriranje
    Z = [PPG_noisy(i); PPG_noisy(i)]; % Merenje [DC, AC]
    K = P_pred * H' / (H * P_pred * H' + R); % Kalmanov dobitak
    X = X_pred + K * (Z - H * X_pred); % Ažuriranje stanja
    P = (eye(2) - K * H) * P_pred; % Ažuriranje kovarijanse

    % Sačuvaj AC komponentu (pulsirajući deo)
    PPG_filtered(i) = X(2);
end

%% Prikaz rezultata na odvojenim grafikonima
figure;

% Grafikon 1: Originalni PPG signal i šumni signal
subplot(3, 1, 1);
plot(t, PPG_clean, 'g', 'LineWidth', 1.5); hold on;
plot(t, PPG_noisy, 'r'); 
title('Original and Noisy PPG Signal');
xlabel('Time (s)'); ylabel('Amplitude');
legend('Clean PPG', 'Noisy PPG');
grid on;

% Grafikon 2: Standardna metoda (Butterworth filtar)
subplot(3, 1, 2);
plot(t, PPG_clean, 'g', 'LineWidth', 1.5); hold on;
plot(t, PPG_standard, 'b', 'LineWidth', 1.5);
title('Standard Method: Butterworth Filter');
xlabel('Time (s)'); ylabel('Amplitude');
legend('Clean PPG', 'Filtered PPG (Butterworth)');
grid on;

% Grafikon 3: Predložena metoda (Kalmanov filter)
subplot(3, 1, 3);
plot(t, PPG_clean, 'g', 'LineWidth', 1.5); hold on;
plot(t, PPG_filtered, 'k', 'LineWidth', 1.5);
title('Proposed Method: Kalman Filter');
xlabel('Time (s)'); ylabel('Amplitude');
legend('Clean PPG', 'Filtered PPG (Kalman)');
grid on;

%% Dodatna analiza: MAPE (Mean Absolute Percentage Error)
MAPE_standard = mean(abs((PPG_clean - PPG_standard) ./ PPG_clean)) * 100;
MAPE_kalman = mean(abs((PPG_clean - PPG_filtered) ./ PPG_clean)) * 100;

disp('Performance Comparison:');
fprintf('MAPE (Standard Method - Butterworth): %.2f%%\n', MAPE_standard);
fprintf('MAPE (Proposed Method - Kalman Filter): %.2f%%\n', MAPE_kalman);
