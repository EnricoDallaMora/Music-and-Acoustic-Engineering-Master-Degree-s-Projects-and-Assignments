clear; close all; clc
profile on;
%% Simscape File for Ground Truth
load('ssc_output.mat')

%% Sampling Frequency
fs = 96e3;

%% Sampling Period
Ts = 1/fs;

%% Simulation Duration
stop_time = 1;  % [seconds]

%% Input Signal
% Fundamental Frequency
f0 = 440;

% Time Axis
t = 0:Ts:stop_time;

% Signal Amplitude
A = 1.5;
vin = A * sin(2*pi*f0*t);

% Music Signal Test (Uncomment the following lines for test)
% [vin, fs_rec] = audioread('guitar_input.wav');
% G = 5;                                  % Gain factor 
% vin = G * (vin(:, 1) + vin(:, 2))/2;    % Convert the signal to MONO
% vin = resample(vin, fs, fs_rec);        % Resampling the signal from 44.1 kHz to 96 kHz
% t = 0:Ts:Ts*(length(vin)-1);

%% Circuit Parameters
% Resistive Elements
Rin = 1;
R1 = 1e4;
Rout = 1e4;

% Dynamic Elements
C1 = 1e-6;
C2 = 1e-9;

%% Setting of Free Parameters (Adaptation Conditions)
Z_8 = R1;
Z_4 = Rout;
Z_6 = Rin;
Z_7 = Ts/(2*C1);
Z_1 = Ts/(2*C2);
Z_5 = Z_6 + Z_7 + Z_8;
Z_2 = Z_5;
Z_3 = 1/(Z_1^(-1) + Z_2^(-1) + Z_4^(-1));

%% Computing Scattering Matrices
B_s = [1,1,1,1];
Z_s = diag([Z_5, Z_6, Z_7, Z_8]);
S_s = eye(4,4) - 2*Z_s*B_s'*inv(B_s*Z_s*B_s')*B_s;

Z_p = diag([Z_1, Z_2, Z_3, Z_4]);
Q_p = [1,1,1,-1];
S_p = 2 * Q_p' * inv((Q_p*inv(Z_p)*Q_p'))*Q_p*inv(Z_p) - eye(4,4);

%% Initialization of Waves
a = zeros(8, length(t));
b = zeros(8, length(t));

a(6,:) = vin;

%% Initialization of Output Signals
vout = zeros(1, length(t));

%% Simulation Algorithm

for n = 1 : length(t)

    % Forward Scan
    a(4,n) = 0;
    a(8,n) = 0;
    
    % For the first iteration, the capacitors' wave variable is set to zero
    if (n<2)
        a(7,n) = 0;
        a(1,n) = 0;
    else
        a(7,n) = b(7,n-1);
        a(1,n) = b(1,n-1);
    end

    b(5,n) = S_s(1,:) * [a(5, n); a(6,n); a(7,n); a(8,n)];
    a(2,n) = b(5,n);

    b(3,n) = S_p(3,:) * [a(1,n); a(2,n); a(3, n); a(4,n)];

    % Local Root Scattering
    % Hint: Use the function 'antiparallel_diodes' to compute the Local Root
    % Scattering
    a(3,n) = antiparallel_diodes(b(3,n), Z_3);

    % Backward Scan
    b(1,n) = S_p(1,:) * [a(1,n); a(2,n); a(3,n); a(4,n)];
    b(2,n) = S_p(2,:) * [a(1,n); a(2,n); a(3,n); a(4,n)];
    b(4,n) = S_p(4,:) * [a(1,n); a(2,n); a(3,n); a(4,n)];
    a(5,n) = b(2,n);

    b(6,n) = S_s(2,:) * [a(5,n); a(6,n); a(7,n); a(8,n)];
    b(7,n) = S_s(3,:) * [a(5,n); a(6,n); a(7,n); a(8,n)];
    b(8,n) = S_s(4,:) * [a(5,n); a(6,n); a(7,n); a(8,n)];


    % Read Output
    vout(1,n) = (a(4,n) + b(4,n)) / 2;

end

% Uncomment the following line to hear the Diode Clipper
% sound(vout, fs)

%% Output Plots

plot_lim = 5/f0; % Limit the plot to just 5 periods of the output signal

figure
set(gcf, 'Color', 'w');
plot(gt(1, :), gt(2, :), 'r', 'Linewidth', 2);
hold on;
plot(t, vout, 'b--', 'Linewidth', 2);
hold on;
plot(t, vin, 'Linewidth', 2);
grid on;
xlim([0, plot_lim]);
xlabel('Time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{out}}$ [V]','Fontsize',16,'interpreter','latex');
legend('Simscape','WDF','Fontsize',16,'interpreter','latex');
title('Output Signal','Fontsize',18,'interpreter','latex');

%% Error Plots

figure
set(gcf, 'Color', 'w');
hold on;
plot(t, vout - gt(2, :), 'k', 'Linewidth', 2);
grid on;
xlim([0, plot_lim]);
xlabel('Time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{out}}$ [V]','Fontsize',16,'interpreter','latex');
title('Error Signal','Fontsize',18,'interpreter','latex');

%% Compute Mean Squared Error (MSE)

mse = mean((vout - gt(2, :)).^2);
disp('MSE = ')
disp(mse)
profile viewer;
