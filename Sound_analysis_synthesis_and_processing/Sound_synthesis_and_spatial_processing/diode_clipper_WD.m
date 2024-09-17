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
Z_9 = R1;
Z_6 = Rout;
Z_12 = Rin;
Z_11 = Ts/(2*C1);
Z_5 = Ts/(2*C2);
Z_10 = Z_11 + Z_12;
Z_8 = Z_10;
Z_7 = Z_8 + Z_9;
Z_2 = Z_7;
Z_4 =1/(Z_5^(-1) + Z_6^(-1));
Z_1 = Z_4;
Z_3 = 1/(Z_1^(-1) + Z_2^(-1));

%% Computing Scattering Matrices
Z_pr = diag([Z_4, Z_5, Z_6]);
Q_pr = [1,1,1];
S_pr = 2 * Q_pr' * inv(Q_pr*inv(Z_pr)*Q_pr')*Q_pr*inv(Z_pr) - eye(3,3);

B_sl = [1,1,1];
Z_sl = diag([Z_10, Z_11, Z_12]);
S_sl = eye(3,3) - 2*Z_sl*B_sl'*inv(B_sl*Z_sl*B_sl')*B_sl;

B_sr = [1,1,1];
Z_sr = diag([Z_7, Z_8, Z_9]);
S_sr = eye(3,3) - 2*Z_sr*B_sr'*inv(B_sr*Z_sr*B_sr')*B_sr;

Z_pl = diag([Z_1, Z_2, Z_3]);
Q_pl = [1,1,1];
S_pl = 2 * Q_pl' * inv((Q_pl*inv(Z_pl)*Q_pl'))*Q_pl*inv(Z_pl) - eye(3,3);

%% Initialization of Waves
a = zeros(12, length(t));
b = zeros(12, length(t));

a(12,:) = vin;

%% Initialization of Output Signals
vout = zeros(1, length(t));

%% Simulation Algorithm

for n = 1 : length(t)

    % Forward Scan
    a(6,n) = 0;
    a(9,n) = 0;
    if (n<2)
        a(5,n) = 0;
        a(11,n) = 0;
    else
        a(5,n) = b(5,n-1);
        a(11,n) = b(11,n-1);
    end

    b(10,n) = S_sl(1,:) * [a(10, n); a(11,n); a(12,n)];
    a(8,n) = b(10,n);
    b(7,n) = S_sr(1,:) * [a(7, n); a(8,n); a(9,n)];
    a(2,n) = b(7,n);

    b(4,n) = S_pr(1,:) * [a(4, n); a(5,n); a(6,n)];
    a(1,n) = b(4,n);
    b(3,n) = S_pl(3,:) * [a(1,n); a(2,n); a(3, n)];

    % Local Root Scattering
    % Hint: Use the function 'antiparallel_diodes' to compute the Local Root
    % Scattering
    a(3,n) = antiparallel_diodes(b(3,n), Z_3);

    % Backward Scan
    b(1,n) = S_pl(1,:) * [a(1,n); a(2,n); a(3,n)];
    b(2,n) = S_pl(2,:) * [a(1,n); a(2,n); a(3,n)];
    a(4,n) = b(1,n);
    a(7,n) = b(2,n);

    b(5,n) = S_pr(2,:) * [a(4,n); a(5,n); a(6,n)];
    b(6,n) = S_pr(3,:) * [a(4,n); a(5,n); a(6,n)];
    
    b(8,n) = S_sr(2,:) * [a(7,n); a(8,n); a(9,n)];
    b(9,n) = S_sr(3,:) * [a(7,n); a(8,n); a(9,n)];
    a(10,n) = b(8,n);

    b(11,n) = S_sl(2,:) * [a(10,n); a(11,n); a(12,n)];
    b(12,n) = S_sl(3,:) * [a(10,n); a(11,n); a(12,n)];

    % Read Output
    vout(1,n) = (a(6,n) + b(6,n)) / 2;

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
legend('Simscape','WDF','$V_{in}$','Fontsize',16,'interpreter','latex');
title('Output Signal','Fontsize',18,'interpreter','latex');

%% Error Plots

figure
set(gcf, 'Color', 'w');
hold on;
plot(t, vout - gt(2,:), 'k', 'Linewidth', 2);
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
