%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exercise 6
% Complete model of a guitar
% This script call the simulink file that implements the simple electrical
% equivalent of a guitar. The simulation is executed and the resulting
% current is resampled with the desired sampling frequency. Finally the
% sound is plotted in time and frequency domain and saved on disk. 
%
% Musical Acoustic Course
% Mirco Pezzoli
% 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;

%% Setup
fs = 44100;                         % Sampling frequency
signalLen = 3;                      % Signal length
t = 0:1/fs:signalLen-1/fs;        % Time axis

fileName = 'simple_guitar.wav';     % Audio file paths

%% Simulation
% run the simulink simulation using the command sim (see doc sim).
sim('SimpGuitarModel.slx'); 
% The variable I contains non constant time intervals between samples.
% Resample the data using resample function in order to obtain an equally
% sampled signal I1
I1 = resample(ans.I,t);
%% Plot the resampled signal in time

figure(1)
plot(I1.time,I1.data);
grid on;
title("Model without the string")
xlabel("$Time \ [s]$", Interpreter='Latex');
ylabel("$Current \ [A]$", Interpreter = 'Latex');

% Normalize the signal
out = (I1.Data./max(abs(I1.Data)));


%% Plot and play
% Plot the signal frequency content as magnitude and phase

OUT = fft(I1.data,fs);
freq = 0:fs-1;

figure(2)
sgtitle("Frequency spectrum");
subplot(2,1,1);
plot(freq,db(abs(OUT)),LineWidth=2);
grid on;
xlabel("$Frequency \ [Hz]$", Interpreter='Latex');
xlim([0 1500]);
ylabel("$Magnitude \ [dB]$", Interpreter='Latex');


subplot(2,1,2);
plot(freq,angle(OUT)*180/pi,LineWidth=2);
grid on;
xlabel("$Frequency \ [Hz]$", Interpreter='Latex');
xlim([0 1500]);
ylabel("$Phase \ [rad]$", Interpreter='Latex');



%% Play the sound

sound(out,fs);

%% Save on disk

disp('Save file on disk...')                    
audiowrite(fileName, out, fs);

%% Spectrogram
windowwidth = 0.04;
dt = 1/fs;
Nw = round(windowwidth/dt);
Nov = round(Nw/5);

figure(3)
spectrogram(out, Nw, Nov, [], fs);
xlim([0,8]);
ylim([0 3]);
title("Model without the string");
