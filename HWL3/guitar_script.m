%% Setup

clear;
close all;
clc;

addpath('Functions')
simulink_folder = './';       % Simulink projects folder
addpath(simulink_folder);

fs = 44100;                         % Sampling frequency
signalLen = 5;                      % Signal length
t = [0:1/fs:signalLen-1/fs];        % Time axis

%% Plate simulation
% run simulink simulation
sim('guitar_plate'); 

% Resampling
I1 = resample(I, t);

% Plot in time
figure(1)
plot(I1.time, I1.data);

% Normalize the signal
soundWave1 = I1.data;
soundWave1 = soundWave1./max(abs(soundWave1));

%% Plot and play
% Plot magnitude and phase
[S1, magS1, angleS1, f1, df1] = myFFT(soundWave1, fs);

h1 = figure(2);
plotFFT_linearFreqScale(magS1, angleS1, f1, df1, fs, 2000, h1);

figure(3);
plot(t, soundWave1);
title('Signal in time'), xlabel('Time [s]')

sound(soundWave1, fs);                           % Play the sound

fileName1 = 'guitar_plate.wav';     % Audio file path
disp('Save file on disk...')                    % Save on disk
%audiowrite(fileName1, soundWave1, fs);

%% Plate+string simulation

% run simulink simulation
sim('guitar_plate_string'); 

% Resampling
I2 = resample(I, t);

% Plot in time
figure(4)
plot(I2.time, I2.data);

% Normalize the signal
soundWave2 = I2.data;
soundWave2 = soundWave2./max(abs(soundWave2));

%% Plot and play
% Plot magnitude and phase
[S2, magS2, angleS2, f2, df2] = myFFT(soundWave2, fs);

h2 = figure(5);
plotFFT_linearFreqScale(magS2, angleS2, f2, df2, fs, 2000, h2);

figure(6);
plot(t, soundWave2);
title('Signal in time'), xlabel('Time [s]')

sound(soundWave2, fs);                           % Play the sound

fileName2 = 'guitar_plate+string.wav';     % Audio file path
disp('Save file on disk...')                    % Save on disk
audiowrite(fileName2, soundWave2, fs);
