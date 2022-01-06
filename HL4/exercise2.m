%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework 4 exercise 2
% Reflection inspection
% By a visual inspection of the recorded signals we evaluate the first
% reflection time from which the distance of the micrphones from the
% reflectors can be inferred.
%
% Musical Acoustic Course
% Mirco Pezzoli
% 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc 

%% First reflection estimation using autocorrelation
% Collecting autocorrelations and plotting them in order to analyze
% reflections. Try different takes and input signals

addpath('Functions');

nMic = 1;                                  % Number of microphones
typeOfSignal = 'sweep/';                   % Noise or sweep
dir = ['./recordings/' typeOfSignal];      % Recordings directory
fileName = strcat(dir, '1.wav');
speed_of_sound = 343.8;                    % [m]/[s]

%% Plot the signal autocorrelation

% Load the signal
[x, Fs] = audioread(fileName);
% Time length
time_length = length(x)/Fs;
% Time axis
t = linspace(0,time_length,length(x))';
% Auto correlation
[xc, lags] = xcorr(x, 'normalized');
% Half of the autocorrelation
xc = xc(round((length(xc)/2)):end);
lags = lags(round((length(lags)/2)):end);

% Plot 
figure(1)
subplot(2,1,1)
plot(t, x)
xlabel('Time [s]', 'fontsize', 20)
title('Original signal', 'fontsize', 20)
subplot(2,1,2)
plot(t, xc);
xlabel('Time [s]', 'fontsize', 20)
title('Autocorrelation', 'fontsize', 20)
sgtitle(['Mic: ', num2str(nMic)], 'fontsize', 20)


%% MIC TO REFLECTORS DISTANCE COMPUTATION
% Put here the difference between first reflection
index_reflex = find(abs(t-0.1)==min(abs(t-0.1)));
xc_reflex = xc(1:index_reflex);
t_reflex = t(1:index_reflex);
figure()
plot(t_reflex,xc_reflex);
[max,locs] = findpeaks(abs(xc_reflex));
locs_times = t_reflex(locs);
figure()
stem(locs_times,max)
distances = locs_times*speed_of_sound;
%delay = 0.06;                        %[s]
%distance = delay*speed_of_sound;       %[m]

% fprintf(sprintf('Average distance between mic and first reflection %f m\n', distance));
