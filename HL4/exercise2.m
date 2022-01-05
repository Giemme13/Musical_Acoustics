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

nMic = 1;% Number of microphones
typeOfSignal = 'sweep/';% Noise or sweep
dir = ['./recordings/' typeOfSignal];      % Recordings directory
fileName=strcat(dir, '1.wav');
speed_of_sound = 343.8; % [m]/[s]

%% Plot the signal autocorrelatio

% Load the signal
[x, Fs]=audioread(fileName);
% Time length
time_length=length(x)/Fs;
% Time axis
t=1:length(x);
% Auto correlation
xc=xcorr(x);
%Half of the autocorrelation
xc=xc((length(xc)/2):end);
    
%% Plot the autocorrelation
xlimit=length(t)/32;
figure(1)
plot(t, xc);
title(['Mic: ', num2str(n)]);
axis([0 xlimit -30 50]);    % Limit the axis
xlabel('Time (sec)');
%% MIC TO REFLECTORS DISTANCE COMPUTATION
% Put here the difference between first reflection 
%delay = %[s]

%distance = ;%[m]

%fprintf(sprintf('Average distance between mic and first reflection %f m\n', distance));
