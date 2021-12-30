%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework 4 exercise 3a
% Evaluation of reflections in the signals
% We compute the impulse response using the white noise source signal
% from which the first direct path is computed. In addition by a visual 
% inspection we estimate the first reflection time of arrival. Both the 
% distance from the source and the distance from the reflectors are 
% computed using the estimated times.
%
% Musical Acoustic Course
% Mirco Pezzoli
% 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc 

%% Setup
addpath('Functions')

nMic = ;              % Number of microphones
speed_of_sound = 343.8; % [m]/[s]

typeOfSignal =  % Noise
dir = % File directory

inputSignalDir = ;            % Source signal directory
inputSignalFileName = ;    % Source signal name

%% Impulse response estimation

% The input signal must be known, load the input signal
[x, fs] = ;
nfft = fs;              % Number of fft points
t = (0:1/fs:nfft/fs);   % Time axis
t = t(1:end-1);         
figure;

for n = 1:nMic            % For each microphone signal
    % Load the signal
    
    % Cut the recordings accoding to input signal length
    
    % Compute the impulse response using the function extractirnoise
    [ir] =

    % Find the first (and highest) impulse of the impulse response

    % Double check if we are finding back the correct direct path length
    directPathTimeOfArrival = 
    directPathLength = 
    
    % Plot the estimate impulse response
    nexttile
    plot(t, ir);
    xlim([0 0.05]);
    xlabel('Time (sec)');
    title(['Mic: ',num2str(n)]);
    
end

%% SPEAKER TO MIC DISTANCE (DIRECT PATH LENGTH)
% Print on screen the estimated distance from the source
fprintf(sprintf('Direct path length %f m\n', directPathLength));

%% MIC TO REFLECTORS DISTANCE COMPUTATION
% Put here the difference between first reflection and direct sound time of
% arrivals (TOA)

% Inspecting the impulse responses determine the delay of the first
% reflection
delay = ; %[s]
% Compute the distance from the reflector
distance = 

% Print on screen the estimated distance from the reflector
fprintf(sprintf('Average distance between first path and reflector %f m\n', ...
    distance));
