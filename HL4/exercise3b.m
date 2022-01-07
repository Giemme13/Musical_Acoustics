%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework 4 exercise 3b
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

nMeasures = 24; 
typeOfSignal =  'sweep/';                   % Sweep
dir = ['./recordings/' typeOfSignal];  % File directory

fs = 48000;                                 % Sampling frequency
speed_of_sound = 343.8;                     %[m]/[s]            
duration = 10;                               %[s] duration of sweep signal
f1=50;
f2=22000;
% This is the input signal. The following function computes 1/FFT(sweep). 
% We will thus need to MULTIPLY this signal for FFT(output) in order to 
% apply deconvolution (i.e. computing the radiation pattern)
% Use the provided synthSweep function.
[sweep invsweepfft sweepRate] = synthSweep(duration,fs,f1,f2);


% Setting up time scale for computed ir
nfft = fs;              % Number of fft points
t = (0:1/fs:nfft/fs);   % Time axis
t = t(1:end-1);  

%% working on single file

% Load the signal
fileName=strcat(dir, '1.wav');
[x]=audioread(fineName);
    
% Comput the impulse response using the function extractirsweep
[ir] = 
    
% Find the first impulse of the impulse response
[peak,loc] = max(abs(ir));    

% Double check if we are finding back the correct direct path length
toa=t(loc);
directPathTimeOfArrival = toa;
directPathLength = directPathTimeOfArrival*speed_of_sound;


%% cycle on all files
for n = 1:nMeasures    % For each microphone signal
    % Load the signal

    
    % Comput the impulse response using the function extractirsweep
    [ir] =

    % Setting up time scale for computed ir
    t = 
    
    % Find the first impulse of the impulse response
    

    % Double check if we are finding back the correct direct path length
    directPathTimeOfArrival = 
    directPathLength(n) = 
    
    % Plot the estimated impulse response
    nexttile
    plot(t, ir);
    xlim([0 0.05]);
    xlabel('Time (sec)');
    title(['Mic: ', num2str(n)]);
    
end

%% SPEAKER TO MIC DISTANCE (DIRECT PATH LENGTH)
% Print on screen the estimated distance from the source
figure;
plot(1:nMeasures, directPathLength)
xlabel('Measurement'), ylabel('Distance highest peak')

fprintf(sprintf('Direct path length %f m\n', directPathLength));

%% MIC TO REFLECTOR DISTANCE COMPUTATION
% Put here the difference between first reflection and direct sound time of
% arrivals (TOA)

% Inspecting the impulse responses determine the delay of the first
% reflection
delay = ; %[s]
% Compute the distance from the reflectors
distance = ;

% Print on screen the estimated distance from the walls
fprintf(sprintf('Average distance between first path and reflector %f m\n', ...
    distance));
