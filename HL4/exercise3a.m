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

nMic = 24;              % Number of microphones
speed_of_sound = 343.8; % [m]/[s]

typeOfSignal = 'noise/'; % Noise
dir = ['./recordings/' typeOfSignal]; % File directory

inputSignalDir = './input signals/';            % Source signal directory
inputSignalFileName = strcat(inputSignalDir, 'noise.wav');    % Source signal name

%% Impulse response estimation

% The input signal must be known, load the input signal
[x, fs] = audioread(inputSignalFileName);
nfft = fs;              % Number of fft points
t = (0:1/fs:nfft/fs);   % Time axis
t = t(1:end-1);         


%% cycle on all files

directPathTimeOfArrival = zeros(1, nMic);
directPathLength = zeros(1, nMic);
firstReflTimeOfArrival = zeros(1, nMic);
firstReflPathLength = zeros(1, nMic);

figure(1)
tiledlayout('flow', 'padding', 'tight');

for i = 1:nMic            % For each microphone signal
    % Load the signal
    fileName = strcat(dir, num2str(i), '.wav'); 
    [y, fs] = audioread(fileName);
    % Cut the recordings accoding to input signal length
    y = y(1:length(x));
    % Compute the impulse response using the function extractirnoise
    [ir] = extractirnoise(x, y, nfft);
    ir = abs(ir);
    %plot(t, abs(ir))
    % Find the first (and highest) impulse of the impulse response
    [peak,loc] = max(ir);
    % Double check if we are finding back the correct direct path length
    TOA = t(loc);
    directPathTimeOfArrival(i) = TOA;
    directPathLength(i) = directPathTimeOfArrival(i)*speed_of_sound;
    %Plot the estimate impulse response
    nexttile
    plot(t, ir);
    grid on
    xlim([0.007 0.013]);
    xlabel('Time (sec)');
    title(['Mic: ',num2str(i)]);
end

%% SPEAKER TO MIC DISTANCE (DIRECT PATH LENGTH)
% Print on screen the estimated distance from the source
directPathLength = mean(directPathLength);
fprintf(sprintf('Direct path length %f m\n', directPathLength));

%% MIC TO REFLECTORS DISTANCE COMPUTATION
% Put here the difference between first reflection and direct sound time of
% arrivals (TOA)
firstReflTimeOfArrival = [0.011875, 0.0111875, 0.0112083, 0.0112083, 0.0114167, 0.0114167, ...
    0.0114583, 0.0115417, 0.0115833, 0.0116458, 0.0117708, 0.0118125, ...
    0.0117917, 0.0118125, 0.0118333, 0.0116458, 0.0113958, 0.0113542, ...
    0.0113125, 0.0112708, 0.0111667, 0.0112083, 0.0111875, 0.0111875];
delays = firstReflTimeOfArrival - directPathTimeOfArrival;
TOA_direct = mean(directPathTimeOfArrival);
TOA_reflex = mean(firstReflTimeOfArrival);
% Inspecting the impulse responses determine the delay of the first
% reflection
delay = mean(delays); %[s]
fprintf(sprintf('Average delay between direct sound and first reflection: %f m\n', ...
    delay));
% Compute the distance from the reflector
distance = delay*speed_of_sound;

% Print on screen the estimated distance from the reflector
fprintf(sprintf('Average distance between reflector and microphone: %f m\n', ...
    distance));
