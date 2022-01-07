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

nMeasures = 24;              % Number of microphones
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

%% working on 1 file
% Load the signal
    fileName=strcat(dir,'24.wav'); 
    [y, fs]=audioread(fileName);
    % Cut the recordings accoding to input signal length
    y=y(1:length(x));
    % Compute the impulse response using the function extractirnoise
    [ir] = extractirnoise(x, y, nfft);
    plot(t, abs(ir))
    % Find the first (and highest) impulse of the impulse response
    ir=abs(ir);
    [peaks,loc] = maxk(ir, 2);
    % Double check if we are finding back the correct direct path length
    toa=t(loc(1));
    directPathTimeOfArrival = toa;
    directPathLength = directPathTimeOfArrival*speed_of_sound;
    %[second_peak, loc2]=max(ir(ir<peak));
    toa2=t(loc(2));
    delays=toa2;
    firstReflPath=delays*speed_of_sound;

%% cycle on all files

directPathTimeOfArrival=zeros(1,nMeasures);
directPathLength=zeros(1,nMeasures);
delays=zeros(1,nMeasures);
firstReflPath=zeros(1,nMeasures);

for i = 1:nMeasures            % For each microphone signal
    % Load the signal
    fileName=strcat(dir, num2str(i), '.wav'); 
    [y, fs]=audioread(fileName);
    % Cut the recordings accoding to input signal length
    y=y(1:length(x));
    % Compute the impulse response using the function extractirnoise
    [ir] = extractirnoise(x, y, nfft);
    ir=abs(ir);
    %plot(t, abs(ir))
    % Find the first (and highest) impulse of the impulse response
    [peak,loc] = maxk(ir, 2);
    % Double check if we are finding back the correct direct path length
    toa=t(loc(1));
    directPathTimeOfArrival(i) = toa;
    directPathLength(i) = directPathTimeOfArrival(i)*speed_of_sound;
    %finding first reflection path
    toa2=t(loc(2));
    delays(i)=toa2;
    firstReflPath(i)=delays(i)*speed_of_sound;
    
    % Plot the estimate impulse response
    %nexttile
    %plot(t, ir);
    %xlim([0 0.05]);
    %xlabel('Time (sec)');
    %title(['Mic: ',num2str(n)]);
    
end

%% SPEAKER TO MIC DISTANCE (DIRECT PATH LENGTH)
% Print on screen the estimated distance from the source
fprintf(sprintf('Direct path length %f m\n', directPathLength));

%% MIC TO REFLECTORS DISTANCE COMPUTATION
% Put here the difference between first reflection and direct sound time of
% arrivals (TOA)

% Inspecting the impulse responses determine the delay of the first
% reflection
%delay = ; %[s]
% Compute the distance from the reflector
distance = firstReflPath-directPathLength;

% Print on screen the estimated distance from the reflector
%fprintf(sprintf('Average distance between first path and reflector %f m\n', ...
 %   distance));