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

nMic = 24;                      % Number of microphones
speed_of_sound = 343.8;         % [m]/[s]

typeOfSignal = 'noise/';        % Noise
dir = ['./recordings/' typeOfSignal]; % File directory

inputSignalDir = './input signals/';            % Source signal directory
inputSignalFileName = strcat(inputSignalDir, 'noise.wav');    % Source signal name

%% Impulse response estimation

% The input signal must be known, load the input signal
[x, fs] = audioread(inputSignalFileName);
nfft = fs;              % Number of fft points
t = (0:1/fs:nfft/fs);   % Time axis
t = t(1:end-1);         

directPathTimeOfArrival = zeros(1, nMic);
directPathLength = zeros(1, nMic);
envHigh = cell(1,nMic);

figure(1);
tiledlayout('flow', 'padding', 'tight');

for i = 1:nMic            % For each microphone signal
    % Load the signal
    fileName = strcat(dir, num2str(i), '.wav'); 
    [y] = audioread(fileName);
    
    % Cut the recordings accoding to input signal length
    y = y(1:length(x));
    
    % Compute the impulse response using the function extractirnoise
    [ir] = extractirnoise(x, y, nfft);
    [envHigh{1,i}, envLow] = envelope(ir, 30, 'peak');
    
    % Find the first (and highest) impulse of the impulse response
    [peak,loc] = max(ir);
    
    % Double check if we are finding back the correct direct path length
    directPathTimeOfArrival(i) = t(loc);
    directPathLength(i) = directPathTimeOfArrival(i)*speed_of_sound;
    
    %Plot the estimate impulse response
    nexttile
    hold on
    plot(t, ir);
    plot(t, envHigh{1,i})
    hold off
    grid on
    xlim([0 0.05]);
    xlabel('Time (sec)');
    title(['Mic: ',num2str(i)]);
    if i == 1
        legend('impulse response', 'envelope')
    end
end

%% SPEAKER TO MIC DISTANCE (DIRECT PATH LENGTH)
% Print on screen the estimated distance from the source
figure(2)
plot(1:nMic, directPathLength)
xlabel('Measurement', 'fontsize', 20), ylabel('Distance highest peak', 'fontsize', 20)

fprintf(sprintf('Direct path length %f m\n', mean(directPathLength)));

%% MIC TO REFLECTORS DISTANCE COMPUTATION
% Put here the difference between first reflection and direct sound time of
% arrivals (TOA)

firstReflTimeOfArrival = zeros(1, nMic);
for i = 1:nMic
    env_ir = envHigh{1,i};
    [peaks,locs] = findpeaks(env_ir(1:length(env_ir)/2));
    [reflValue, reflIndex] = maxk(peaks, 2);
    firstReflTimeOfArrival(i) = t(locs(reflIndex(2)));
end

delays = firstReflTimeOfArrival - directPathTimeOfArrival;

% Inspecting the impulse responses determine the delay of the first
% reflection
delay = mean(delays);   %[s]
fprintf(sprintf('Average delay between direct sound and first reflection: %f m\n', ...
    delay));

% Compute the distance from the reflector
distance = delay*speed_of_sound;

% Print on screen the estimated distance from the reflector
fprintf(sprintf('Average distance between reflector and microphone: %f m\n', ...
    distance));
