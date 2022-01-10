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

nMic = 24;                         % Number of microphone signals
typeOfSignal = '/sweep/';          % Sweep
dir = ['./recordings', typeOfSignal]; % File directory

fs = 48000;                        % Sampling frequency
speed_of_sound = 343.8;            %[m]/[s]
duration = 10;                     %[s] duration of sweep signal
f1 = 50;
f2 = 22000;

%% Impulse response estimation

% This is the input signal. The following function computes 1/FFT(sweep). 
% We will thus need to MULTIPLY this signal for FFT(output) in order to 
% apply deconvolution (i.e. computing the radiation pattern)
% Use the provided synthSweep function.
[sweep, invsweepfft, sweepRate] = synthSweep(duration,fs,f1,f2);

directPathTimeOfArrival = zeros(1, nMic);
directPathLength = zeros(1, nMic);
envHigh = cell(1,nMic);

figure(1);
tiledlayout('flow', 'padding', 'tight');

for i = 1:nMic    % For each microphone signal
    % Load the signal
    fileName = strcat(dir, num2str(i), '.wav');
    [y] = audioread(fileName);
    
    % Comput the impulse response using the function extractirsweep
    [ir] = extractirsweep(y, invsweepfft);
    [envHigh{1,i}, envLow] = envelope(ir, 30, 'peak');

    % Setting up time scale for computed ir
    nfft = fs;                      % Number of fft points
    t = (0:1/fs:length(ir)/fs);     % Time axis
    t = t(1:end-1);
    
    % Find the first impulse of the impulse response
    [peak,loc] = max(ir);
    
    % Double check if we are finding back the correct direct path length
    directPathTimeOfArrival(i) = t(loc);
    directPathLength(i) = directPathTimeOfArrival(i)*speed_of_sound;
    
    % Plot the estimated impulse response
    nexttile
    hold on
    plot(t, ir);
    plot(t, envHigh{1,i})
    hold off
    grid on
    xlim([0 0.05]);
    xlabel('Time (sec)');
    title(['Mic: ', num2str(i)]);
    if i == 1
        legend('impulse response', 'envelope')
    end
end

%% SPEAKER TO MIC DISTANCE (DIRECT PATH LENGTH)
% Print on screen the estimated distance from the source
figure(2)
plot(1:nMic, directPathLength)
xlabel('Measurement'), ylabel('Distance highest peak')

fprintf(sprintf('Direct path length %f m\n', mean(directPathLength)));

%% MIC TO REFLECTOR DISTANCE COMPUTATION
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
delay = mean(delays);    %[s]
fprintf(sprintf('Average delay between direct sound and first reflection: %f s\n', ...
    delay));

% Compute the distance from the reflectors
distance = delay*speed_of_sound;

% Print on screen the estimated distance from the walls
fprintf(sprintf('Average distance between reflector and microphone %f m\n', ...
    distance));
