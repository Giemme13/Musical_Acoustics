%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework 4 exercise 4a
% Radiance pattern estimation directly from the signals
% We estimate the radiance pattern using the white noise source signal.
% The microphone signals are windowed accordingly to the first reflection
% time estimated in the previous exercises. 
%
% Musical Acoustic Course
% Mirco Pezzoli
% 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
close all
clc

%% Setup
fs = 48000;         % Sampling frequency
nfft = fs;          % Number of fft points
nMic = ;          % Number of microphones

R =             % Distance between source and microphones

typeOfSignal = % Noise
dir = % Signals directory

%% RADIANCE COMPUTATION THROUGH THE USE OF THE DIRECT RECORDED SIGNAL
% Window the signals according to the reflection time and estimate the
% radiance pattern. 

sig = [];   % Signal structure

% Early reflections attenuation. We consider only 2*ns+1 samples of the signal
ns = 350; % Ideal ns should be on the direct path
t = 
t = t(1:end-1);

TOA_directSignal =% TOA
TOA_firstReflection =  % First reflection TOA

w =       % Window
figure;
tiledlayout('flow');

for n = 1:nMic            % For each microphone signal
    % Load the signal
    
    
    
    % The window is applied to the signal from the 0 time instant.
    % Check the effect of differen window sizes.
    y_w = 
    
    % Plot the estimate impulse response
    nexttile
        % Plot the mic signal
    hold on
        % Plot the window over the signal
    % Plot the TOA using stem (see doc stem)
    
%     % Plot the first reflection TOA using stem (see doc stem)
    
    % Plot the windowed signal with thicker line
    
    hold off
    % Add a legend
    legend('signal', 'window', 'direct TOA', 'First reflection TOA', ...
        'windowed signal');
    xlim([0 0.05]);             % Limit the plot btw 0 and 0.05s
    xlabel('Time (sec)');
    title(['Mic: ', num2str(n)]);
    sig = % Add current signal to the structure sig
end

%% Radiance estimation

SIG = % FFT of the windowed signal

rad_patt = % Compute the radiance pattern magnitude

%% Radiance pattern plot

% Frequencies must be centered at a frequency bins
ctr_freqs = [150 250 500 750 1000 1500 2000 5000 7500 10000];

% Microphone angle wrt the speaker frontal face
angs = 

% Plot the estimated radiance pattern using the provided function
% radianceplot

radianceplot(ctr_freqs, rad_patt, angs, 'noise direct: ');

