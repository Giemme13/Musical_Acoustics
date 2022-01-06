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
addpath('Functions');

fs = 48000;          % Sampling frequency
nfft = fs;           % Number of fft points
nMic = 1;            % Number of microphones

R = 2.57;            % Distance between source and microphones

typeOfSignal = 'noise/';                % Noise
dir = ['./recordings/' typeOfSignal];   % Signals directory

%% RADIANCE COMPUTATION THROUGH THE USE OF THE DIRECT RECORDED SIGNAL
% Window the signals according to the reflection time and estimate the
% radiance pattern. 

sig = [];   % Signal structure

% Early reflections attenuation. We consider only 2*ns+1 samples of the signal
ns = 50; % Ideal ns should be on the direct path
t = (0:1/fs:((2*ns)+1)/fs);
t = t(1:end-1);

TOA_directSignal = 0;                % TOA
TOA_firstReflection = 0.0024;        % First reflection TOA

w = hann(length(t));             % Window
figure;
tiledlayout('flow');

for n = 1:nMic            % For each microphone signal
    % Load the signal
    fileName = strcat(dir, num2str(1), '.wav');  % i-th file name
    [y, Fs] = audioread(fileName);               % read i-th audio
    y = y./max(y);
    
    % The window is applied to the signal from the 0 time instant.
    % Check the effect of differen window sizes.
    y = y(1:length(t));
    y_w = y.*w;
    
    % Plot the estimate impulse response    
    nexttile
    % Plot the mic signal
    plot(t, y(1:length(t)))
    hold on
    % Plot the window over the signal
    plot(t, w)
    % Plot the TOA using stem (see doc stem)
    index = find(abs(t-TOA_directSignal)==min(abs(t-TOA_directSignal)));
    stem(t(index), y(index))
    % Plot the first reflection TOA using stem (see doc stem)
    index = find(abs(t-TOA_firstReflection)==min(abs(t-TOA_firstReflection)));
    stem(t(index), y(index))
    % Plot the windowed signal with thicker line
    plot(t, y_w, 'linewidth', 2)
    hold off
    % Add a legend
    legend('signal', 'window', 'direct TOA', 'First reflection TOA', ...
        'windowed signal');
    xlim([0 0.05]);             % Limit the plot btw 0 and 0.05s
    xlabel('Time (sec)');
    title(['Mic: ', num2str(n)]);
    sig = y_w;      % Add current signal to the structure sig
end

%% Radiance estimation

SIG = fft(sig, nfft);   % FFT of the windowed signal

rad_patt = abs(SIG);    % Compute the radiance pattern magnitude

%% Radiance pattern plot

% Frequencies must be centered at a frequency bins
ctr_freqs = [150 250 500 750 1000 1500 2000 5000 7500 10000];

% Microphone angle wrt the speaker frontal face
angs = 0;

% Plot the estimated radiance pattern using the provided function
% radianceplot

radianceplot(ctr_freqs, rad_patt, angs, 'noise direct: ');



%% i miei dubbi sono:
% - signal structure: non so cosa farci
% - t lo scrivo già sapendo che prendo ns samples? o devo prendere l'esatta
%   durata del segnale i-esimo? in questo caso dovrei fare un ciclo for per
%   raccogliere tutti i file e per ciascuno fare tutto ciò che c'è in
%   questo eserciizio?
% - i TOA li metto io sulla base dei risultati degli esercizi precedenti?
% - il file caricato lo lascio della sua lunghezza o lo accorcio a t?
% - come si fa il radiation pattern magnitude? è semplicemente abs()?