%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework 4 exercise 4b
% Radiance pattern estimation directly from the signals
% We estimate the radiance pattern using the sine sweep source signal.
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

fs = 48000;          % Sampling frequency
nfft = fs;           % Number of fft points
nMic = 24;           % Number of microphone signals

% Distance between source and microphones
R = [2.7575624, 2.7647250, 2.7647250, 2.7718875, 2.7862124, 2.8005376, ...
     2.8864875, 2.9079750, 2.9366250, 2.9581125, 2.9796000, 3.0082500, ...
     3.0154126, 2.9939251, 2.9724374, 2.9509499, 2.9223001, 2.9008124, ...
     2.8148625, 2.8005376, 2.7862124, 2.7718875, 2.7647250, 2.7647250];

typeOfSignal = 'sweep/';                % Sweep
dir = ['./recordings/' typeOfSignal];   % Signals directory

%% RADIANCE COMPUTATION THROUGH THE USE OF THE DIRECT RECORDED SIGNAL
% Window the signals according to the reflection time and estimate the
% radiance pattern.
sig = [];   % Signal structure

% Early reflections attenuation. We consider only 2*ns+1 samples of the signal
ns = 350; % Ideal ns should be on the direct path
t = (0:1/fs:((2*ns)+1)/fs);
t = t(1:end-1);

% TOA direct sound
TOA_directSignal =  [0.0080208331, 0.0080416668, 0.0080416668, ...
    0.0080624996, 0.0081041669, 0.0081458334, 0.0083958330, 0.0084583331, ...
    0.0085416669, 0.0086041670, 0.0086666662, 0.0087500000, 0.0087708337, ...
    0.0087083336, 0.0086458335, 0.0085833333, 0.0085000005, 0.0084375003, ...
    0.0081874998, 0.0081458334, 0.0081041669, 0.0080624996, 0.0080416668, ...
    0.0080416668];
% TOA first reflection
TOA_firstReflection = [0.011083334, 0.011104167, 0.011083334, ...
    0.011125000, 0.011145833, 0.011333333, 0.011250000, 0.011270833, ...
    0.011708333, 0.011729167, 0.011708333, 0.011729167, 0.011708333, ...
    0.011687500, 0.011687500, 0.011750000, 0.011479166, 0.011500000, ...
    0.011291667, 0.011229167, 0.011166667, 0.011166667, 0.011104167, ...
    0.011083334];

w = ones(length(t),1);      % Window
w(1) = 0;
w(end) = 0;

figure(1)
tiledlayout('flow', 'padding', 'tight');

for i = 1:nMic            % For each microphone signal
    % Load the signal
    fileName = strcat(dir, num2str(i), '.wav');  % i-th file name
    [y, Fs] = audioread(fileName);               % read i-th audio
    assert(Fs==fs)
    t_y = linspace(0,length(y)/Fs, length(y));
    
    % The window is applied to the signal from the 0 time instant.
    % ------------ Check the effect of differen window sizes. -------------
    y_w = y(1:length(t)).*w;
    
    % Plot the estimate impulse response
    %nexttile
    
    % Plot the mic signal
    hold on
    plot(t_y, y./max(y))
    % Plot the window over the signal
    plot(t, w, 'linewidth', 2)
    % Plot the windowed signal with thicker line
    plot(t, y_w./max(y), 'linewidth', 2)
    % Plot the TOA using stem (see doc stem)
    index = find(abs(t-TOA_directSignal(i))==min(abs(t-TOA_directSignal(i))));
    stem(t(index), y(index)./max(y))
    % Plot the first reflection TOA using stem (see doc stem)
    index = find(abs(t-TOA_firstReflection(i))==min(abs(t-TOA_firstReflection(i))));
    stem(t(index), y(index)./max(y))
    hold off
    if i == 1
        % Add a legend
        legend('signal', 'window', 'direct TOA', 'First reflection TOA', ...
            'windowed signal');
    end
    xlim([0 0.05]);             % Limit the plot btw 0 and 0.05s
    xlabel('Time (sec)');
    title(['Mic: ', num2str(i)]);
    
    sig(:,i) = y_w;         % Add current signal to the structure sig
end
%% Radiance estimation

SIG = fft(sig, nfft);  % FFT of the windowed signal

rad_patt = abs(R.*SIG);% Compute the radiance pattern magnitude

%% Radiance pattern plot
% Frequencies must be centered at a frequency bins
ctr_freqs = [150 250 500 750 1000 1500 2000 5000 7500 10000];

% Microphone angle wrt the speaker frontal face
angs = linspace(0, 360, nMic);

% Plot the estimated radiance pattern using the provided function
% radianceplot

radianceplot(ctr_freqs, rad_patt, angs, 'noise direct: ');
