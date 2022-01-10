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

fs = 48000;          % Sampling frequency
nfft = fs;           % Number of fft points
nMic = 24;           % Number of microphone signals

% Distance between source and microphones
R = [2.7647250, 2.7647250, 2.7718875, 2.7790501, 2.7862124, 2.8005376, ...
     2.8220251, 2.8435125, 2.8650000, 2.8864875, 2.9151375, 2.9294624, ...
     3.0225749, 2.9294624, 2.9079750, 2.8864875, 2.8650000, 2.8435125, ...
     2.8148625, 2.8005376, 2.7862124, 2.7790501, 2.7647250, 2.7647250];

typeOfSignal = 'noise/';                % Noise
dir = ['./recordings/' typeOfSignal];   % Signals directory

%% RADIANCE COMPUTATION THROUGH THE USE OF THE DIRECT RECORDED SIGNAL
% Window the signals according to the reflection time and estimate the
% radiance pattern. 
sig = [];   % Signal structure

% Early reflections attenuation. We consider only 2*ns+1 samples of the signal
ns = 350;       % Ideal ns should be on the direct path
t = (0:1/fs:((2*ns)+1)/fs);
t = t(1:end-1);

% TOA direct sound
TOA_directSignal = [0.0080416668, 0.0080416668, 0.0080624996, ...
    0.0080833333, 0.0081041669, 0.0081458334, 0.0082083335, 0.0082708336, ...
    0.0083333338, 0.0083958330, 0.0084791668, 0.0085208332, 0.0087916665, ...
    0.0085208332, 0.0084583331, 0.0083958330, 0.0083333338, 0.0082708336, ...
    0.0081874998, 0.0081458334, 0.0081041669, 0.0080833333, 0.0080416668, ...
    0.0080416668];
% TOA first reflection
TOA_firstReflection = [0.011208333, 0.011291667, 0.011250000, ...
    0.011291667, 0.011354167, 0.011395833, 0.011437500, 0.011458334, ...
    0.011541666, 0.011479166, 0.011895834, 0.011791667, 0.011854167, ...
    0.011666667, 0.011791667, 0.011604167, 0.011604167, 0.011395833, ...
    0.011333333, 0.011166667, 0.011270833, 0.011354167, 0.011250000, ...
    0.0059166667];

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
    t_y = linspace(0, length(y)/fs, length(y));
    
    % The window is applied to the signal from the 0 time instant.
    % ------------ Check the effect of different window sizes. ------------
    y_w = y(1:length(t)).*w;
    
    % Plot the estimate impulse response
    nexttile
    
    % Plot the mic signal
    hold on
    plot(t_y, y./max(y))
    % Plot the window over the signal
    plot(t, w, 'linewidth', 2)
    % Plot the windowed signal with thicker line
    plot(t, y_w./max(y), 'linewidth', 2)
    % Plot the TOA using stem (see doc stem)
    index = find(abs(t-TOA_directSignal(i))==min(abs(t-TOA_directSignal(i))));
    stem(t(index), y(index)./max(y), 'filled', 'linewidth', 2)
    % Plot the first reflection TOA using stem (see doc stem)
    index = find(abs(t-TOA_firstReflection(i))==min(abs(t-TOA_firstReflection(i))));
    stem(t(index), y(index)./max(y), 'filled', 'linewidth', 2)
    hold off
    if i == 1
        % Add a legend
        legend('signal', 'window', 'windowed signal', ...
            'direct TOA', 'First reflection TOA');
    end
    xlim([0 0.05]);             % Limit the plot btw 0 and 0.05s
    ylim([-1 1]);
    xlabel('Time (sec)');
    title(['Mic: ', num2str(i)]);
    
    sig(:,i) = y_w;         % Add current signal to the structure sig
end

%% Radiance estimation

SIG = fft(sig, nfft);   % FFT of the windowed signal

rad_patt = abs(R.*SIG); % Compute the radiance pattern magnitude

%% Radiance pattern plot

% Frequencies must be centered at a frequency bins
ctr_freqs = [150 250 500 750 1000 1500 2000 5000 7500 10000];

% Microphone angle wrt the speaker frontal face
angs = linspace(0, 360, nMic);

% Plot the estimated radiance pattern using the provided function
% radianceplot

radianceplot(ctr_freqs, rad_patt, angs, 'noise direct: ');

