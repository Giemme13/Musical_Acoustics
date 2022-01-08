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
nMic = 24;           % Number of microphone signals

% Distance between source and microphones
R = [2.764725, 2.764725, 2.779050, 2.786213, 2.786213, 2.800538, ...
     2.822025, 2.843513, 2.865000, 2.886488, 2.915138, 2.929463, ...
     3.022575, 2.929463, 2.907975, 2.886488, 2.865000, 2.843513, ...
     2.814863, 2.800538, 2.786213, 2.779050, 2.764725, 2.764725];

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
TOA_directSignal = [0.0080417, 0.0080417, 0.0080833, 0.0081042, ...
    0.0081042, 0.0081458, 0.0082083, 0.0082708, 0.0083333, 0.00839583, ...
    0.0084792, 0.0085208, 0.0087917, 0.0085208, 0.0084583, 0.00839583, ...
    0.0083333, 0.0082708, 0.0081875, 0.0081458, 0.0081042, 0.00808333, ...
    0.0080417, 0.0080417];
% TOA first reflection
TOA_firstReflection = [0.0118750, 0.0111875, 0.0112083, 0.0112083, ...
    0.0114167, 0.0114167, 0.0114583, 0.0115417, 0.0115833, 0.0116458, ...
    0.0117708, 0.0118125, 0.0117917, 0.0118125, 0.0118333, 0.0116458, ...
    0.0113958, 0.0113542, 0.0113125, 0.0112708, 0.0111667, 0.0112083, ...
    0.0111875, 0.0111875];

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
    y = y./max(y);
    t_axis = linspace(0,length(y)/Fs, length(y));
    
    % The window is applied to the signal from the 0 time instant.
    % ------------ Check the effect of differen window sizes. -------------
    y_w = y(1:length(t));
    y_w = y_w.*w;
    
    % Plot the estimate impulse response
    Y = fft(y_w, nfft);
    [x, Fs_noise] = audioread('./input signals/noise.wav');
    assert(Fs_noise==fs)
    X = fft(x, nfft);
    H = Y./X;
    h = real(ifft(H, length(x)));
    t_h = linspace(0, length(h)/Fs_noise, length(h));
    nexttile
    plot(t_h, h)
    % Plot the mic signal
    hold on
    plot(t_axis, y)
    % Plot the window over the signal
    plot(t, w)
    % Plot the windowed signal with thicker line
    plot(t, y_w, 'linewidth', 2)
    % Plot the TOA using stem (see doc stem)
    index = find(abs(t-TOA_directSignal(i))==min(abs(t-TOA_directSignal(i))));
    stem(t(index), y(index), 'filled', 'linewidth', 2)
    % Plot the first reflection TOA using stem (see doc stem)
    index = find(abs(t-TOA_firstReflection(i))==min(abs(t-TOA_firstReflection(i))));
    stem(t(index), y(index), 'filled', 'linewidth', 2)
    hold off
    if i == 1
        % Add a legend
        legend('impulse response', 'signal', 'window', 'windowed signal', ...
            'direct TOA', 'First reflection TOA');
    end
    xlim([0 0.05]);             % Limit the plot btw 0 and 0.05s
    xlabel('Time (sec)');
    title(['Mic: ', num2str(i)]);
    
    sig(:,i) = y_w;         % Add current signal to the structure sig
end

%% Radiance estimation

SIG = fft(sig, nfft);              % FFT of the windowed signal

rad_patt = abs(R.*(SIG./X));       % Compute the radiance pattern magnitude

%% Radiance pattern plot

% Frequencies must be centered at a frequency bins
ctr_freqs = [150 250 500 750 1000 1500 2000 5000 7500 10000];

% Microphone angle wrt the speaker frontal face
angs = linspace(0, 360, nMic);

% Plot the estimated radiance pattern using the provided function
% radianceplot

radianceplot(ctr_freqs, rad_patt, angs, 'noise direct: ');

