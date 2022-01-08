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
addpath('Functions');

fs = 48000;          % Sampling frequency
nfft = fs;           % Number of fft points
nMic = 24;           % Number of microphone signals

% Distance between source and microphones
R = [2.7575625, 2.7647250, 2.7647250, 2.7718875, 2.7862125, 2.8005375, ...
     2.8864875, 2.9079750, 2.9366250, 2.9581125, 2.9796000, 3.0082500, ...
     3.0154125, 2.9939250, 2.9724375, 2.9509500, 2.9223000, 2.9008125, ...
     2.8148625, 2.8005375, 2.7862125, 2.7718875, 2.7647250, 2.7647250];

typeOfSignal = 'sweep/';                % Sweep
dir = ['./recordings/' typeOfSignal];   % Signals directory

% data usefull for sine sweep synthesis
duration = 10;
f1 = 50;
f2 = 22000;
%% RADIANCE COMPUTATION THROUGH THE USE OF THE DIRECT RECORDED SIGNAL
% Window the signals according to the reflection time and estimate the
% radiance pattern. 
sig = [];   % Signal structure

[sweep, invsweepfft, sweepRate] = synthSweep(duration,fs,f1,f2);

% Early reflections attenuation. We consider only 2*ns+1 samples of the signal
ns = 350; % Ideal ns should be on the direct path
t = (0:1/fs:((2*ns)+1)/fs);
t = t(1:end-1);

% TOA direct sound
TOA_directSignal =  [0.0080208, 0.0080417, 0.0080417, 0.0080625, ...
    0.0081041, 0.0081458, 0.0083958, 0.0084583, 0.0085417, 0.0086042, ...
    0.0086667, 0.0087500, 0.0087708, 0.0087083, 0.0086458, 0.0085833, ...
    0.0085000, 0.0084375, 0.0081875, 0.0081458, 0.0081042, 0.0080625, ...
    0.0080417, 0.0080417];
% TOA first reflection
TOA_firstReflection = [0.0111250, 0.0111875, 0.0111875, 0.0112083, ...
    0.0111875, 0.0112708, 0.0113125, 0.0113542, 0.0115833, 0.0116250, ...
    0.0117292, 0.0117708, 0.0117917, 0.0117708, 0.0116875, 0.0116458, ...
    0.0113750, 0.0113333, 0.0112917, 0.0112500, 0.0112292, 0.0112083, ...
    0.0111875, 0.0111250];

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
    [ir] = extractirsweep(y_w, invsweepfft);
    t_ir = linspace(0, length(ir)/fs, length(ir));
%     [x, Fs_noise] = audioread('./input signals/noise.wav');
%     assert(Fs_noise==fs)
%     X = fft(x, nfft);
%     H = Y./X;
%     h = real(ifft(H, length(x)));
%     t_h = linspace(0, length(h)/Fs_noise, length(h));
    nexttile
    plot(t_ir, ir)
    % Plot the mic signal
    hold on
    plot(t_axis, y)
    % Plot the window over the signal
    plot(t, w)
    % Plot the windowed signal with thicker line
    plot(t, y_w, 'linewidth', 2)
    % Plot the TOA using stem (see doc stem)
    index = find(abs(t-TOA_directSignal(i))==min(abs(t-TOA_directSignal(i))));
    stem(t(index), y(index))
    % Plot the first reflection TOA using stem (see doc stem)
    index = find(abs(t-TOA_firstReflection(i))==min(abs(t-TOA_firstReflection(i))));
    stem(t(index), y(index))
    hold off
    if i == 1
        % Add a legend
        legend('impulse response', 'signal', 'window', 'direct TOA', 'First reflection TOA', ...
            'windowed signal');
    end
    xlim([0 0.05]);             % Limit the plot btw 0 and 0.05s
    xlabel('Time (sec)');
    title(['Mic: ', num2str(i)]);
    
    sig(:,i) = y_w;         % Add current signal to the structure sig
end
%% Radiance estimation

SIG = fft(sig, nfft);              % FFT of the windowed signal

rad_patt = abs(R.*(SIG./(fft(sweep, nfft)')));       % Compute the radiance pattern magnitude

%% Radiance pattern plot
% Frequencies must be centered at a frequency bins
ctr_freqs = [150 250 500 750 1000 1500 2000 5000 7500 10000];

% Microphone angle wrt the speaker frontal face
angs = linspace(0, 360, nMic);

% Plot the estimated radiance pattern using the provided function
% radianceplot

radianceplot(ctr_freqs, rad_patt, angs, 'noise direct: ');
