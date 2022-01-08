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
nMic = 1;            % Number of microphones

R = 2.57;            % Distance between source and microphones

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

TOA_directSignal =  0.0083;      % TOA
TOA_firstReflection = 0.0114;  % First reflection TOA

w = hann(length(t));      % Window

figure;
tiledlayout('flow');

for n = 1:nMic            % For each microphone signal
    % Load the signal
    fileName = strcat(dir, num2str(n), '.wav');  % i-th file name
    [y, Fs] = audioread(fileName);               % read i-th audio
    assert(Fs==fs)
    y = y./max(y);
    t_axis = linspace(0,length(y)/Fs, length(y));
    
    % The window is applied to the signal from the 0 time instant.
    % ------------ Check the effect of differen window sizes. -------------
    y_w = y(1:length(t));
    y_w = y_w.*w;
    
    if plot_just == n
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
        % Plot the TOA using stem (see doc stem)
        index = find(abs(t-TOA_directSignal(n))==min(abs(t-TOA_directSignal(n))));
        stem(t(index), y(index))
        % Plot the first reflection TOA using stem (see doc stem)
        index = find(abs(t-TOA_firstReflection(n))==min(abs(t-TOA_firstReflection(n))));
        stem(t(index), y(index))
        % Plot the windowed signal with thicker line
        plot(t, y_w, 'linewidth', 2)
        hold off
        % Add a legend
        legend('impulse response', 'signal', 'window', 'direct TOA', 'First reflection TOA', ...
            'windowed signal');
        xlim([0 0.05]);             % Limit the plot btw 0 and 0.05s
        xlabel('Time (sec)');
        title(['Mic: ', num2str(n)]);
    end
    sig(:,n) = y_w;         % Add current signal to the structure sig
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
