%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework 4 exercise 5a
% Radiance pattern estimation from impulse response
% We compute the impulse response using the sine sweep source signal
% from which the radiance pattern is estimated windowing the signal before
% the occurence of the first reflection
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

fs = 48000;             % Sampling frequency
nfft = fs;              % Number of fft points
speed_of_sound = 343.8; % [m]/[s]
nMic = 24;              % Number of microphones
% Distance between source and microphones
R = [2.7647250, 2.7647250, 2.7718875, 2.7790501, 2.7862124, 2.8005376, ...
     2.8220251, 2.8435125, 2.8650000, 2.8864875, 2.9151375, 2.9294624, ...
     3.0225749, 2.9294624, 2.9079750, 2.8864875, 2.8650000, 2.8435125, ...
     2.8148625, 2.8005376, 2.7862124, 2.7790501, 2.7647250, 2.7647250];

typeOfSignal = 'noise/'; % Noise
dir = ['./recordings/' typeOfSignal];% File directory

inputSignalDir = './input signals/';            % Source signal directory
inputSignalFileName = strcat(inputSignalDir, 'noise.wav');    % Source signal name

%% Radiance estimation using the impulse response

sig = [];

% First reflections attenuation. We consider a small part of the ir
% windowing the impulse response
ns = 50;        % We need a small window     

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

% Source signal must be known. Load the source signal.
[x, fs] = audioread(inputSignalFileName);

figure(1)
tiledlayout('flow', 'padding', 'tight');

sig = zeros(fs, nMic);

% Create a hanning window of 2*ns+1 length
w = hann((2*ns)+1);

for n = 1:nMic              % For each microphone signal
    % Load the signal
    fileName = strcat(dir, num2str(n), '.wav');  % i-th file name
    [y, Fs] = audioread(fileName);               % read i-th audio
      
    % Comput the impulse response using the function extractirnoise
    [ir] = extractirnoise(x, y, nfft);
    ir = abs(ir);
    t = linspace(0, length(ir)/fs, length(ir));

    % Windowing the ir. In order to take into account windows of arbitrary
    % length, the ir is circularly shifted until the first peak is at the
    % middle of the array. This prevents the left boundary of the windows
    % to be out of the vector. See function circshift.
    irCircular = circshift(ir', length(ir)/2); 
    
    % Use the TOA information in order to find the correct poisition of the
    % first path in the shifted signal.
    [peak, loc] = max(irCircular); %find TOA index
    
    % Impulse response windowing
    % we want to maintain time information in the windowed impulse response.
    % We allocate one array of the same length of the original ir, and then
    % add samples from the windowed signal. Finally, we invert the
    % previously applied circShift to return back to the correct impulse
    % response timing.
    irCircular_w = zeros(1, length(ir)); % Initialize the vector
    irCircular_w(loc-ns:loc+ns) = w'.*irCircular(loc-ns:loc+ns); % Window the ir
    ir_w = circshift(irCircular_w, length(ir)-loc+ns); % Shift ir back
    
    % Plot the estimated impulse response
    nexttile
    hold on
    plot(t, ir);
    % Plot the TOA and TOA First over the IR
    stem(TOA_directSignal(n), abs(max(ir_w)));
    stem(TOA_firstReflection(n), abs(max(ir_w)));
    % Plot the windowed IR over with thicker line
    plot(t, ir_w, 'linewidth', 2);
    hold off
    
    if n == 1
        legend('original ir', 'direct TOA', 'First reflection TOA', 'windowed ir');
    end
    
    xlim([0 0.05]);
    xlabel('Time (sec)');
    title(['Mic: ', num2str(n)]);
    
    sig(:,n) = ir_w;
end

%% Radiance estimation

SIG = fft(sig, nfft);     % FFT of the windowed signal

rad_patt = abs(R.*SIG);   % Compute the radiance pattern magnitude

%% Radiance pattern plot

% Frequencies must be centered at a frequency bins
ctr_freqs = [150 250 500 750 1000 1500 2000 5000 7500 10000];

% Microphone angle wrt the speaker frontal face
angs = linspace(0, 360, nMic);

% Plot the estimated radiance pattern using the provided function
% radianceplot

radianceplot(ctr_freqs, rad_patt, angs, 'noise IR: ');

