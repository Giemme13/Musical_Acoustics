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
R = [2.764725, 2.764725, 2.779050, 2.786213, 2.786213, 2.800538, ...
     2.822025, 2.843513, 2.865000, 2.886488, 2.915138, 2.929463, ...
     3.022575, 2.929463, 2.907975, 2.886488, 2.865000, 2.843513, ...
     2.814863, 2.800538, 2.786213, 2.779050, 2.764725, 2.764725];

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

