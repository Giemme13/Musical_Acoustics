%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework 4 exercise 5b
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

fs = 48000;         % Sampling frequency
nfft = fs;          % Number of fft points
duration = 10;      % [s] duration of sweep signal
nMic = 24;          % Number of microphones
% Distance between source and microphones
R = [2.7575624, 2.7647250, 2.7647250, 2.7718875, 2.7862124, 2.8005376, ...
     2.8864875, 2.9079750, 2.9366250, 2.9581125, 2.9796000, 3.0082500, ...
     3.0154126, 2.9939251, 2.9724374, 2.9509499, 2.9223001, 2.9008124, ...
     2.8148625, 2.8005376, 2.7862124, 2.7718875, 2.7647250, 2.7647250];

typeOfSignal ='sweep/'; % Sweep
dir = ['./recordings/', typeOfSignal]; % File directory

%% RADIANCE COMPUTATION THROUGH THE USE OF THE IMPULSE RESPONSE

sig = [];

% First reflections attenuation. We consider a small opart of the ir
% windowing the impulse response
ns = 50;        % We need a small window

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

% This is the input signal. The following function computes 1/FFT(sweep). 
% We will thus need to MULTIPLY this signal for FFT(output) in order to 
% apply deconvolution
% Use the provided synthSweep function with frequency interval between 50Hz
% and 22kHz
f1 = 50;
f2 = 22000;
[sweep, invsweepfft, sweepRate] = synthSweep(duration,fs,f1,f2);

figure;
tiledlayout('flow', 'padding', 'tight');


% Create a hanning window of 2*ns+1 length
    w = hann((2*ns)+1);

for n = 1:nMic                    % For each microphone signal
    % Load the signal
    fileName = strcat(dir, num2str(n), '.wav');  % i-th file name
    [y, Fs] = audioread(fileName);               % read i-th audio
    
    % Comput the impulse response using the function extractirsweep
    [ir] = extractirsweep(y, invsweepfft);
    t = linspace(0,length(ir)/fs, length(ir));
    
    % Windowing the ir. In order to take into account windows of arbitrary
    % length, the ir is circularly shifted until the max peak is at the
    % middle of the array. This prevents the left boundary of the windows
    % to be out of the vector. See function circshift.
    irCircular = circshift(ir', length(ir)/2);
    
    % Use the TOA information in order to find the correct poisition of the
    % first path in the shifted signal.
    index_toa = find(abs(t-TOA_directSignal(n))==min(abs(t-TOA_directSignal(n))));
    index_toa = index_toa + length(ir)/2;
    index_toa_refl = find(abs(t-TOA_firstReflection(n))==min(abs(t-TOA_firstReflection(n))));
    index_toa_refl = index_toa_refl + length(ir)/2;
    
    % Impulse response windowing
    % we want to maintain time information in the windowed impulse response.
    % We allocate one array of the same length of the original ir, and then
    % add samples from the windowed signal. Finally, we invert the
    % previously applied circShift to return back to the correct impulse
    % response timing.
    irCircular_w = zeros(1, length(ir)); % Initialize the vector
    irCircular_w(index_toa-ns:index_toa+ns) = w'.*irCircular(index_toa-ns:index_toa+ns); % Window the ir 
    ir_w = circshift(irCircular_w, length(ir)-index_toa+ns); % Shift ir back
    
    % Plot the estimated impulse response
    nexttile
    hold on
    plot(t, ir)
    % Plot the TOA and TOA First over the IR
    stem(TOA_directSignal(n), ir(index_toa));
    stem(TOA_firstReflection(n), ir(index_toa_refl));
    % Plot the windowed IR over with thicker line
    plot(t, ir_w)
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

SIG = fft(sig, nfft);    % FFT of the windowed signal

rad_patt = abs(R.*SIG);  % Compute the radiance pattern magnitude

%% Radiance pattern plot
% Frequencies must be centered at a frequency bins
ctr_freqs = [150 250 500 750 1000 1500 2000 5000 7500 10000];

% Microphone angle wrt the speaker frontal face
angs = linspace(0, 360, nMic);

% Plot the estimated radiance pattern using the provided function
% radianceplot

radianceplot(ctr_freqs, rad_patt, angs, 'sweep IR: ');

