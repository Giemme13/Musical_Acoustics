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
duration = 10;       % [s] duration of sweep signal
nMic = ;          % Number of microphones
R = ;            % Distance between source and microphones

typeOfSignal = % Sweep
dir =  % File directory

%% RADIANCE COMPUTATION THROUGH THE USE OF THE IMPULSE RESPONSE

sig = [];

% First reflections attenuation. We consider a small opart of the ir
% windowing the impulse response
ns = 50;        % We need a small window

TOA_directSignal = % Direct signal Time Of Arrival
TOA_firstReflection = 

% This is the input signal. The following function computes 1/FFT(sweep). 
% We will thus need to MULTIPLY this signal for FFT(output) in order to 
% apply deconvolution
% Use the provided synthSweep function with frequency interval between 50Hz
% and 22kHz
figure;
tiledlayout('flow');
sig = zeros(fs, nMic);
for n = 1:nMic                    % For each microphone signal
    % Load the signal
    
    % Create a hanning window of 2*ns+1 length
    w = 
    
    % Comput the impulse response using the function extractirsweep
    [ir] = 
    t = % Time axis
        
    % Windowing the ir. In order to take into account windows of arbitrary
    % length, the ir is circularly shifted until the max peak is at the
    % middle of the array. This prevents the left boundary of the windows
    % to be out of the vector. See function circshift.
    irCircular = 
    
    % Use the TOA information in order to find the correct poisition of the
    % first path in the shifted signal.
    
    % Impulse response windowing
    % we want to maintain time information in the windowed impulse response.
    % We allocate one array of the same length of the original ir, and then
    % add samples from the windowed signal. Finally, we invert the
    % previously applied circShift to return back to the correct impulse
    % response timing.
    irCircular_w = % Initialize the vector
    irCircular_w % Window the ir 
    ir_w =  % Shift ir back
    
   % Plot the estimated impulse response
    hold on
    % Plot the TOA and TOA First over the IR
   
   
    % Plot the windowed IR over with thicker line
    nexttile
   
    hold off
    legend('windowed ir');

%     legend('original ir', 'direct TOA', 'First reflection TOA', 'windowed ir');
    xlim([0 0.05]);
    xlabel('Time (sec)');
    title(['Mic: ', num2str(n)]);
    
    sig(:,n) = 
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

radianceplot(ctr_freqs, rad_patt, angs, 'sweep IR: ');

