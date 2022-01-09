%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework 4 exercise 2
% Reflection inspection
% By a visual inspection of the recorded signals we evaluate the first
% reflection time from which the distance of the micrphones from the
% reflectors can be inferred.
%
% Musical Acoustic Course
% Mirco Pezzoli
% 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc 

%% First reflection estimation using autocorrelation
% Collecting autocorrelations and plotting them in order to analyze
% reflections. Try different takes and input signals

addpath('Functions');

nMic = 24;                                 % Number of microphone signals
typeOfSignal = 'sweep/';                   % Noise or sweep
dir = ['./recordings/' typeOfSignal];      % Recordings directory

speed_of_sound = 343.8;                    % [m]/[s]

%% Plot the signal autocorrelation

autocorrelations = cell(1, nMic);
time_axis = cell(1,nMic);

figure(1)
tiledlayout('flow', 'padding', 'tight');
for i = 1:nMic
    % Load the signal
    fileName = strcat(dir, num2str(i), '.wav');
    [x, Fs] = audioread(fileName);
    % Time length
    time_length = length(x)/Fs;
    % Time axis
    t = linspace(0,time_length,length(x));
    time_axis{1, i} = t;
    % Auto correlation
    xc = xcorr(x, 'normalized');
    xc = xc(round((length(xc)/2)):end);
    autocorrelations{1, i} = xc;
    % Plot the autocorrelation
    nexttile
    plot(t, xc);
    title(['Mic: ', num2str(i)]);
    axis([0 0.02 -1 1]);    % Limit the axis
    xlabel('Time (sec)');
end


%% MIC TO REFLECTORS DISTANCE COMPUTATION
% Put here the difference between first reflection 
% values taken from inspection of autocorrelation diagrams
delay_noise = [0.00316668, 0.00316668, 0.00314585, 0.00314585, 0.00314585, ...
    0.00306252, 0.00304168, 0.00304168, 0.00302085, 0.00306252, 0.00302085, ...
    0.00306252, 0.00300002, 0.00304168, 0.00306252, 0.00302085, 0.00302085, ...
    0.00308335, 0.00310418, 0.00312502, 0.00312502, 0.00314585, 0.00314585, ...
    0.00314585];    %[s]
delay_sweep = [0.00658335, 0.00668751, 0.00662501, 0.00679168, 0.00679168, ...
    0.00695835, 0.00656251, 0.00650001, 0.00641668, 0.00652085, 0.00645835, ...
    0.00664585, 0.00672918, 0.00664585, 0.00683335, 0.00691668, 0.00641668, ...
    0.00639584, 0.00664585, 0.00683335, 0.00687501, 0.00677084, 0.00662501, ...
    0.00662501];    %[s]

distance_noise = delay_noise * speed_of_sound;      %[m]
distance_sweep = delay_sweep * speed_of_sound;      %[m]

fprintf(sprintf('Average distance between mic and first reflection (noise) %f m\n',...
    mean(distance_noise)));
fprintf(sprintf('Average distance between mic and first reflection (sweep) %f m\n',...
    mean(distance_sweep)));
