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

nMic = % Number of microphones
typeOfSignal = % Noise or sweep
dir =       % Recordings directory

speed_of_sound = 343.8; % [m]/[s]

%% Plot the signal autocorrelation
figure;
t = tiledlayout('flow');

for n = 1:nMic
    % Load the signal
       % Time length
    % Time axis
    % Auto correlation
    Half of the autocorrelation
    % Plot the autocorrelation
    nexttile
    plot(t, xc);
    title(['Mic: ', num2str(n)]);
    axis([0 0.02 -1 1]);    % Limit the axis
    xlabel('Time (sec)');
end

%% MIC TO REFLECTORS DISTANCE COMPUTATION
% Put here the difference between first reflection 
delay = %[s]

distance = ;%[m]

fprintf(sprintf('Average distance between mic and first reflection %f m\n', distance));
