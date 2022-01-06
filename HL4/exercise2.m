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

nMic = 1;                                  % Number of microphones
typeOfSignal = 'sweep/';                   % Noise or sweep
dir = ['./recordings/' typeOfSignal];      % Recordings directory

speed_of_sound = 343.8;                    % [m]/[s]

%% Plot the signal autocorrelation

signals = cell(1,24);               %stores signals
time_axis = cell(1,24);             %stores time axis of each signal
autocorrelations = cell(1,24);      %stores signals' autocorrelation
plot_just = 6;                      %plots diagrams of just one signal

for i = 1:24
    % Load the signal
    fileName = strcat(dir, num2str(i), '.wav');  % i-th file name
    [x, Fs] = audioread(fileName);               % read i-th audio
    signals{1,i} = x;
    % Time length
    time_length = length(x)/Fs;
    % Time axis
    t = linspace(0,time_length,length(x))';
    time_axis{1,i} = t;
    % Auto correlation
    xc = xcorr(x, 'normalized');
    % Half of the autocorrelation
    xc = xc(round((length(xc)/2)):end);
    autocorrelations{1,i} = xc;
    if plot_just == i
        % Plot 
        figure(1)
        subplot(2,1,1)
        plot(t, x)
        xlabel('Time [s]', 'fontsize', 20)
        title('Original signal', 'fontsize', 20)
        subplot(2,1,2)
        plot(t, xc);
        xlabel('Time [s]', 'fontsize', 20)
        title('Autocorrelation', 'fontsize', 20)
        sgtitle(['Mic: ', num2str(nMic)], 'fontsize', 20)
    end
end


%% MIC TO REFLECTORS DISTANCE COMPUTATION
% Put here the difference between first reflection
delay = zeros(1,24);
distance = zeros(1,24);

for i = 1:24
    xc = autocorrelations{1,i};
    t = time_axis{1,i};
    
    [max,locs] = findpeaks(abs(xc));
    if plot_just == i
        figure(2)
        plot(t, xc);
        xlim([0,0.03])
        xlabel('Time [s]', 'fontsize', 20)
        sgtitle('Autocorrelation', 'fontsize', 20)
    end
    delays = t(locs);
    delay(i) = delays(1);                           %[s]
    distance(i) = delay(i)*speed_of_sound;          %[m]
end

mean_delay = mean(delay)
mean_distance = mean(distance);


fprintf(sprintf('Average distance between mic and first reflection %f m\n', mean_distance));
