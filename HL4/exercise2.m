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

nMic = 24;                                 % Number of microphone signals
typeOfSignal = 'noise/';                   % Noise or sweep
dir = ['./recordings/' typeOfSignal];      % Recordings directory

speed_of_sound = 343.8;                    % [m]/[s]

%% Plot the signal autocorrelation

autocorrelations = cell(1, nMic);
time_axis = cell(1,nMic);
envHigh = cell(1,nMic);

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
    [envHigh{1,i}, envLow] = envelope(xc, 10, 'peak');
    % Plot the autocorrelation
    nexttile
    hold on
    plot(t, xc);
    plot(t, envHigh{1,i});
    hold off
    title(['Mic: ', num2str(i)]);
    axis([0 0.02 -1 1]);    % Limit the axis
    xlabel('Time (sec)');
    if i == 1
        legend('autocorrelation', 'envelope')
    end
end


%% MIC TO REFLECTORS DISTANCE COMPUTATION
% Put here the difference between first reflection 

delays = zeros(1,nMic);

for i = 1:nMic
    xc_env = envHigh{1,i};
    t = time_axis{1,i};
    [peaks,locs] = findpeaks(xc_env);
    delays(i) = t(locs(1));
end

distance = delays * speed_of_sound;

fprintf(sprintf('Average distance between mic and first reflection %f m\n',...
    mean(distance)));
