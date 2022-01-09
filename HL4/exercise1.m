%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework 4 exercise 1
% Recording inspection
% We inspect the recordings for the estimation of the radiance pattern.
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

typeOfSignal = 'sweep/';        % Noise or sweep

dir = ['./recordings/' typeOfSignal];  % Recordings directory

nMic = 24;                      % Number of microphone signals
signalEnergy = zeros(1,nMic);   % Vector of energy

% Labels will be the angle 
labels = cell(1,nMic);
step = 360/nMic;                % one measure every microphone measure
deg = zeros(1,nMic);
for i = 1:nMic
    deg(i)=step*(i-1);
    labels{:,i} = strcat(num2str(deg(i)), '°');
end

%% Load the signals and compute the energy

for i = 1:nMic
    fileName = strcat(dir, num2str(i), '.wav');
    [x, Fs] = audioread(fileName);
    for j = 1:length(x)
        signalEnergy(i) = signalEnergy(i)+(abs(x(j)))^2;
    end    
end

% close the 'circle' placing at 360° the same value we have at 0°
signalEnergy = [signalEnergy, signalEnergy(1)];
deg = [deg, 360];

%% Plot the results

figure(1)
plot(deg, signalEnergy, 'linewidth', 2)
xlabel('Degrees [°]', 'fontsize', 20)
title('Energy of the signal', 'fontsize', 20)
xticks([0,45,90,135,180,225,270,315,360])
xlim([0,360])

theta = deg;
phi = linspace(0,0,25);
figure(2)
patternCustom(signalEnergy, theta, phi,'CoordinateSystem','polar',...
    'Slice','phi','SliceValue',0);
