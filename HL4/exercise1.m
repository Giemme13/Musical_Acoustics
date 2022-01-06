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

nMic = 1;                       % Number of microphones
signalEnergy = zeros(1,24);     % Vector of energy

% Labels will be the angle 
labels = cell(1,24);
step = 15;                      % one measure every 15 degrees
deg = zeros(1,24);
for i = 1:24
   deg(i)=step*(i-1);
   labels{:,i} = strcat(num2str(deg(i)), '°');
end

%% Load the signals and compute the energy

for i = 1:24
   fileName = strcat(dir, num2str(i), '.wav');  % i-th file name
   [x, Fs] = audioread(fileName);               % read i-th audio
   for j = 1:length(x)
       signalEnergy(i) = signalEnergy(i)+(abs(x(j)))^2; % compute energy for i-th signal
   end    
end
%% Plot the results

figure(1)
plot(deg, signalEnergy, 'linewidth', 2)
xlabel('Degrees [°]', 'fontsize', 20)
title('Energy of the signal', 'fontsize', 20)
xticks([0,45,90,135,180,225,270,315,360])
xlim([0,360])



