%% Script for generating the thickness profile of a marimba bar as a function of
%the horizontal position
% Fabio Antonacci, 2021 @ Musical Acoustics, A.A. 2021/2022
dx = 0.01; % Step of the horizontal axis
x = 0:dx:10-dx; % Horizontal axis
N = length(x);
z = zeros(size(x)); % Variable for storing the thickness of the bar as a
%function of the horizontal position
a = 1.23; % To be modified in the range 0.1-->0.9, step 0.1
b = 1;
for n = 1:floor(N/8)-1
z(n) = a-2.5-a*cos((x(n)-1.25)*2*pi*b/5);
end
for n = floor(N/8):floor(N/4) % Constant thickness for the first fourth of the
%length of the bar
 z(n) = -2.5;
end
for n = floor(N/4)+1:3*floor(N/4) % Sinusoidal thickness profile in the central
%part of the bar
 z(n) = a-2.5+a*cos((x(n)-5)*2*pi*b/5);
end
for n = 6*floor(N/8)+1:7*floor(N/8) % Constant thickness for the last fourth of
%the length of the bar
 z(n) = -2.5;
end
for n = 7*floor(N/8)+1:N
z(n) = a-2.5-a*cos(-(x(n)-8.75)*2*pi*b/5);
end
plot(x,z)
axis equal

%%
a1=0.0123;

