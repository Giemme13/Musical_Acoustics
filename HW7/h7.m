
close all
clear
clc


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

%% Inharmonicity

a1 = 0.1:0.1:1.2;

%matrix [12x5] of first 5 eigenfrequencies for each value of
%parameter a. Copy from COMSOL results
f=[8448.14568, 8839.68910, 12557.46521, 20013.34040, 20068.92623;
    7907.51943, 8749.23353, 12068.14288, 19633.61191, 19757.32693;
    7283.16410, 8623.21846, 11478.21556, 19107.43680, 19347.97367;
    6587.74544, 8459.82585, 10785.81062, 18508.85673, 18768.58148;
    5836.37006, 8256.43438,  9991.69984, 17851.28629, 18007.22886;
    5043.81943, 8008.06239,  9097.45666, 17053.93612, 17137.14158;
    4225.47694, 7706.66028,  8106.90850, 15903.28736, 16356.82223;
    3394.27400, 7020.85481,  7337.63313, 14539.14271, 15482.89802;
    2565.63223, 5841.53745,  6877.11348, 12935.93688, 14468.12186;
    1750.28042, 4559.18915,  6275.61646, 11010.66347, 13214.74176;
     971.22527, 3153.43876,  5427.67839,  8569.57419, 11524.49953;
     268.33845, 1495.38814,  3955.68088,  4842.88208,  8662.59627];  

m = 0:1:5;
%matrix containing the inharmonicities of each harmonic with respect to the
%corresponding lower eigenfrequencies, for each value of a
I = zeros(12, 4);

for i = 1:length(a1) %loop all values of parameter a
    for N = 2:5 %loop on harmonics
        for n = 2:N %loop on eigenfrequencies lower than N
            [value, m_n] = min(abs(f(i,n)-m*f(i,n-1)));
            I(i,N-1) = I(i,N-1) + abs((f(i,n)/f(i,n-1)) - m_n);
        end
    end
end

%% Plot

line1 = I(:,1);
line2 = I(:,2);
line3 = I(:,3);
line4 = I(:,4);

figure(1)
hold on
plot(a1, line1, '-o')
plot(a1, line2, '-o')
plot(a1, line3, '-o')
plot(a1, line4, '-o')
hold off
grid on

