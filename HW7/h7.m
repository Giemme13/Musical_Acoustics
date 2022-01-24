
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
figure(1)
plot(x,z)
axis equal

%% Inharmonicity

a = 0.1:0.1:1.2;

a_values = length(a);
f_values = 15;  %this must correspond to the searched frequencies in COMSOL

% f = matrix [12x5] of first 5 eigenfrequencies for each value of a
f = zeros(12,5);

fileID = fopen('./DataFromComsol.txt');

i = 1;
while ~feof(fileID)
    line = fgetl(fileID);
    if (line(1)) ~= '%'
        text_lines{i} = line;
        i = i+1;
    end
end
fclose(fileID);

for i = 1:a_values
    for j = 7:11    %taking only the searched 5 frequencies each value of a
        for k = 51:60  %taking only the correct characters inside each line
            freq(k) = text_lines{1,(i-1)*f_values+j}(k);
        end
        f(i,j-6) = str2num(freq);
    end
end


m = 1:1:5;
%matrix containing the inharmonicities of each harmonic with respect to the
%corresponding lower eigenfrequencies, for each value of a
I = zeros(12, 4);

for i = 1:length(a) %loop all values of parameter a
    for N = 2:5 %loop on harmonics
        for n = 2:N %loop on eigenfrequencies lower than N
            [value, m_n] = min(abs(f(i,n)-m*f(i,n-1)));
            I(i,N-1) = I(i,N-1) + abs((f(i,n)/f(i,n-1)) - m_n);
        end
    end
end

%% Plot

a = 0.1:0.1:1.2;

line1 = I(:,1);
line2 = I(:,2);
line3 = I(:,3);
line4 = I(:,4);

figure(2)
tightfig;
hold on
plot(a, line1, '-o')
plot(a, line2, '-o')
plot(a, line3, '-o')
plot(a, line4, '-o')
hold off
grid on
title('Inharmonicity', 'fontsize', 20)
xlabel('$a$', 'interpreter', 'latex', 'fontsize', 18)
ylabel('$I$', 'interpreter', 'latex', 'fontsize', 18)
legend('1st harmonic', '2nd harmonic', '3rd harmonic', '4th harmonic', ...
    'fontsize', 15, 'location', 'northwest')

