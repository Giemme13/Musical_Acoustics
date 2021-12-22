clear all;
close all;
clc;

%% Data

c = 343;                %speed of sound in air [m/s]
rho = 1.225;            %air density [kg/m^3]
alpha = 0.75;           %cone semiangle [deg]
alpha = deg2rad(alpha); %cone semiangle [rad]
L0 = 0.45;              %length of the resonator [m]

%% 1) RESONATOR

%% a) Find head and foot diameters

F4 = 349.23;            %frequency with closed holes [Hz]
omega_F4 = 2*pi*F4;     %frequency [rad/s]
k_F4 = omega_F4/c;      %wavenumber [m^-1]
lambda_F4 = 2*pi/k_F4;  %wavelength [m]

%write the input impedance as a function of x1
%geometry: moving on x-axis -> vertex v=0, x1>v, x2>x1
x_1 = linspace(-10, 10, 10000);
x_2 = x_1+L0;
theta_1 = (1/k_F4)*atan(k_F4*x_1);
theta_2 = (1/k_F4)*atan(k_F4*x_2);
r_1 = x_1*tan(-alpha);
r_2 = x_2*tan(-alpha);
S_1 = pi*r_1.^2;
S_2 = pi*r_2.^2;

d_L = 0.6 * r_2;        %end correction
L_ = L0+d_L;

Z_IN = (1i*rho*c./S_2).*(sin(k_F4*L_).*sin(k_F4*theta_2))./sin(k_F4*(L_-theta_2));

%% Plot the impedance

figure()
subplot(2,1,1)
plot(x_1, db(abs(Z_IN)), 'linewidth', 1.5);
grid on
xlabel('$x_1\;[m]$', 'interpreter', 'latex', 'fontsize', 17)
ylabel('$|Z_{IN}|\;[dB]$', 'interpreter', 'latex', 'fontsize', 17)
subplot(2,1,2)
plot(x_1, angle(Z_IN), 'linewidth', 1.5);
grid on
xlabel('$x_1\;[m]$', 'interpreter', 'latex', 'fontsize', 17)
ylabel('$\angle Z_{IN}\;[dB]$', 'interpreter', 'latex', 'fontsize', 17)

%% By inspection

[min, locs] = findpeaks(db(abs(Z_IN)));
loc = locs(1);

x1 = x_1(loc);
x2 = x_2(loc);
theta1 = theta_1(loc);
theta2 = theta_2(loc);
S1 = S_1(loc);
S2 = S_2(loc);
dL = abs(d_L(loc));
L = L_(loc);
r1 = r_1(loc);
r2 = r_2(loc);

d1 = 2*r1;                  %diameter of the resonator at the head [m]
d2 = 2*r2;                  %diameter of the resonator at the foot [m]


%% Check on the actual input impedance

freq = linspace(1,2000,10000);
omega = 2*pi*freq;
k = omega/c;

Z_IN = (1i*rho*c/S2) * (sin(k*L).*sin(k*theta2))./sin(k*(L-theta2));

figure()
subplot(2,1,1)
plot(freq, db(abs(Z_IN)), 'linewidth', 1.5)
xlabel('$x_1\;[m]$', 'interpreter', 'latex', 'fontsize', 17)
ylabel('$|Z_{IN}|\;[dB]$', 'interpreter', 'latex', 'fontsize', 17)
subplot(2,1,2)
plot(freq, angle(Z_IN), 'linewidth', 1.5);
xlabel('$x_1\;[m]$', 'interpreter', 'latex', 'fontsize', 17)
ylabel('$\angle Z_{IN}\;[dB]$', 'interpreter', 'latex', 'fontsize', 17)


%% b) Position of the last finger hole

G4 = 392;               %frequency with last hole opened [Hz]
omega_G4 = 2*pi*G4;     %frequency [rad/s]
k_G4 = omega_G4/c;      %wavenumber [m^-1]
lambda_G4 = 2*pi/k_G4;  %wavelength [m]

%write the input impedance as a function of D
D = linspace(0,L0,1000);        %distance of the hole from the foot
delta = D + dL^2./(D+2*dL);     %difference between equivalent resonators

L_G4 = L - delta;               %length of a the pipe if the hole is opened

Z_IN = (1i*rho*c./S2).*(sin(k_G4*L_G4).*sin(k_G4*theta2))./sin(k_G4*(L_G4-theta2));

%% Plot the impedance

figure()
subplot(2,1,1)
plot(D, db(abs(Z_IN)), 'linewidth', 1.5);
xlabel('$x_1\;[m]$', 'interpreter', 'latex', 'fontsize', 17)
ylabel('$|Z_{IN}|\;[dB]$', 'interpreter', 'latex', 'fontsize', 17)
subplot(2,1,2)
plot(D, angle(Z_IN), 'linewidth', 1.5);
xlabel('$x_1\;[m]$', 'interpreter', 'latex', 'fontsize', 17)
ylabel('$\angle Z_{IN}\;[dB]$', 'interpreter', 'latex', 'fontsize', 17)

%% By inspection

[min, locs] = findpeaks(db(abs(Z_IN)));
loc = locs(1);
D_G4 = D(loc);                  %distance of the hole from the foot
x_G4 = L0 - D_G4;               %distance of the hole from the head
delta_G4 = delta(loc);
dL_1 = D_G4 + dL - delta_G4;    %end correction of the equivalent resonator


%% c) Position of the second to last finger hole

A4 = 440;               %frequency with last hole opened [Hz]
omega_A4 = 2*pi*A4;     %frequency [rad/s]
k_A4 = omega_A4/c;      %wavenumber [m^-1]
lambda_A4 = 2*pi/k_A4;  %wavelength [m]

%write the input impedance as a function of D
D = linspace(0,x_G4,1000);        %distance of the hole from the foot
delta = D + (dL_1)^2./(D+2*dL_1); %difference between equivalent resonators

L_A4 = L - delta_G4 - delta;   %length of a the pipe if the two holes are opened

Z_IN = (1i*rho*c./S2).*(sin(k_A4*L_A4).*sin(k_A4*theta2))./sin(k_A4*(L_A4-theta2));

%% Plot the impedance

figure()
subplot(2,1,1)
plot(D, db(abs(Z_IN)), 'linewidth', 1.5);
xlabel('$x_1\;[m]$', 'interpreter', 'latex', 'fontsize', 17)
ylabel('$|Z_{IN}|\;[dB]$', 'interpreter', 'latex', 'fontsize', 17)
subplot(2,1,2)
plot(D, angle(Z_IN), 'linewidth', 1.5);
xlabel('$x_1\;[m]$', 'interpreter', 'latex', 'fontsize', 17)
ylabel('$\angle Z_{IN}\;[dB]$', 'interpreter', 'latex', 'fontsize', 17)

%% By inspection

[min, locs] = findpeaks(db(abs(Z_IN)));
loc = locs(1);
D_2_prime = D(loc);            %distance between holes
D_A4 = D_2_prime + D_G4;        %distance of the hole from the foot
x_A4 = L0 - D_A4;               %distance of the hole from the head
delta_A4 = delta(loc);
dL_2 = D_2_prime + dL_1 - delta_A4;  %end correction of the equivalent resonator


%% Plot the flute
%we already have the diameters at head and foot -> d1, d2
r_head = d1/2; r_foot = d2/2;
%let's place the input at coord x=0
x1 = 0;
x2 = L0;
%the last hole is positioned at x_G4 and has section (pi*r2^2)
x_lh_right = x_G4+r2;
x_lh_left = x_G4-r2;
r_lh_right = (x2-x_lh_right)*tan(alpha)+r2;
r_lh_left = (x2-x_lh_left)*tan(alpha)+r2;
%the second to last hole is positioned at x_A4 and has section (pi*r2^2)
x_slh_right = x_A4+r2;
x_slh_left = x_A4-r2;
r_slh_right = (x2-x_slh_right)*tan(alpha)+r2;
r_slh_left = (x2-x_slh_left)*tan(alpha)+r2;

figure()
hold on
yline(0, '-.')
%line([0,L0], [r1,r2], 'color', 'k')
line([0,x_slh_left], [r1,r_slh_left], 'linewidth', 1.5)
line([x_slh_right,x_lh_left], [r_slh_right,r_lh_left], 'linewidth', 1.5)
line([x_lh_right,L0], [r_lh_right,r2], 'linewidth', 1.5)
line([0,L0], [-r1,-r2], 'linewidth', 1.5)
hold off
grid on
xlim([0,L0])
ylim([-L0/3,L0/3])