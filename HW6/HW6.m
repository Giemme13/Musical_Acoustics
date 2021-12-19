clear;
close all;
clc;

%% Data

c = 343;          %speed of sound in air [m/s]
rho = 1.225;      %air density [kg/m^3]
alpha = 0.75;     %cone semiangle [deg]
L0 = 0.45;        %length of the resonator [m]

%% 1) RESONATOR

%% a) Find head and foot diameters

F4 = 349.23;            %frequency with closed holes [Hz]
omega_F4 = 2*pi*F4;     %frequency [rad/s]
k_F4 = omega_F4/c;      %wavenumber [m^-1]
lambda_F4 = 2*pi/k_F4;  %wavelength [m]

%write the input impedance as a function of x1
x_axis = linspace(0,10,1000);

theta1 = atan(k_F4*x_axis)/k_F4;
S1 = pi*(x_axis*tand(alpha)).^2;
dL = 8/(3*pi) * (x_axis+L0)*tand(alpha);    %end correction
L = L0+dL;

Z_IN = (1i*rho*c./S1).*(sin(k_F4*L).*sin(k_F4*theta1))./sin(k_F4*(L+theta1));

%% Plot the impedance

figure()
subplot(2,1,1)
plot(x_axis, db(abs(Z_IN)), 'linewidth', 1.5);
xlabel('$x_1\;[m]$', 'interpreter', 'latex', 'fontsize', 17)
ylabel('$|Z_{IN}|\;[dB]$', 'interpreter', 'latex', 'fontsize', 17)
subplot(2,1,2)
plot(x_axis, angle(Z_IN), 'linewidth', 1.5);
xlabel('$x_1\;[m]$', 'interpreter', 'latex', 'fontsize', 17)
ylabel('$\angle Z_{IN}\;[dB]$', 'interpreter', 'latex', 'fontsize', 17)

%% By inspection

[min, locs] = findpeaks(-db(abs(Z_IN)));
x1 = x_axis(locs(2));

theta1 = theta1(locs(2));
S1 = S1(locs(2));
dL = dL(locs(2));
L = L(locs(2));

d1 = 2*x1*tand(alpha);           %diameter of the resonator at the head [m]
d2 = 2*(x1+L0)*tand(alpha);      %diameter of the resonator at the foot [m]


%% b) Position of the last finger hole

G4 = 392;               %frequency with last hole opened [Hz]
omega_G4 = 2*pi*G4;     %frequency [rad/s]
k_G4 = omega_G4/c;      %wavenumber [m^-1]
lambda_G4 = 2*pi/k_G4;  %wavelength [m]

%write the input impedance as a function of D
D = linspace(0,L0,1000);
delta = D + dL^2./(D+dL);

L_G4 = L-delta;

Z_IN = (1i*rho*c./S1).*(sin(k_G4*L_G4).*sin(k_G4*theta1))./sin(k_G4*(L_G4+theta1));

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

[min, locs] = findpeaks(-db(abs(Z_IN)));
D_G4 = D(locs);
x_G4 = L0 - D_G4;


%% c) Position of the second to last finger hole

A4 = 440;               %frequency with last hole opened [Hz]
omega_A4 = 2*pi*A4;     %frequency [rad/s]
k_A4 = omega_A4/c;      %wavenumber [m^-1]
lambda_A4 = 2*pi/k_A4;  %wavelength [m]

%write the input impedance as a function of D
D = linspace(0,L0,1000);
delta = D + (dL+D_G4)^2./(D+dL+D_G4);

L_A4 = L-delta;

Z_IN = (1i*rho*c./S1).*(sin(k_A4*L_A4).*sin(k_A4*theta1))./sin(k_A4*(L_A4+theta1));

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

[min, locs] = findpeaks(-db(abs(Z_IN)));
D_A4 = D(locs);
x_A4 = L0 - D_A4;

