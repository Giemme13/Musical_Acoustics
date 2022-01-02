clear all;
close all;
clc

%% CONSTANTS

c = 343;                %speed of sound in air [m/s]
rho = 1.225;            %air density [kg/m^3]
nu = 1.5e-5;            %kinematic viscosity [m^2/s]

%% 1) RESONATOR

alpha = 0.75;           %cone semiangle [deg]
alpha = deg2rad(alpha); %cone semiangle [rad]
L0 = 0.45;              %length of the whole instrument (and resonator) [m]

%% a) Find head and foot diameters

F4 = 349.23;            %frequency with closed holes [Hz]
omega_F4 = 2*pi*F4;     %frequency [rad/s]
k_F4 = omega_F4/c;      %wavenumber [m^-1]

%write the input impedance as a function of x1
%geometry: vertex on x=0, x_1>0, x_2>x_1, L>0
%narrow end
r_1 = linspace(0,0.05,10000);
x_1 = r_1/tan(alpha);
theta_1 = atan(k_F4*x_1)/k_F4;
S_1 = pi*r_1.^2;
%wide end
x_2 = x_1 + L0;
r_2 = x_2*tan(alpha);
theta_2 = atan(k_F4*x_2)/k_F4;
S_2 = pi*r_2.^2;

dL_ = 0.61*r_1;            %end correction
L_ = L0 + dL_;

%Input impedance as in Fletcher's book, page 451
Delta_L = 0.04;                 %typical value for alto recorder
M_ = (rho*Delta_L)./S_2;        %inertance of the mouth window

Z_IN = ((1i*rho*c)./S_2) .* ((sin(k_F4*L_).*sin(k_F4*theta_1))./sin(k_F4*(L_+theta_1))) + 1i*omega_F4*M_;


%% Plot the impedance

[min, locs] = findpeaks(-db(abs(Z_IN)));

figure()
subplot(2,1,1)
hold on
plot(r_1, db(abs(Z_IN)), 'linewidth', 1.5);
plot(r_1(locs), -min, 'or', 'linewidth', 2);
grid on
xlabel('$x_1\;[m]$', 'interpreter', 'latex', 'fontsize', 17)
ylabel('$|Z_{IN}|\;[dB]$', 'interpreter', 'latex', 'fontsize', 17)
subplot(2,1,2)
hold on
plot(r_1, angle(Z_IN), 'linewidth', 1.5);
xline(r_1(locs), '--r', 'linewidth', 2)
hold off
grid on
xlabel('$x_1\;[m]$', 'interpreter', 'latex', 'fontsize', 17)
ylabel('$\angle Z_{IN}\;[dB]$', 'interpreter', 'latex', 'fontsize', 17)

%% By inspection

loc = locs(1);

r1 = r_1(loc);
r2 = r_2(loc);
S1 = S_1(loc);
S2 = S_2(loc);
x1 = x_1(loc);
x2 = x_2(loc);
theta1 = theta_1(loc);
theta2 = theta_2(loc);
dL = dL_(loc);
L = L_(loc);
M = M_(loc);

d1 = 2*r1;                  %diameter of the resonator at the head [m]
d2 = 2*r2;                  %diameter of the resonator at the foot [m]


%% Check on the actual input impedance

freq = linspace(1,2000,10000);
omega = 2*pi*freq;
k = omega/c;

theta_1 = (1./k).*atan(k*x1);
theta_2 = (1./k).*atan(k*x2);
M_ = (rho*tan(k*Delta_L))./(k*S2); 

Z_IN = ((1i*rho*c)./S2) .* ((sin(k*L).*sin(k.*theta_1))./sin(k.*(L+theta_1))) + 1i*omega*M;

figure()
subplot(2,1,1)
plot(freq, db(abs(Z_IN)), 'linewidth', 1.5)
grid on
xlabel('$x_1\;[m]$', 'interpreter', 'latex', 'fontsize', 17)
ylabel('$|Z_{IN}|\;[dB]$', 'interpreter', 'latex', 'fontsize', 17)
subplot(2,1,2)
plot(freq, angle(Z_IN), 'linewidth', 1.5)
grid on
xlabel('$x_1\;[m]$', 'interpreter', 'latex', 'fontsize', 17)
ylabel('$\angle Z_{IN}\;[dB]$', 'interpreter', 'latex', 'fontsize', 17)


%% b) Position of the last finger hole
%approssimo il cono ad un cilindro

G4 = 392;               %frequency with last hole opened [Hz]
omega_G4 = 2*pi*G4;     %frequency [rad/s]
k_G4 = omega_G4/c;      %wavenumber [m^-1]
lambda_G4 = 2*pi/k_G4;  %wavelength [m]

%write the input impedance as a function of D
D = linspace(0,L0,1000);        %distance of the hole from the foot
delta = D + dL^2./(D+2*dL);     %difference between equivalent resonators

L_G4 = L - delta;               %length of a the pipe if the hole is opened

%Input impedance as a function of L_G4 -> delta -> D
Z_IN = (1i*rho*c./S2).*(sin(k_G4*L_G4).*sin(k_G4*theta1))./sin(k_G4*(L_G4+theta1)) + 1i*omega_G4*M;

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

Z_IN = (1i*rho*c./S2).*(sin(k_A4*L_A4).*sin(k_A4*theta1))./sin(k_A4*(L_A4+theta1)) + 1i*omega_G4*M;

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
loc = locs(1);
D_2_prime = D(loc);            %distance between holes
D_A4 = D_2_prime + D_G4;        %distance of the hole from the foot
x_A4 = L0 - D_A4;               %distance of the hole from the head
delta_A4 = delta(loc);
dL_2 = D_2_prime + dL_1 - delta_A4;  %end correction of the equivalent resonator


%% Plot the flute
%we already have the diameters at head and foot, let's rotate the flute
r_head = r2; r_foot = r1;
%let's place the input at coord x=0
x_head = 0;
x_foot = L0;
%the last hole is positioned at x_G4 and has section (pi*r_foot^2)
x_lh_right = x_G4+r_foot;
x_lh_left = x_G4-r_foot;
r_lh_right = (x_foot-x_lh_right)*tan(alpha)+r_foot;
r_lh_left = (x_foot-x_lh_left)*tan(alpha)+r_foot;
%the second to last hole is positioned at x_A4 and has section (pi*r2^2)
x_slh_right = x_A4+r_foot;
x_slh_left = x_A4-r_foot;
r_slh_right = (x_foot-x_slh_right)*tan(alpha)+r_foot;
r_slh_left = (x_foot-x_slh_left)*tan(alpha)+r_foot;

figure()
hold on
yline(0, '-.')
%line([x_head,x_foot], [r_head,r_foot], 'color', 'k') %for check
line([x_head,x_slh_left], [r_head,r_slh_left], 'linewidth', 1.5)
line([x_slh_right,x_lh_left], [r_slh_right,r_lh_left], 'linewidth', 1.5)
line([x_lh_right,x_foot], [r_lh_right,r_foot], 'linewidth', 1.5)
line([x_head,x_foot], [-r_head,-r_foot], 'linewidth', 1.5)
xline(x_G4, '-.')
xline(x_A4, '-.')
hold off
grid on
xlim([-0.025,0.475])
ylim([-0.25,0.25])


%% 2) FLUE CHANNEL AND MOUTH

centr = 2000;           %spectrum centroid [Hz]
omega_c = 2*pi*centr;
Delta_P = 62;           %pressure difference between channel and mouth [Pa]

%% d) Channel thickness, jet velocity, Reynolds number

Uj = sqrt(2*Delta_P/rho);       %central velocity of the jet [m/s]
h = 0.3*Uj/centr;               %thickness of the jet
Re = Uj*h/nu;                   %Reynolds number

if Re < 2000
    disp('The jet remains laminar for a short distance after the flue channel exit')
elseif Re > 3000
    disp('The jet becomes turbolent immediately downstream the flue exit')
else
    disp('The regime is transient between laminar and turbolent')
end

%% e) Thickness of the boundary layer at the channel exit

Lchannel = 0.02;            %length of the flue channel [m]

bl = sqrt(nu*Lchannel/Uj);  %boundary layer thickness at channel exit [m]

%% f) Magnitude of oscillations of the jet at the labium

W = 0.004;              %distance between channel exit and labium [m]
SPL = 50;               %sound pressure level at the flue channel exit [dB]

f = F4;                 %frequency impressed on the jet by the resonator
omega = 2*pi*F4;
k = omega/c;
Str = f*h/Uj;

M = rho*0.04/S2;
Z_m = 1i*omega*M;
theta1 = atan(k*x1)/k;
Z_IN = (1i*rho*c/S2) * (sin(k*L)*sin(k*theta1))/sin(k*(L+theta1)) + Z_m;

Delta_P = 20e-6*10^(SPL/20);    %pressure source at the flue channel exit
Vac = Delta_P/abs(Z_IN);        %acoustic field velocity at the labium

alpha_i = 0.4/h;                %obtained from diagram inspection

delta_j = Lchannel/sqrt(Re);
delta_ac = sqrt(2*nu/omega);

eta = abs((Vac*h*delta_j)/(Uj*delta_ac)*exp(alpha_i*W));

