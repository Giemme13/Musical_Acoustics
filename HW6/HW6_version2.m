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

x_foot = ((pi*c/omega_F4)-L0)/(0.6*tan(alpha)+1);   %closest to vertex
x_head = x_foot+L0;

r_foot = x_foot*tan(alpha);
d_foot = 2*r_foot;
r_head = r_foot+L0*tan(alpha);
d_head = 2*r_head;

dL = r_foot *0.6;
L = L0 + dL;

%% Check on impedance

f = linspace(1,2000,10000);
omega = 2*pi*f;
k = omega/c;

S1 = pi*r_foot^2;
theta1 = atan(k*x_foot)./k;

Z_IN = 1i*rho*c/S1 *(sin(k*L).*sin(k.*theta1))./sin(k.*(L+theta1));

%%
plot(f,db(abs(Z_IN)))

