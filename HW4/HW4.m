%%
clear all
close all
clc

%% DATA

rho = 1.225;    %air density [kg/m^3]
c = 344;        %sound velocity in air [m/s]

D = 0.07;       %diameter of the venting tube [m]
L = 0.247;      %length of the venting tube [m]
V = 0.02;       %volume of the cabinet [m^3]
m = 0.0455;     %mass of the diaphragm [kg]
fd = 49;        %damped resonance frequency not considering the cabinet [Hz]
Q = 50;         %quality factor at resonance [\]

%% Mechanical characterization

R = 2*pi*fd*m/Q;     %equivalent damping coefficient of the mass-spring-damper system
%omega_d = sqrt(k/m - R^2/(4m^2))
k = m*((2*pi*fd)^2 + R^2/(4*m^2));   %equivalent stiffness coefficient of
                                     %mass-spring-damped system
f0 = sqrt(k/m)/(2*pi);      %natural resonance frequency [Hz]

%% Equivalent circuit

S = pi*D^2/4;               %cross sectional area of the venting tube [m^2]
delta_L = 0.85*D/2;         %end correction on one end of the venting pipe

L2 = rho*(L+2*delta_L)/S;   %inductor for the mass of the venting tube
R2 = rho*c/S;               %resistance for the dissipation inside the venting tube
C2 = V/(rho*c^2);           %condenser for the air volume
L1 = m;                     %inductor for the mass of the loudspeaker
C1 = 1/k;                   %condenser for the equivalent spring of the loudspeaker
R1 = R;                     %resistance for the damping of the loudspeaker


%% FRF of mechanical circuit

num_mech = fft(out.out_velocity.Data);
den_mech = fft(out.in_force.Data);
FRF_mech = num_mech./den_mech;

f_mech = linspace(1,1000,length(out.out_velocity.Data));

figure()
subplot(2,1,1)
plot(f_mech, db(abs(FRF_mech)))
subplot(2,1,2)
plot(f_mech, angle(FRF_mech))

%% FRF of electrical circuit

num_elec = fft(out.out_current.Data);
den_elec = fft(out.in_voltage.Data);
FRF_elec = num_elec./den_elec;

f_elec = linspace(1,1000,length(out.out_current.Data));

figure()
subplot(2,1,1)
plot(f_elec, db(abs(FRF_elec)))
subplot(2,1,2)
plot(f_elec, angle(FRF_elec))


%% Comparison

figure()
subplot(2,1,1)
hold on
plot(f_mech, db(abs(FRF_mech)))
plot(f_elec, db(abs(FRF_elec)))
hold off
subplot(2,1,2)
hold on
plot(f_mech, angle(FRF_mech))
plot(f_elec, angle(FRF_elec))
hold off