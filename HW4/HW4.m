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

S = pi*D^2/4;               %cross sectional area of the venting tube [m^2]
delta_L = 0.85*D/2;         %end correction on one end of the venting pipe [m]

%% Equivalent electrical circuit

% diaphragm
L1 = m;                     %inductor for the mass of the loudspeaker [kg]
k1 = k;
C1 = 1/k1;                  %condenser for the equivalent spring of the loudspeaker [s^2/kg]
R1 = R;                     %resistance for the damping of the loudspeaker [kg/s]

% cavity
C2 = V/(rho*c^2*S^2);       %condenser for the air volume [s^2/kg]
k2 = 1/C2;

%venting tube
L2 = rho*(L+2*delta_L)*S;   %inductor for the mass of the venting tube [kg]
R2 = rho*c*S;               %resistance for the dissipation inside the venting tube [kg/s]


%% Admittance computation

omega = linspace(0,10000*2*pi, 100000);

%diaphragm
Zmd = 1i.*omega.*L1;
Zrd = R1;
Zcd = 1./(1i.*omega*C1);

Zd = Zmd + Zrd + Zcd;

%cavity
Zcc = 1./(1i.*omega*C2);

%venting tube
Zmt = 1i.*omega*L2;
Zrt = R2;

Zt = Zmt + Zrt;

Z = Zd + (1./(1./Zcc + 1./Zt));
Y = 1./Z;
f = omega./(2*pi);

figure(1)
subplot(2,1,1)
plot(f, db(abs(Y)), 'linewidth', 2.5)
grid on
xlabel('Frequency [Hz]', 'fontsize', 17)
ylabel('$|Y|\,[dB]$', 'interpreter', 'latex', 'fontsize', 17)
title('Magnitude', 'fontsize', 20)
subplot(2,1,2)
plot(f, angle(Y), 'linewidth', 2.5)
grid on
xlabel('Frequency [Hz]', 'fontsize', 17)
ylabel('$\angle Z\,[deg]$', 'interpreter', 'latex', 'fontsize', 17)
yticks([-pi,-pi/2,0,pi/2,pi])
yticklabels({'-180','-90','0','90','180'})
title('Phase', 'fontsize', 20)


%% FRF of mechanical circuit

vel = fft(out.out_velocity.Data);
force = fft(out.in_force.Data);
FRF_mech = vel./force;
FRF_mech = FRF_mech(1:(length(FRF_mech)-1)/2);

time_step = out.SimulationMetadata.ModelInfo.SolverInfo.FixedStepSize;
fs = 1/time_step;
f_mech = linspace(0, fs/2, length(FRF_mech));

figure(2)
subplot(2,1,1)
plot(f_mech, db(abs(FRF_mech)), 'linewidth', 2.5)
grid on
xlabel('Frequency [Hz]', 'fontsize', 17)
ylabel('$|Z|\,[dB]$', 'interpreter', 'latex', 'fontsize', 17)
title('Magnitude', 'fontsize', 20)
subplot(2,1,2)
plot(f_mech, angle(FRF_mech), 'linewidth', 2.5)
grid on
xlabel('Frequency [Hz]', 'fontsize', 17)
ylabel('$\angle Z\,[deg]$', 'interpreter', 'latex', 'fontsize', 17)
yticks([-pi,-pi/2,0,pi/2,pi])
yticklabels({'-180','-90','0','90','180'})
title('Phase', 'fontsize', 20)

%% FRF of electrical circuit

current = fft(out.out_current.Data);
voltage = fft(out.in_voltage.Data);
FRF_elec = current./voltage;
FRF_elec = FRF_elec(1:(length(FRF_elec)-1)/2);

time_step = out.SimulationMetadata.ModelInfo.SolverInfo.FixedStepSize;
fs = 1/time_step;
f_elec = linspace(0, fs/2, length(FRF_elec));

figure(3)
subplot(2,1,1)
plot(f_elec, db(abs(FRF_elec)), 'linewidth', 2.5)
grid on
xlabel('Frequency [Hz]', 'fontsize', 17)
ylabel('$|Z|\,[dB]$', 'interpreter', 'latex', 'fontsize', 17)
title('Magnitude', 'fontsize', 20)
subplot(2,1,2)
plot(f_elec, angle(FRF_elec), 'linewidth', 2.5)
grid on
xlabel('Frequency [Hz]', 'fontsize', 17)
ylabel('$\angle Z\,[deg]$', 'interpreter', 'latex', 'fontsize', 17)
yticks([-pi,-pi/2,0,pi/2,pi])
yticklabels({'-180','-90','0','90','180'})
title('Phase', 'fontsize', 20)


%% Comparison

figure(4)
subplot(2,1,1)
hold on
plot(f, db(abs(Y)), 'linewidth', 2.5)
plot(f_mech, db(abs(FRF_mech)), '-.', 'linewidth', 2.5)
plot(f_elec, db(abs(FRF_elec)), '--', 'linewidth', 2.5)
hold off
grid on
xlim([0,500])
legend('Analytical solution', 'Mechanical circuit', 'Electrical circuit', 'fontsize', 20)
xlabel('Frequency [Hz]', 'fontsize', 17)
ylabel('$|Y|\,[dB]$', 'interpreter', 'latex', 'fontsize', 17)
title('Magnitude', 'fontsize', 20)
subplot(2,1,2)
hold on
plot(f, angle(Y), 'linewidth', 2.5)
plot(f_mech, angle(FRF_mech), '-.', 'linewidth', 2.5)
plot(f_elec, angle(FRF_elec), '--', 'linewidth', 2.5)
hold off
grid on
xlim([0,500])
xlabel('Frequency [Hz]', 'fontsize', 17)
ylabel('$\angle Y\,[deg]$', 'interpreter', 'latex', 'fontsize', 17)
yticks([-pi,-pi/2,0,pi/2,pi])
yticklabels({'-180','-90','0','90','180'})
title('Phase', 'fontsize', 20)