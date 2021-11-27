clear all 
close all
clc

%% Data

V = 0.1;    %air volume [m^3]
l = 0.1;    %length of the neck [m]
S = 100;    %section of the neck [m^2]
c = 343;    %speed of sound in air [m/s]
rho = 1.2;  %density of the air [kg/m^3]


%% virtual neck

a = sqrt(S/pi);     %radius of the neck
dl = 0.85 * a;      %end correction
l_tot = l+2*dl;     %length + virtual elongation of the neck


%% Equivalent circuit parameters

C = V/(rho*c^2);        %condenser for the air volume
L = rho*(l_tot)/S;      %inductor for the mass of the venting tube
R = rho*c/S;            %resistance for the venting tube


%% A) Transfer function (admittance) in simscape

U = fft(out.U.Data);
p = fft(out.p.Data);
FRF = U./p;
FRF = FRF(1:(length(FRF)-1)/2);

time_step = out.SimulationMetadata.ModelInfo.SolverInfo.FixedStepSize;
fs = 1/time_step;
f = linspace(0, fs/2, length(FRF));

figure(1)
subplot(2,1,1)
plot(f, db(abs(FRF)), 'linewidth', 2.5)
grid on
xlabel('Frequency [Hz]', 'fontsize', 17)
ylabel('$|Y|\,[dB]$', 'interpreter', 'latex', 'fontsize', 17)
title('Magnitude', 'fontsize', 20)
subplot(2,1,2)
plot(f, angle(FRF), 'linewidth', 2.5)
grid on
xlabel('Frequency [Hz]', 'fontsize', 17)
ylabel('$\angle Y\,[deg]$', 'interpreter', 'latex', 'fontsize', 17)
yticks([-pi,-pi/2,0,pi/2,pi])
yticklabels({'-180','-90','0','90','180'})
title('Phase', 'fontsize', 20)


%% B)  Compute the natural frequency

%from simscape
[pks, locs] = findpeaks(db(abs(FRF)));
f_0s = f(locs(1));

%analytically
f_0a = (c/(2*pi)) * sqrt(S/(V*l_tot));

