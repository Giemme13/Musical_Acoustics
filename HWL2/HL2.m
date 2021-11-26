clear all 
close all
clc

%% Data

V=0.1; %air volume
l=0.1;  %length of the neck
S=100;  %section of the neck
c=343;  
rho=1.2;

a=sqrt(S/pi); %radius
l_tot=l+2*0.85*a; %virtual elongation

f_an=(c/(2*pi))*sqrt(S/(V*l_tot)); %natural frequency of the resonator

%% Equivalent circuit parameters

C = V/(rho*c^2);       %condenser for the air volume 

L = rho*(l_tot)/S;   %inductor for the mass of the venting tube 
R = rho*c/S;               %resistance for the venting tube 

%% Transfer function

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
ylabel('$|Z|\,[dB]$', 'interpreter', 'latex', 'fontsize', 17)
title('Magnitude', 'fontsize', 20)
subplot(2,1,2)
plot(f, angle(FRF), 'linewidth', 2.5)
grid on
xlabel('Frequency [Hz]', 'fontsize', 17)
ylabel('$\angle Z\,[deg]$', 'interpreter', 'latex', 'fontsize', 17)
yticks([-pi,-pi/2,0,pi/2,pi])
yticklabels({'-180','-90','0','90','180'})
title('Phase', 'fontsize', 20)
