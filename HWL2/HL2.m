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
dl = 8/(3*pi) * a;  %end correction
l_tot = l+2*dl;     %length + virtual elongation of the neck


%% Equivalent circuit parameters

C = V/(rho*c^2);        %condenser for the air volume
L = rho*(l_tot)/S;      %inductor for the mass of the venting tube
R = rho*c/S;            %resistance for the venting tube


%% Analytical admittance computation

omega = linspace(0,10000*2*pi, 100000);
Z = 1i*omega*L + R + 1./(1i*omega*C);
Y = 1./Z;

f_axis = omega./(2*pi);

figure()
subplot(2,1,1)
plot(f_axis, db(abs(Y)), 'linewidth', 2.5)
grid on
xlim([0,3000])
ylim([-140,0])
xlabel('Frequency [Hz]', 'fontsize', 17)
ylabel('$|Y|\,[dB]$', 'interpreter', 'latex', 'fontsize', 17)
title('Magnitude', 'fontsize', 20)
subplot(2,1,2)
plot(f_axis, angle(Y), 'linewidth', 2.5)
grid on
xlim([0,3000])
xlabel('Frequency [Hz]', 'fontsize', 17)
ylabel('$\angle Y\,[deg]$', 'interpreter', 'latex', 'fontsize', 17)
yticks([-pi,-pi/2,0,pi/2,pi])
yticklabels({'-180','-90','0','90','180'})
title('Phase', 'fontsize', 20)

% resonance frequency
f_res_analytical = (c/(2*pi)) * sqrt(S/(V*l_tot));
%f_res_analytical = 1/(2*pi*(L*C)^(1/2));


%% A) Transfer function (admittance) in simscape

U = fft(out.U.Data);
p = fft(out.p.Data);
FRF = U./p;
FRF = FRF(1:(length(FRF)-1)/2);

time_step = out.SimulationMetadata.ModelInfo.SolverInfo.FixedStepSize;
fs = 1/time_step;
f = linspace(0, fs/2, length(FRF));

figure()
subplot(2,1,1)
plot(f, db(abs(FRF)), 'linewidth', 2.5)
grid on
xlim([0,3000])
ylim([-140,0])
xlabel('Frequency [Hz]', 'fontsize', 17)
ylabel('$|Y|\,[dB]$', 'interpreter', 'latex', 'fontsize', 17)
title('Magnitude', 'fontsize', 20)
subplot(2,1,2)
plot(f, angle(FRF), 'linewidth', 2.5)
grid on
xlim([0,3000])
xlabel('Frequency [Hz]', 'fontsize', 17)
ylabel('$\angle Y\,[deg]$', 'interpreter', 'latex', 'fontsize', 17)
yticks([-pi,-pi/2,0,pi/2,pi])
yticklabels({'-180','-90','0','90','180'})
title('Phase', 'fontsize', 20)

% resonance frequency
[pks, locs] = findpeaks(db(abs(FRF(1:3000))));
f_res_simscape = f(locs);



%% Tree admittance computation as in the article

K = 4;
N = 2;

Z = 1i*omega*L + R + 1./(1i*omega*C);
Zden = zeros(1,length(Z));

for k = 1:K-1
    for n = 1:N
        Zden = Zden + 1./Z;
    end
    Z = 1i*omega*L + R + 1./(1i*omega*C + Zden);
    Zden = zeros(1,length(Z));
end

Ya = 1./Z;

figure()
subplot(2,1,1)
plot(f_axis, db(abs(Ya)), 'linewidth', 2.5)
grid on
xlim([0,3000])
xlabel('Frequency [Hz]', 'fontsize', 17)
ylabel('$|Y|\,[dB]$', 'interpreter', 'latex', 'fontsize', 17)
title('Magnitude', 'fontsize', 20)
subplot(2,1,2)
plot(f_axis, angle(Ya), 'linewidth', 2.5)
grid on
xlim([0,3000])
xlabel('Frequency [Hz]', 'fontsize', 17)
ylabel('$\angle Y\,[deg]$', 'interpreter', 'latex', 'fontsize', 17)
yticks([-pi,-pi/2,0,pi/2,pi])
yticklabels({'-180','-90','0','90','180'})
title('Phase', 'fontsize', 20)

[pks,locs] = findpeaks(db(abs(Ya)));
f_res_tree_a = f_axis(locs);


%% Tree admittance computation as in the HW

K = 4;
N = 2;

Z1 = 1i*omega*L + R + 1./(1i*omega*C);
Zden = 1./Z1;
Zprev = zeros(1,length(Z));

for k = 1:K-1
    for n = 1:N-1
        Zden = Zden + 1./Z1;
    end
    Zden = Zden + Zprev;
    Z = 1i*omega*L + R + 1./(1i*omega*C + Zden);
    Zprev = 1./Z;
    Zden = zeros(1,length(Z));
end

Yb = 1./Z;

figure()
subplot(2,1,1)
plot(f_axis, db(abs(Yb)), 'linewidth', 2.5)
grid on
xlim([0,3000])
xlabel('Frequency [Hz]', 'fontsize', 17)
ylabel('$|Y|\,[dB]$', 'interpreter', 'latex', 'fontsize', 17)
title('Magnitude', 'fontsize', 20)
subplot(2,1,2)
plot(f_axis, angle(Yb), 'linewidth', 2.5)
grid on
xlim([0,3000])
xlabel('Frequency [Hz]', 'fontsize', 17)
ylabel('$\angle Y\,[deg]$', 'interpreter', 'latex', 'fontsize', 17)
yticks([-pi,-pi/2,0,pi/2,pi])
yticklabels({'-180','-90','0','90','180'})
title('Phase', 'fontsize', 20)

[pks,locs] = findpeaks(db(abs(Yb)));
f_res_tree_b = f_axis(locs);


%% comparison

figure()
subplot(2,1,1)
hold on
plot(f, db(abs(FRF)), 'linewidth', 2)
plot(f_axis, db(abs(Ya)), '--', 'linewidth', 2)
plot(f_axis, db(abs(Yb)), '--', 'linewidth', 2)
hold off
grid on
legend('simscape', 'article', 'hw')
xlim([0,3000])
xlabel('Frequency [Hz]', 'fontsize', 17)
ylabel('$|Y|\,[dB]$', 'interpreter', 'latex', 'fontsize', 17)
title('Magnitude', 'fontsize', 20)
subplot(2,1,2)
hold on
plot(f, angle(FRF), 'linewidth', 2)
plot(f_axis, angle(Ya), '--', 'linewidth', 2)
plot(f_axis, angle(Yb), '--', 'linewidth', 2)
hold off
grid on
xlim([0,3000])
xlabel('Frequency [Hz]', 'fontsize', 17)
ylabel('$\angle Y\,[deg]$', 'interpreter', 'latex', 'fontsize', 17)
yticks([-pi,-pi/2,0,pi/2,pi])
yticklabels({'-180','-90','0','90','180'})
title('Phase', 'fontsize', 20)