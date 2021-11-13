%%
clear all
close all
clc

%% DATA

rho = 1.225;     %air density    [kg/m^3]
c = 343;         %speed in air   [m/s]
Zair = rho*c;    %air impedance
a0 = 0.01;       %radius of exponential horn's section: a = a0*exp(m*x) [m]
m = 4;
L = 0.4;         %total length of the horn [m]
fmin = 1;        %frequency study range [fmin, fmax] [Hz]
fmax = 2000;
freq_res = 10000;  %resolution of the frequency interval for the impedances' plots

f = linspace(fmin,fmax,freq_res);  %frequency axis [Hz]
omega = 2*pi*f;         %frequency axis [rad/s]
k = omega/c;            %wavenumber [m^-1]
b = sqrt(k.^2-m^2);     %geometrical quantity [m^-1]
theta = atan(m./b);     %geometrical quantity [/]

%cross sections
S1 = pi*a0^2;               %throat section [m^2]
S2 = pi*(a0*exp(L*m))^2;    %mouth section [m^2]

res = 1000;       %number of points used in the axial sections' plots


%% EXPONENTIAL HORN

x = linspace(0,L,res);
a = a0*exp(m*x);
numEXP = Zair*cos(b*L+theta)+1i*(rho*c/S2)*sin(b*L);
denEXP = 1i*Zair*sin(b*L)+(rho*c/S2)*cos(b*L-theta);
ZinEXP = (rho*c/S1)*numEXP./denEXP;

figure(1)
plot(x, a, 'k', 'linewidth', 2)
hold on
plot(x, -a, 'k', 'linewidth', 2)
hold off
yline(0, '-.')

figure(2)
subplot(2,1,1)
plot(f, db(abs(ZinEXP)))
subplot(2,1,2)
plot(f, angle(ZinEXP))


%% 1 CONICAL HORN
a1 = a0*exp(m*0);
a2 = a0*exp(m*L);
% x2 : a2 = L : (a2-a1)
x2 = (a2*L)/(a2-a1);
x1 = x2-L;
theta1 = atan(k*x1);
theta2 = atan(k*x2);

x = linspace(0,L,res);
slope = (a2-a1)/L;
a = slope*x + a1;
Zload = Zair;
numCON = 1i*Zload*(sin(k*L-theta2)/sin(theta2))+(rho*c/S2)*sin(k*L);
denCON = Zload*(sin(k*L+theta1-theta2)./(sin(theta1).*sin(theta2)))-(1i*rho*c/S2).*(sin(k*L+theta1)./sin(theta1));
ZinCON_1 = (rho*c/S1)*numCON./denCON;


figure(3)
plot(x, a, 'k', 'linewidth', 2)
hold on
plot(x, -a, 'k', 'linewidth', 2)
hold off
yline(0, '-.')

figure(4)
subplot(2,1,1)
plot(f, db(abs(ZinCON_1)))
subplot(2,1,2)
plot(f, angle(ZinCON_1))


%% FROM 1 TO n CONICAL HORNS

n = 20;

deltas = zeros(1,n);

e1 = zeros(1, n);
e2 = zeros(1, n);

for j = 1:n
    delta = L/j;
    deltas(j) = delta;
    points = round(res/j);
    
    Zload = Zair;
    for i = 1:j
        a2 = a0*exp(m*(j-i+1)*delta);
        a1 = a0*exp(m*(j-i)*delta);
        x1 = (a1*delta)/(a2-a1);
        x2 = x1+delta;
        theta1 = atan(k*x1);
        theta2 = atan(k*x2);
        S1 = pi*(a1^2);
        S2 = pi*(a2^2);
        x = linspace(0, delta, points);
        slope = (a2-a1)/delta;
        a = slope*x + a1;

        numCON = 1i*Zload.*(sin(k*delta-theta2)./sin(theta2))+(rho*c/S2)*sin(k*delta);
        denCON = Zload.*(sin(k*delta+theta1-theta2)./(sin(theta1).*sin(theta2)))-(1i*rho*c/S2).*(sin(k*delta+theta1)./sin(theta1));
        ZinCON_n = (rho*c/S1)*numCON./denCON;

        Zload = ZinCON_n;
        
        if j == n
            figure(7)
            plot((j-i)*delta+x, a, 'k', 'linewidth', 2)
            hold on
            plot((j-i)*delta+x, -a, 'k', 'linewidth', 2)
            yline(0, '-.')
            xline((j-i)*delta, '--', 'color', [0.8500, 0.3250, 0.0980], 'linewidth', 0.2)
        end
    end
    hold off
    
    if j == n
        figure(8)
        subplot(2,1,1)
        plot(f,db(abs(ZinCON_n)))
        subplot(2,1,2)
        plot(f, angle(ZinCON_n))
    end
    
    for i = 1:length(f)
        e1(j) = e1(j) + (abs(ZinCON_n(i) - ZinEXP(i)))^2;
    end
    e1(j) = e1(j)/(2*pi*(fmax-fmin));
    
    
    [pks_CON, locs_CON] = findpeaks(real(ZinCON_n));
    [pks_EXP, locs_EXP] = findpeaks(real(ZinEXP));
    
    argmax_omega_CON = f(locs_CON);
    argmax_omega_EXP = f(locs_EXP);
    
    for i = 1:5
        e2(j) = e2(j) + abs(argmax_omega_CON(i)-argmax_omega_EXP(i));
    end
end

figure(9)
stem(deltas,e1)
figure(10)
stem(deltas,e2)



%% DIFFERENT SAMPLING STRATEGY

%uniform sampling on y-axis
n = n;

%geometry
a1 = a0*exp(m*0);
a2 = a0*exp(m*L);
H = a2-a1;

e1_ysamp = zeros(1, n);
e2_ysamp = zeros(1, n);

for j = 1:n
    
    height = H/j;
    
    Zload = Zair;
    
    len = L;
    
    for i = 1:j
        a1 = (j-i)*height + a0;             %y_value of the throat of the (j-i+1)-th cone
        a2 = (j-i+1)*height + a0;           %y_value of the mouth of the (j-i+1)-th cone
        d1 = (1/m)*log(a1/a0);              %x_value of the throat of the (j-i+1)-th cone
        d2 = (1/m)*log(a2/a0);              %x_value of the mouth of the (j-i+1)-th cone
        delta = d2-d1;                      %distance between throat and mouth
        x1 = (a1*delta)/(a2-a1);
        x2 = x1+delta;
        theta1 = atan(k*x1);
        theta2 = atan(k*x2);
        S1 = pi*(a1^2);
        S2 = pi*(a2^2);
        
        ratio = L/delta;                    %number of points on x_axis for the plot of the axial section
        points = round(res/ratio);
        x = linspace(0, delta, points);
        slope = (a2-a1)/delta;
        a = slope*x + a1;

        numCON = 1i*Zload.*(sin(k*delta-theta2)./sin(theta2))+(rho*c/S2)*sin(k*delta);
        denCON = Zload.*(sin(k*delta+theta1-theta2)./(sin(theta1).*sin(theta2)))-(1i*rho*c/S2).*(sin(k*delta+theta1)./sin(theta1));
        ZinCON_ysamp = (rho*c/S1)*numCON./denCON;
        
        Zload = ZinCON_ysamp;
        
        if j == n
            figure(11)
            plot((len-delta)+x, a, 'k', 'linewidth', 2)
            hold on
            plot((len-delta)+x, -a, 'k', 'linewidth', 2)
            yline(0, '-.')
            xline(len-delta, '--', 'color', [0.8500, 0.3250, 0.0980], 'linewidth', 0.2)
        end
        len = len-delta;        
    end
    hold off
    
    if j == n
        figure(12)
        subplot(2,1,1)
        plot(f, db(abs(ZinCON_ysamp)))
        subplot(2,1,2)
        plot(f, angle(ZinCON_ysamp))
    end

    
    for i = 1:length(f)
        e1_ysamp(j) = e1_ysamp(j) + (abs(ZinCON_ysamp(i) - ZinEXP(i)))^2;
    end
    e1_ysamp(j) = e1_ysamp(j)/(2*pi*(fmax-fmin));
    
    
    [pks_CON, locs_CON] = findpeaks(real(ZinCON_ysamp));
    [pks_EXP, locs_EXP] = findpeaks(real(ZinEXP));
    
    argmax_omega_CON = f(locs_CON);
    argmax_omega_EXP = f(locs_EXP);
    
    for i = 1:5
        e2_ysamp(j) = e2_ysamp(j) + abs(argmax_omega_CON(i)-argmax_omega_EXP(i));
    end
end

figure(13)
stem(linspace(1,n,n),e1_ysamp)
figure(14)
stem(linspace(1,n,n),e2_ysamp)


%% COMPARISON

figure(15)
plot(f,db(abs(ZinEXP)), 'linewidth', 2)
hold on
plot(f,db(abs(ZinCON_1)))
plot(f,db(abs(ZinCON_n)), '-.')
plot(f,db(abs(ZinCON_ysamp)), '-.')
legend()
hold off

%% Air load
%first approximation whit e2=0 is for L/10
n = 10;
delta = L/n;

%calculating air impedance load
a = a0*exp(m*L);        %radius of the unflanged cylinder
Sp = pi*a^2;            %cross section of the unflanged cylinder
ZL0 = 0.25*(rho.*(omega.^2)./(pi*c))+1i*0.61*rho.*omega./(pi*a);

%considering only the 10th cone, before the cylinder
a2 = a0*exp(m*(10*delta));      %end radius
a1 = a0*exp(m*(9*delta));       %start radius
flare = atan((a2-a1)/delta);    %flare angle of the 10th cone
Ss = 2*Sp/(1+cos(flare));       %area of the wavefront leaving te cone

ZL = ZL0*(Sp/Ss);

%calculating new Zin
for i = 1:n
    a2 = a0*exp(m*(n-i+1)*delta);
    a1 = a0*exp(m*(n-i)*delta);
    x1 = (a1*delta)/(a2-a1);
    x2 = x1+delta;
    theta1 = atan(k*x1);
    theta2 = atan(k*x2);
    S1 = pi*(a1^2);
    S2 = pi*(a2^2);

    numCON = 1i*ZL.*(sin(k*delta-theta2)./sin(theta2))+(rho*c/S2)*sin(k*delta);
    denCON = ZL.*(sin(k*delta+theta1-theta2)./(sin(theta1).*sin(theta2)))-(1i*rho*c/S2).*(sin(k*delta+theta1)./sin(theta1));
    ZinCON_10 = (rho*c/S1)*numCON./denCON;
    
    ZL = ZinCON_10;
end


figure(16)
subplot(2,1,1)
plot(f, db(abs(ZinCON_10)))
subplot(2,1,2)
plot(f, angle(ZinCON_10))

%% COMPOUND HORN

%Data
L_cyl = 0.6;
a_cyl = a0;
S_cyl = pi*(a_cyl^2);
Z0 = rho*c/S_cyl;

numCYL = ZinCON_10.*cos(k*L_cyl)+1i*Z0*sin(k*L_cyl);
denCYL = 1i*ZinCON_10.*sin(k*L_cyl)+Z0*cos(k*L_cyl);
ZinCYL = Z0*(numCYL./denCYL);

[pks_COMP, locs_COMP] = findpeaks(abs(ZinCYL));
argmax_omega_COMP = f(locs_COMP);
tenMax = zeros(1,10);

for i=1:10
    tenMax(i) = argmax_omega_COMP(i);
end

figure(17)
subplot(2,1,1)
plot(f, db(abs(ZinCYL)))
subplot(2,1,2)
plot(f, angle(ZinCYL))