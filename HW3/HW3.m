%%
clear all
close all
clc

%% DATA

rho = 1.225;     %air density    [kg/m^3]
c = 343;         %speed in air   [m/s]
Zair = 0;        %air impedance %per ora 0
a0 = 0.01;       %radius of exponential horn's section: a = a0*exp(m*x) [m]
m = 4;
L = 0.4;         %total length of the horn [m]
fmin = 1;        %frequency study range [fmin, fmax] [Hz]
fmax = 2000;
freq_res = 10000;  %resolution of the frequency interval

f = linspace(fmin,fmax,freq_res);  %frequency axis [Hz]
omega = 2*pi*f;         %frequency axis [rad/s]
k = omega/c;            %wavenumber [m^-1]
b = sqrt(k.^2-m^2);     %geometrical quantity [m^-1]
theta = atan(m./b);     %geometrical quantity [/]

%cross sections
S1 = pi*a0^2;               %throat section [m^2]
S2 = pi*(a0*exp(L*m))^2;    %mouse section [m^2]

res = 1000;       %number of points used in the plots


%% EXPONENTIAL HORN

x = linspace(0,L,res);
a = a0*exp(m*x);
numEXP = Zair*cos(b*L+theta)+1i*(rho*c/S2)*sin(b*L);
denEXP = 1i*Zair*sin(b*L)+(rho*c/S2)*cos(b*L-theta);
ZinEXP = (rho*c/S1)*numEXP./denEXP;

%figure(1)
%plot(x, a, 'k', 'linewidth', 2)
%hold on
%plot(x, -a, 'k', 'linewidth', 2)
%hold off
%yline(0, '-.')

%figure(2)
%subplot(2,1,1)
%plot(f, db(abs(ZinEXP)))
%subplot(2,1,2)
%plot(f, angle(ZinEXP))


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


%figure(3)
%plot(x, a, 'k', 'linewidth', 2)
%hold on
%plot(x, -a, 'k', 'linewidth', 2)
%hold off
%yline(0, '-.')

%figure(4)
%subplot(2,1,1)
%plot(f, db(abs(ZinCON_1)))
%subplot(2,1,2)
%plot(f, angle(ZinCON_1))


%% n CONICAL HORNS

n = 5;
delta = L/n;

points = round(res/n);

Zload = Zair;

for i = 1:n
    a2 = a0*exp(m*(n-i+1)*delta);
    a1 = a0*exp(m*(n-i)*delta);
    x1 = (a1*delta)/(a2-a1);
    x2 = x1+delta;
    theta1 = atan(k*x1);
    theta2 = atan(k*x2);
    S1 = pi*(a1^2);
    S2 = pi*(a2^2);
    disp([a1,a2])
    x = linspace(0, delta, points);
    slope = (a2-a1)/delta;
    a = slope*x + a1;
    
    numCON = 1i*Zload.*(sin(k*delta-theta2)./sin(theta2))+(rho*c/S2)*sin(k*delta);
    denCON = Zload.*(sin(k*delta+theta1-theta2)./(sin(theta1).*sin(theta2)))-(1i*rho*c/S2).*(sin(k*delta+theta1)./sin(theta1));
    ZinCON_n = (rho*c/S1)*numCON./denCON;
    
    Zload = ZinCON_n;
    
    %figure(5)
    %plot((n-i)*delta+x, a, 'k', 'linewidth', 2)
    %hold on
    %plot((n-i)*delta+x, -a, 'k', 'linewidth', 2)
    %yline(0, '-.')
end

%hold off

%figure(6)
%subplot(2,1,1)
%plot(f,db(abs(ZinCON_n)))
%subplot(2,1,2)
%plot(f, angle(ZinCON_n))


%% FROM 1 TO n CONICAL HORNS

n = 10;

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

        %figure(4+j)
        %plot((j-i)*delta+x, a, 'k', 'linewidth', 2)
        %hold on
        %plot((j-i)*delta+x, -a, 'k', 'linewidth', 2)
        %yline(0, '-.')
    end
    %hold off
    
    %figure(4+n+j)
    %subplot(2,1,1)
    %plot(f,db(abs(ZinCON_n)))
    %subplot(2,1,2)
    %plot(f, angle(ZinCON_n))
    
    for i = 1:length(f)
        e1(j) = e1(j) + (abs(ZinCON_n(i) - ZinEXP(i)))^2;
    end
    e1(j) = e1(j)/(2*pi*(fmax-fmin));
    
    
    [pks_CON, locs_CON] = findpeaks(abs(ZinCON_n));
    [pks_EXP, locs_EXP] = findpeaks(abs(ZinEXP));
    for i = 1:5
        e2(j) = e2(j) + abs(locs_CON(i)-locs_EXP(i));
    end
end

%figure(n+j+4+1)
%stem(deltas,e1)
%figure(n+j+4+2)
%stem(deltas,e2)


%%

%figure(50)
%plot(f,db(abs(ZinEXP)), 'linewidth', 2)
%hold on
%plot(f,db(abs(ZinCON_1)))
%plot(f,db(abs(ZinCON_n)), '-.')
%legend()
%hold off

%% Air load
%first approximation whit e2=0 is for L/10
n=10;
%calculating air impedance load
ZL0=0.25*(rho.*(omega.^2)./(pi*c))+1i*0.61*rho.*omega./(pi*a0*exp(m*L));
flare=atan(a0*(exp(m*L)-1)/L);
ZL=ZL0*2/(1+cos(flare));

%calculating new Zin
delta=1/n;
for i=1:n
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
end