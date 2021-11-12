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
fmin = 0;        %frequency study range [fmin, fmax] [Hz]
fmax = 2000;

f = linspace(fmin,fmax,10000);  %frequency axis [Hz]
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
    
    figure(5)
    plot((n-i)*delta+x, a, 'k', 'linewidth', 2)
    hold on
    plot((n-i)*delta+x, -a, 'k', 'linewidth', 2)
    yline(0, '-.')
end

hold off

figure(6)
subplot(2,1,1)
plot(f,db(abs(ZinCON_n)))
subplot(2,1,2)
plot(f, angle(ZinCON_n))


%%

figure(10)
plot(f,db(abs(ZinEXP)), 'linewidth', 2)
hold on
plot(f,db(abs(ZinCON_1)))
plot(f,db(abs(ZinCON_n)))
legend()
hold off

%% n CONICAL HORNS

%recursive attempt
n = 1:10;
delta = L./n;
ZL = 0;
Zin = 0;

%empty arrays for errors
e1=zeros(1,length(n));
e2=zeros(1,length(n));

for i=1:length(n)
    for j=1:i
        %parameters
        %index=i+1-j;
        if j==1
            ZL=Zair;
            x2=(exp(m*L))*L/(exp(m*L)-1);
            x1=x2-delta(i);
        else
            %x2=(exp(m*(L-(j-1)*delta(i))))*(L-(j-1)*delta(i))/(exp(m*(L-(j-1)*delta(i)))-1);
            x2=x1;
            x1=x2-(j-1)*delta(i);
        end
        S1=pi*(max([a0, a0*exp(x1*m)]))^2;
        S2=pi*(max([a0, a0*exp(x2*m)]))^2;
        theta1=atan(k*x1);
        theta2=atan(k*x2);
        %impedance
        Zin=(rho*c./S1).*(1i*ZL.*(sin(k.*delta(i)-theta2)./sin(theta2))+(rho*c/S2)*sin(k.*delta(i)))./(ZL.*(sin(k*delta(i)+theta1-theta2)./(sin(theta1).*sin(theta2)))-(1i*rho*c/S2).*(sin(k*delta(i)+theta1)./sin(theta1)));
        ZL=Zin;
        %fprintf('%f \n', S2,S1)
    end
    e1Temp=0;
    e2Temp=0;
    for h=1:(length(f)-1)
       e1Temp=e1Temp+(abs(Zin(h)-ZinEXP(h)))^2;
    end    
    %calculate errors e1 and e2 for given number of subdivisions
    e1(i)=e1Temp/(fmax-1);
    %e2=...;
end

%plot(n, e1)