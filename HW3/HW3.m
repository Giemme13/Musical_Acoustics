rho=1.225; %air density
c=434; %speed in air
ZL=0; %air impedance %per ora 0
a0=0.01;
m=4;
L=0.4;

f=0:10:2000; %frequency range
omega=2*pi*f;
k=omega/c;
b=sqrt(k.^2-m^2);
theta=atan(m./b);

%cross sections
S1=pi*a0^2;
S2=pi*(a0*exp(L*m))^2;

%exponential horn
numEXP=ZL*cos(b.*L+theta)+1i*(rho*c/S2)*sin(b.*L);
denEXP=1i*ZL*sin(b.*L)+rho*c./(S2*cos(b.*L-theta));
ZinEXP=numEXP./denEXP;
plot(f, db(abs(ZinEXP)))

%1 conical horn
x2=(exp(m*L)*L)/(exp(m*L)-1);
x1=x2-L;
theta1=atan(k*x1);
theta2=atan(k*x2);

numCON=1i*ZL*(sin(k*L-theta2)/sin(theta2))+(rho*c/S2)*sin(k*L);
denCON=ZL*(sin(k*L+theta1-theta2)./(sin(theta1).*sin(theta2)))-(1i*rho*c/S1).*(sin(k*L+theta1)./sin(theta1));
ZinCON=numCON./denCON;
plot(f,db(abs(ZinCON)))


