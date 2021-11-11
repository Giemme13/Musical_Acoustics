rho=1.225; %air density
c=434; %speed in air
Zair=0; %air impedance %per ora 0
a0=0.01;
m=4;
L=0.4;
fmax=2000;

f=1:10:fmax; %frequency range
omega=2*pi*f;
k=omega/c;
b=sqrt(k.^2-m^2);
theta=atan(m./b);

%cross sections
S1=pi*a0^2;
S2=pi*(a0*exp(L*m))^2;

%exponential horn
numEXP=Zair*cos(b.*L+theta)+1i*(rho*c/S2)*sin(b.*L);
denEXP=1i*Zair*sin(b.*L)+rho*c./(S2*cos(b.*L-theta));
ZinEXP=(rho*c/S1)*numEXP./denEXP;
%plot(f, db(abs(ZinEXP)))

%1 conical horn
%x2=(exp(m*L)*L)/(exp(m*L)-1);
%x1=x2-L;
%theta1=atan(k*x1);
%theta2=atan(k*x2);
%numCON=1i*Zair*(sin(k*L-theta2)/sin(theta2))+(rho*c/S2)*sin(k*L);
%denCON=Zair*(sin(k*L+theta1-theta2)./(sin(theta1).*sin(theta2)))-(1i*rho*c/S2).*(sin(k*L+theta1)./sin(theta1));
%ZinCON=numCON./denCON;
%plot(f,db(abs(ZinCON)))

%recursive attempt
subdivision=1:15;
delta=L./subdivision;
ZL=0;
Zin=0;

%empty arrays for errors
e1=zeros(1,length(subdivision));
e2=zeros(1,length(subdivision));

for i=1:length(subdivision)
    for j=1:i
        if j==1
            ZL=Zair;
        end
        index=i+1-j;
        %parameters
        x2=(exp(m*delta(index)))/(exp(m*delta(index))-1);
        x1=x2-delta(index);
        S1=pi*(a0*exp(x1*m))^2;
        S2=pi*(a0*exp(x2*m))^2;
        theta1=atan(k*x1);
        theta2=atan(k*x2);
        %impedance
        Zin=(rho*c./S1).*(1i*ZL.*(sin(k.*delta(index)-theta2)./sin(theta2))+(rho*c/S2)*sin(k.*delta(index)))./(ZL.*(sin(k*delta(index)+theta1-theta2)./(sin(theta1).*sin(theta2)))-(1i*rho*c/S2).*(sin(k*delta(index)+theta1)./sin(theta1)));
        ZL=Zin;
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

plot(subdivision, e1)