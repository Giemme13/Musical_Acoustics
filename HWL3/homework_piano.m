%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modeling of musical instruments homework.                               %
% Numerical simulation of piano strings                                   %
% Physical model for a struck string using finite difference.             %
%                                                                         %
% Musical Acoustics course                                                %
% Mirco Pezzoli                                                           %
% 2021                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

%% Setup
% define the string parameters and the simulation variables defined
% according to the provided values and the numerical implementation.
% We want to implement the finite difference scheme of a piano string tuned
% as C2.

% Temporal sampling parameters
duration = 8;       %signal is 8 second long [s]
fs = 4*44100;         %sampling frequency [Hz]
t_axis = linspace(0,duration,duration*fs);
T = 1/fs;           %time resolution

% Fundamental note
f0 = 65.4;          %[Hz]

% Boundary          
bound_b = 1000;     %normalized bridge impedance [Ohm/kg m^2 s]
bound_l = 1e20;     %normalized impedance of the left end [Ohm/kg m^2 s]

% String parameters
L = 1.92;           %length [m]
mass = 35e-3;       %mass [kg]
Te = 750;           %tension [N]
b1 = 0.5;           %air damping coefficient [Hz]
b2 = 6.25e-9;       %string internal friction coefficient [s]
k = 7.5e-6;         %stiffness [kg/s^2]
rho = mass/L;       %linear density [kg/m]
c = sqrt(Te/rho);   %speed stiffness [m/s]

% Spatial sampling parameters
gamma = fs/(2*f0);
Nmax = sqrt((-1+sqrt(1+16*k*gamma^2))/(8*k)); %maximum number of spatial samples

% Aliasing condition
% Number of maximum spatial steps
% prompt = 'Please, choose a number of spatial samples between 5 and ';
% N = 0;
% while N<5 || N>Nmax
%         N = input(strcat([prompt, num2str(floor(Nmax)), ':\n']));
% end

% Integer values
% Spatial sampling
N = floor(Nmax);                   %spacial samples
x_axis = linspace(0, L, N);
X = 1/N;                        %spacial resolution


% FD parameters
mu = k^2/(c^2*X^2);
nu = 2*b2*T/(X^2);
lambda = c*T/X;

% Hammer parameters
Mh = 4.9e-3;    %mass of the hammer [kg]
p = 2.3;        %stiffness exponent [\]
bh = 1e-4;      %fluid damping coeff [Hz]
K = 4e8;        %stiffness of the hammer [kg/s^2]
a = 0.12;       %position of the striking point in percentage
V_h0 = 2.5;     %initial velocity of hammer [m/s]

% Hammer contact window definition
w = 0.2;                        %width of window g
g_samples = floor(w*N);         %length of the hammer's window g(x,x0)
if mod(g_samples,2) ~= 0        %we want an odd number of samples
    g_sample = g_samples +1;    %to center the window on x0
end
delta = (g_samples-1)/2;       %half length of the window
g = hann(g_samples);            %window definition


%PDE Coefficients:

% string coefficients
aden = (1+b1*T);
a1 = (-lambda^2*mu) / aden;
a2 = (lambda^2+4*mu*lambda^2+nu) / aden;
a3 = (2-2*lambda^2-6*mu*lambda^2-2*nu) / aden;
a4 = (-1+b1*T+2*nu) / aden;
a5 = (-nu) / aden;
aF = (T^2/rho) / aden;

% Bridge boundary coefficients
bRden = (1+b1*T+bound_b*lambda);
bR1 = (2-2*mu*lambda^2-2*lambda^2) / bRden;
bR2 = (4*mu*lambda^2+2*lambda^2) / bRden;
bR3 = (-2*mu*lambda^2) / bRden;
bR4 = (-1+b1*T+bound_b*lambda) / bRden;
bRF = (T^2/rho) / bRden;

% Left hand (hinged string end) boundary coefficients
bLden = (1+b1*T+bound_l*lambda);
bL1 = (2-2*mu*lambda^2-2*lambda^2) / bLden;
bL2 = (4*mu*lambda^2+2*lambda^2) / bLden;
bL3 = (-2*mu*lambda^2) / bLden;
bL4 = (-1+b1*T+bound_b*lambda) / bLden;
bLF = (T^2/rho) / bLden;

% Hammer felt parameters
dden = (1+bh*T/(2*Mh));
d1 = 2 / dden;
d2 = (-1+bh*T/(2*Mh)) / dden;
dF = (-(T^2)/Mh) / dden;


%% Computation of the FD scheme
% Initialization

%displacement of the string
Y = zeros(length(t_axis), length(x_axis));
%force excerted by the hammer
F = zeros(length(t_axis), length(x_axis));
FH = zeros(length(t_axis),1);
%application of the window to the force axis
x0 = a*L;   %striking position [m]
m0 = find(abs(x_axis-x0)==min(abs(x_axis-x0)));
G = zeros(1,length(x_axis));
G(m0-delta:m0+delta) = G(m0-delta:m0+delta) + g';
%displacement of the hammer
eta = zeros(length(t_axis),1);
eta(2) = V_h0*T;


% Computation loop
for n = 2:(size(Y, 1) -1)      %iteration over time
    
            %force of the hammer
    if eta(n)>Y(n,m0)          %the hammer leaves the string
       FH(n) = K*abs(eta(n)-Y(n,m0))^p;
    else
        FH(n) = 0;
    end
    eta(n+1) = d1*eta(n) + d2*eta(n-1) + dF*FH(n);      %displacement of the hammer
    
    
    F(n,:) = FH(n)*G;
    
    clear m
    for m = 1:size(Y, 2)       %iteration over space for the displacement
        if m==1
            Y(n+1,m) = bL1*Y(n,m) + bL2*Y(n,m+1) + ...
                bL3*Y(n,m+2) + bL4*Y(n-1,m) + bLF*F(n,m);
        elseif m==2
            Y(n+1,m) = a1*(Y(n,m+2)-Y(n,m)+2*Y(n,m-1)) + a2*(Y(n,m+1)+Y(n,m-1)) + ...
                a3*Y(n,m) + a4*Y(n-1,m) + a5*(Y(n-1,m+1)+Y(n-1,m-1)) + aF*F(n,m);
        elseif m==(size(Y, 2))
            Y(n+1,m) = bR1*Y(n,m) + bR2*Y(n,m-1) + ...
                bR3*Y(n,m-2) + bR4*Y(n-1,m) + bRF*F(n,m);
        elseif m==(size(Y, 2) -1)
            Y(n+1,m) = a1*(2*Y(n,m+1)-Y(n,m)+Y(n,m-2)) + a2*(Y(n,m+1)+Y(n,m-1)) + ...
                a3*Y(n,m) + a4*Y(n-1,m) + a5*(Y(n-1,m+1)+Y(n-1,m-1)) + aF*F(n,m);
        else
            Y(n+1,m) = a1*(Y(n,m+2)+Y(n,m-2)) + a2*(Y(n,m+1)+Y(n,m-1)) + ...
                a3*Y(n,m) + a4*Y(n-1,m) + a5*(Y(n-1,m+1)+Y(n-1,m-1)) + aF*F(n,m);
        end
    end
end
%% Plot the displacement in time

figure(1);
i = 1;
while 1
    plot(x_axis, Y(i,:));
    xlabel('$x\,[m]$', 'interpreter', 'latex', 'fontsize', 17);
    xlim([0,L]);
    ylabel('$Displacement\,[m]$', 'interpreter', 'latex', 'fontsize', 17);
    ylim([-3e-4,3e-4]);
    title(strcat(['Time: ', num2str(i*T), ' s']), 'fontsize', 20)
    pause(T);
    if T*i > 8
        break
    end
    i = i+10;
end
%% Plot the synthesized signal play it and save it on the disk

% Play the sound



% Save on disk










