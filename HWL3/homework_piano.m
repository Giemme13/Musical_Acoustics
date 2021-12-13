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
length=8; %signal is 8 second long
T=44100; %time resolution
t=linspace(0,length,length*T);

% Fundamental note
f=65.4;

% Boundary          


% String parameters
m_s=35e-3;%mass 
k=7.5e-6; %stiffness
l=1.92; %length
Te=750; %tension
u=m_s/l; %linear density
c=sqrt(Te/u); %speed stiffness, 

% Spatial sampling parameters
w=0.2; %width of window g
gamma=T/(2*f);
Nmax=sqrt((-1+sqrt(1+16*k*gamma^2))/(8*k));
X=floor(Nmax); %spacial resolution

% Aliasing condition
% Number of maximum spatial steps

% Integer values
% Spatial sampling



% FD parameters
b1=0.5; %air damping coeff.
b2=6.25e-9; %internal friction of the string
bh=1e-4; %fluid damping coeff 
mu=k^2/(c^2*X^2);
v=2*b2*T/(X^2);
lambda=c*T/X;
a1=(-lambda^2*mu)/(1+b1*T);
a2=(lambda^2+4*mu*lambda^2+v)/(1+b1*T);
a3=(2-2*lambda^2-6*mu*lambda^2-2*v)/(1+b1*T);
a4=(-1+b1*T+2*v)/(1+b1*T);
a5=(-v)/(1+b1*T);
aF=0.12;

% Hammer parameters
Mh=4.9e-3; %mass of the hammer
K=4e8; %stiffness of the hammer
V_h0=2.5; %initial velocity of hammer

% Hammer contact window definition
d1=2/(1+bh*T/(2*Mh));
d2=(-1+bh*T/(2*Mh))/(1+bh*T/(2*Mh));
dF=(-(T^2)/Mh)/(1+bh*T/(2*Mh));

%PDE Coefficients:

% Bridge boundary coefficients
bound_b=1000;
bR1=(2-2*mu*lambda^2-2*lambda^2)/(1+b1*T+bound_b*lambda);
bR2=(4*mu*lambda^2+2*lambda^2)/(1+b1*T+bound_b*lambda);
bR3=(-2*mu*lambda^2)/(1+b1*T+bound_b*lambda);
bR4=(-1-b1*T+bound_b*lambda)/(1+b1*T+bound_b*lambda);
bRF=(T^2/u)/(1+b1*T+bound_b*lambda);

% Left hand (hinged string end) boundary coefficients
bound_l=1e20;
bL1=(2-2*mu*lambda^2-2*lambda^2)/(1+b1*T+bound_l*lambda);
bL2=(4*mu*lambda^2+2*lambda^2)/(1+b1*T+bound_l*lambda);
bL3=(-2*mu*lambda^2)/(1+b1*T+bound_l*lambda);
bL4=(-1-b1*T+bound_l*lambda)/(1+b1*T+bound_l*lambda);
bLF=(T^2/u)/(1+b1*T+bound_l*lambda);

% Hammer felt parameters
p=2.3; %stiffness exponent

%% Computation of the FD scheme
% Initialization

% Computation loop
%% Plot the displacement in time

%% Plot the synthesized signal play it and save it on the disk

% Play the sound



% Save on disk










