%%
clear all
close all
clc

%% DATA

m=0.0455;
f0=49;
D=0.07;
L=0.247;
V=0.20;
Q=50;
rho=1.225;
c=343;

%% Mechanical characterization

k=m*(2*pi*f0)^2;
R=2*pi*f0*m/Q;

%% Equivalent circuit

L1=rho*(L+0.85*D); %inductor for the mass of the tube
C1=V/(rho*c^2); %condenser for the air volume
L2=M; %inductor for the mass of the loudspeaker
C2=1/k; %condenser for the equivalent spring of the loudspeaker
R2=R; %resistance for the damping of the loudspeaker