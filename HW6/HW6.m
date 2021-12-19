clear;
close all;
clc;

%% Data

c=343;
theta=0.75; %cone semiangle
L=0.45; %length
F4=349.23; %frequency with hole closed
G4=392; %frequency for last hole open
A4=440; %frequency for both holes open

r_foot=L*tand(theta); 
d_foot=2*r_foot; %diameter at foot of resonator

L_tot=L+0.6*r_foot; %end correction of tube

x=c/(2*F4)-L_tot; %distance of mouth from cone vertex

r_mouth=x*tand(theta);
d_mouth=2*r_mouth; %diameter at mouth of resonator
