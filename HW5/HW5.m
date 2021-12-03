close all
clear all
clc

%% Data for the plate

L_x=1;
L_y=1.4;
h=0.01;

%Longitudinal parameters for Sitka spruce
mu=416.21; % density (average) [kg/m^3] http://www.conforg.fr/isma2014/cdrom/data/articles/000110.pdf
E=10800*(10^6); %Young modulus https://amesweb.info/Materials/Youngs-Modulus-of-Wood.aspx
%nu=; %poisson ratio 