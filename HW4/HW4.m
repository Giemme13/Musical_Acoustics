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

%% Mechanical characterization

k=m*(2*pi*f0)^2;
%Q to find resistance