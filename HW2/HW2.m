%%
clear all
close all
clc

%% Circular Membrane

% Data
a = 0.2;        %radius  [m]
T = 10;         %tension [N/m]
sigma = 0.1;    %surface weight [kg/m^2]

% a) propagation speed
c = sqrt(T/sigma);

% b) bessel function values
J_zeros = [2.4048, 3.8317, 5.1356, 5.5201, 6.3802];      % vector of bessel function zero-values
memb_freqs = zeros(length(J_zeros),1);           % vector of natural frequencies

for i = 1:length(J_zeros)
    memb_freqs(i) = J_zeros(i)*c/(2*pi*a);       % computation of the natural frequencies
end

% Mode shapes
m = [0 1 2 0 3];    % m values
n = [1 1 1 2 1];    % n values
res = 50;           % resolution of plots

modes_membrane = cell(1,5);  % collection of mode shapes

for i = 1:5
    modes_membrane{1,i} = zeros(res,res);  % initialization
end

radius = linspace(0, a, res);                % radius interval [0, 0.2] (for the plots)
theta = linspace(0, 2*pi, res);              % angle interval [0, 2pi] (for plots and computation)
[Radius, Theta] = meshgrid(radius, theta);   % grids for the mesh of the plots

for i=1:5
    r = linspace(0, J_zeros(i), res);    % radius vector (for the bessel function)
    
    [R, T] = meshgrid(r, theta);         % grids for the computation of the mode shapes
    
    modes_membrane{1,i} = besselj(m(i), R) .* cos(m(i)*T);   % computation of the natural modes
                        % we use the grid R so that we consider the
                        % boundary condition: z = 0 at radius = r
    
    [X,Y,Z] = pol2cart(Theta, Radius, modes_membrane{1,i});     % convertion from polar to cartesian coordinates
                        % we use the grid Radius so that the plots will be
                        % contained in a circle of radius 0.2
    
    figure(i)
    subplot(1,2,1)
    pcolor(X,Y,Z)
    colorbar
    subplot(1,2,2)
    surf(X,Y,Z, 'FaceAlpha',0.8)
    sgtitle(strcat(['Mode ', num2str(i), ' - ', num2str(memb_freqs(i)), ' Hz']), 'fontsize', 30)
end

%% Circular Plate

% Data
a = 0.2;      % radius  [m]
h = 0.001;    % thickness [m]
E = 69e9;     % Young's modulus [GPa]
rho = 2700;   % volume density [kg/m^3]
nu = 0.334;   % Poisson's ratio [-]

% c) propagation speed
cQL = sqrt(E/(rho*(1-nu^2)));               % quasi-longitudinal waves
cL  = sqrt(E*(1-nu)/(rho*(1+nu)*(1-2*nu)));  % longitudinal waves

% d) propagation speed as function of frequency 
f = 0:0.1:500;
v = sqrt(1.8*h*cQL.*f);

figure(6)
plot(f,v)
xlabel('$f [Hz]$', 'interpreter', 'latex', 'fontsize', 20)
ylabel('$c_B [\frac{m}{s}]$', 'interpreter', 'latex', 'fontsize', 20)
title('Dispersion of bending waves velocity', 'fontsize', 20)
grid on

%e) frequencies

coeff = [0.4694, 2.08, 3.41, 3.89, 5.00];  % coefficients multiplying f01
plate_freqs = zeros(1, 5);                 % computation of the natural frequencies

for i = 1:5
    if i == 1
        plate_freqs(i) = coeff(i)*cQL*h/(a^2);
    else
        plate_freqs(i) = coeff(i)*plate_freqs(1);
    end
end

k = [3.189, 4.612, 5.904, 6.308, 7.143];

%f) Mode shapes
m = [0 1 2 0 3];    % m values
n = [1 1 1 2 1];    % n values

res = 50;           % resolution of plots

modes_plate = cell(1,5);  % collection of mode shapes

for i = 1:5
    modes_plate{1,i} = zeros(res,res);  % initialization
end

radius = linspace(0, a, res);                             % radius interval [0, 0.2] (for the plots)
theta = linspace(0, 2*pi, res);                           % angle interval [0, 2pi] (for plots and computation)
[Radius, Theta] = meshgrid(radius, theta);                % grids for the mesh of the plots

      % mode shape m,n (r,phi) = cos(m*phi + alpha) * (A*(Jm(kr) - j^(-m)*Jm(jkr)*(J'm(ka)/I'm(ka)) , impose A = 1 & alpha = 0.

for i = 1:5
    r = linspace(0, k(i), res);    % radius vector (for the bessel function)
    
    [R, T] = meshgrid(r, theta);         % grids for the computation of the mode shapes
    
    ratio = besselj(m(i), k(i))/besseli(m(i), k(i));
    
    modes_plate{1,i} = ( besselj(m(i), R) - ratio*besseli(m(i), R)) .* cos(m(i)*T);   % computation of the natural modes
                        % we use the grid R so that we consider the
                        % boundary condition: z = 0 at radius = a
    
                        
    [X,Y,Z] = pol2cart(Theta, Radius, modes_plate{1,i});     % convertion from polar to cartesian coordinates
                        % we use the grid Radius so that the plots will be
                        % contained in a circle of radius 0.2
    
                        
    figure(6+i)
    subplot(1,2,1)
    pcolor(X,Y,Z)
    colorbar
    subplot(1,2,2)
    surf(X,Y,Z, 'FaceAlpha',0.8)
    sgtitle(strcat(['Mode ', num2str(i), ' - ', num2str(plate_freqs(i)), ' Hz']), 'fontsize', 20)
end

%% Coupled system

% Data
r = 0.2;                   % plate radius [m]
h = 0.001;                 % thickness [m]
rho_p = 2700;              % plate volume density [kg/m^3]
Q = 50;                    % plate merit factor
r_str = 0.001;             % string section radius [m]
L = 0.4;                   % string length [m]
rho_s = 5000;              % string volume density [kg/m^3]
mu = rho_s*pi*(r_str^2);   % string linear density [kg/m]

% g) tension on the string assuming both ends fixed
string_f0 = plate_freqs(1);  % the string is tunes to the fundamental of the plate
cS = 2*L*string_f0;          % wave speed on the string [m/s]
T = (cS^2)*mu;               % tension [N]

% h) coupling frequencies   (condition m/(n^2*M) = pi^2/(4Q^2)
m = L*mu;                 % string mass
M = rho_p*(h*(r^2)*pi);   % plate mass
RH = (pi^2)/(4*(Q^2));    % right hand side of condition
str_freqs = [string_f0, string_f0*2, string_f0*3, string_f0*4, string_f0*5]; % string frequencies
LH = zeros(1,5);          % left hand side of condition (to be computed)


diff = zeros(1,5);        % normalized differences between resonance frequencies of the string and of the plate
for i=1:5
    LH(i) = m /((i^2)*M);   % left hand side of condition computation
    diff(i) = (str_freqs(i) - plate_freqs(i))/plate_freqs(i);  % calculation of the normalized differences
    disp(strcat(['Mode ', num2str(i), ' ratio: ', num2str(diff(i))]))
end


% LH<RH weak coupling, LH>RH strong
for i = 1:5
    if LH(i)<RH    % weak coupling, frequencies unaltered
        disp(strcat(['Mode ', num2str(i), ': weak coupling']))
    else           % strong coupling, frequencies altered
        disp(strcat(['Mode ', num2str(i), ': strong coupling']))
    end
end

% by inspection on graph - values on x-axis
scale_high = [0.02, 0.009, 0.002, 0.04, 0];
scale_low = [-0.02, -0.045, -0.12, -0.011, -0.01];

% computation of coupling frequencies
coup_low = scale_low .* plate_freqs + plate_freqs;
coup_high = scale_high .*plate_freqs + plate_freqs;
