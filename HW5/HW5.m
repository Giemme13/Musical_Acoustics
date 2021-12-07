close all
clear all
clc

%% Data for the plate

%dimensions of the plate
L_x = 1;
L_y = 1.4;
h = 0.01;

%Longitudinal parameters for Sitka spruce
rho = 400;                  % density [kg/m^3]
E = 12*(10^9);              % Young modulus [Pa]
nu = 0.37;                  % poisson ratio [\]
B = (E*h^3)/(12*(1-nu^2));  % bending stiffness

% frequencies of the pairs of strings
f_F2 = 349.23;    %[Hz]
f_A4 = 440;
f_C5 = 523.25;
f_E5 = 659.25;
f_G5 = 783.99;
notes = [f_F2, f_A4, f_C5, f_E5, f_G5];

%% Soundboard characterization

space_res = 100;             % resolution of x and y axis
x_axis = linspace(0,L_x,space_res);
y_axis = linspace(0,L_y,space_res);
[x_grid,y_grid] = meshgrid(x_axis,y_axis);

f = linspace(10,4000,10000);        % frequency axis [Hz]
omega = f*(2*pi);                   % [rad/s]

% Modal approach
max_nm = 10;          % max value to give to m and n
n = 1:1:max_nm;
m = 1:1:max_nm;
% wavenumbers
k_xn = pi*n/L_x;
k_ym = pi*m/L_y;
[K_xn, K_ym] = meshgrid(k_xn, k_ym);
k_nm = sqrt(K_xn.^2+K_ym.^2);
% resonances
omega_nm = k_nm.^2*sqrt(B/(rho*h));
f_nm = omega_nm/(2*pi);
% mode shapes initialization
Phi_nm = cell(max_nm,max_nm);
toll = 1e-5;    % minimum value of a displacement of a point
                % in order for it to be considered non zero

% computation of the mode shapes for every point of the soundboard
for i = 1:max_nm
    for j = 1:max_nm
        Phi_nm{i,j} = sin(K_xn(i,j).*x_grid).*sin(K_ym(i,j).*y_grid);
    end
end
% don't know why, but if i place the following inside
% the previous cycle is doesn't work for m = 10
for i = 1:max_nm
    for j = 1:max_nm
        Phi_nm{i,j}(abs(Phi_nm{j,j})<toll) = 0;
    end
end
%%
surf(x_grid,y_grid,Phi_nm{2,4}) % plot a mode shape


%% Input impedance computation of single point as a function of frequency

%considered point
x_index = 60;
y_index = 50;

%damping matrix with reyleigh damping
alpha = 0.0;
beta = 0.0000;
%damping matrix with isotropic loss factor
eta = 0.3;

Y = zeros(1, length(f));

for i = 1:max_nm
    for j = 1:max_nm
        m_mod = rho*h*L_x*L_y/4;
        k_mod = m_mod*omega_nm(i,j)^2;
        c_mod = alpha*m_mod + beta*k_mod;
        Y_num = 1i*omega*Phi_nm{i,j}(x_index,y_index)^2;
        Y_den = -omega.^2*m_mod + 1i*omega*c_mod + k_mod;
        %Y_den = -omega.^2*m_mod + (1+1i*eta)*k_mod;
        Y = Y + Y_num./Y_den;
    end
end
Z = 1./Y;
%%
figure()
subplot(2,1,1)
plot(f, db(abs(Z)), 'linewidth', 1.5)
grid on
xlabel('Frequency [Hz]', 'fontsize', 17)
ylabel('$|Z_{INPUT}|\,[dB]$', 'interpreter', 'latex', 'fontsize', 17)
title('Magnitude', 'fontsize', 20)
subplot(2,1,2)
plot(f, angle(Z), 'linewidth', 1.5)
grid on
xlabel('Frequency [Hz]', 'fontsize', 17)
ylabel('$\angle Z_{INPUT}\,[deg]$', 'interpreter', 'latex', 'fontsize', 17)
yticks([-pi,-pi/2,0,pi/2,pi])
yticklabels({'-180','-90','0','90','180'})
title('Phase', 'fontsize', 20)


%% Input impedance computation of single frequency as a function of coordinates

%considered frequency
f_value = notes(4);
f_index = find(abs(f-f_value)==min(abs(f-f_value)),1);

Y = zeros(length(x_axis), length(y_axis));

for i = 1:max_nm
    for j = 1:max_nm
        m_mod = rho*h*L_x*L_y/4;
        k_mod = m_mod*omega_nm(i,j)^2;
        c_mod = alpha*m_mod + beta*k_mod;
        Y_num = 1i*omega(f_index)*Phi_nm{i,j}'.*Phi_nm{i,j};
        Y_den = -omega(f_index)^2*m_mod + 1i*omega(f_index)*c_mod + k_mod;
        %Y_den = -omega.^2*m_mod + (1+1i*eta)*k_mod;
        Y = Y + Y_num./Y_den;
    end
end
Z = 1./Y;
%% surf plot
figure()
subplot(2,1,1)
surf(x_grid,y_grid, db(abs(Z)))
xlabel('$x\,[m]$', 'interpreter', 'latex', 'fontsize', 17)
ylabel('$y\,[m]$', 'interpreter', 'latex', 'fontsize', 17)
zlabel('$|Z_{INPUT}|\,[dB]$', 'interpreter', 'latex', 'fontsize', 17)
title('Magnitude', 'fontsize', 20)
colorbar
subplot(2,1,2)
surf(x_grid,y_grid, angle(Z))
xlabel('$x\,[m]$', 'interpreter', 'latex', 'fontsize', 17)
ylabel('$y\,[m]$', 'interpreter', 'latex', 'fontsize', 17)
zlabel('$\angle Z_{INPUT}\,[deg]$', 'interpreter', 'latex', 'fontsize', 17)
zticks([-pi,-pi/2,0,pi/2,pi])
zticklabels({'-180','-90','0','90','180'})
title('Phase', 'fontsize', 20)

%% pcolor plot
figure()
pcolor(x_grid,y_grid, db(abs(Z)))
xlabel('$x\,[m]$', 'interpreter', 'latex', 'fontsize', 17)
ylabel('$y\,[m]$', 'interpreter', 'latex', 'fontsize', 17)
title('$|Z_{INPUT}|\,[dB]$', 'interpreter', 'latex', 'fontsize', 20)
colorbar

