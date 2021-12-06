close all
clear all
clc

%% Data for the plate

%dimensions of the plate
L_x = 1;
L_y = 1.4;
h = 0.01;

%Longitudinal parameters for Sitka spruce
rho = 400;                  % density (average) [kg/m^3] http://www.conforg.fr/isma2014/cdrom/data/articles/000110.pdf
E = 12*(10^9);              % Young modulus https://amesweb.info/Materials/Youngs-Modulus-of-Wood.aspx
nu = 0.37;                  % poisson ratio
B = (E*h^3)/(12*(1-nu^2));  % bending stiffness

% frequencies of the pairs of strings
f_F2 = 349.23;
f_A4 = 440;
f_C5 = 523.25;
f_E5 = 659.25;
f_G5 = 783.99;


%% Soundboard characterization

space_res = 100;
x = linspace(0,L_x,space_res);      % x_axis
y = linspace(0,L_y,space_res);      % y_axis
[X,Y] = meshgrid(x,y);          % grids

f = linspace(20,10000,10000);    % frequency axis
omega = f*(2*pi);

% Modal approach
max_mn = 10;       % max value to give to m and n
n = 1:1:max_mn;
m = 1:1:max_mn;
% wavenumbers
k_xn = pi*n/L_x;
k_ym = pi*m/L_y;
[K_xn, K_ym] = meshgrid(k_xn, k_ym);
k_nm = sqrt(K_xn.^2+K_ym.^2);
% resonances
omega_nm = k_nm.^2*sqrt(B/(rho*h));
f_nm = omega_nm/(2*pi);
% mode shapes initialization
Phi_nm = cell(max_mn,max_mn);

% computation of the mode shapes for every point of the soundboard
for i = 1:max_mn
    for j = 1:max_mn
        Phi_nm{i,j} = sin(K_xn(i,j).*X).*sin(K_ym(i,j).*Y);
    end
end
%%
surf(X,Y,Phi_nm{2,2}) % plotta uno dei mode shapes


%% Input impedance computation of single point

x_coord = 50;
y_coord = 70;
alpha = 0.0;
beta = 0.000;
eta = 0.69;
Z = zeros(1, length(f));

for i = 1:max_mn
    for j = 1:max_mn
        m_mod = rho*h*L_x*L_y/4;
        k_mod = m_mod*omega_nm(i,j)^2;
        c_mod = alpha*m_mod + beta*k_mod;
        Z_num = -omega.^2*m_mod + 1i*omega*c_mod + k_mod;
        %Z_num = -omega.^2*m_mod + (1+1i*eta)*k_mod;
        Z_den = 1i*omega*Phi_nm{i,j}(x_coord,y_coord)^2;
        Z = Z + Z_num./Z_den;
    end
end
%%
figure(1)
subplot(2,1,1)
plot(f, db(abs(Z)))
subplot(2,1,2)
plot(f, angle(Z))


%% Input impedance computation of single frequency

Z = cell(1,length(f));

for k = 1:length(f)
    for i = 1:res_num
        m_mod = rho*h*L_x*L_y/4;
        k_mod = m_mod*omega_nm(i)^2;
        csi = 0.7;
        c_mod = csi*2*m_mod*omega_nm(i);
        Z_num = -omega(k).^2*m_mod + 1i*omega(k)*c_mod + k_mod;
        Z_den = 1i*omega(k)*Phi_nm{i}'.*Phi_nm{i};
        Z{1,k} = Z_num./Z_den;
    end
end
%%
figure(2)
subplot(2,1,1)
surf(X,Y, db(abs(Z{1,900})))
subplot(2,1,2)
surf(X,Y, angle(Z{1,401}))
