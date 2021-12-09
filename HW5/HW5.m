close all
clear all
clc


%% Data for the plate

%dimensions of the plate
L_x = 1.4;
L_y = 1;
h = 0.01;

%Longitudinal parameters for Sitka spruce
rho = 404.05;                  % density [kg/m^3]
E = 12571*(10^6);              % Young modulus [Pa]
nu = 0.03; %0.37               % poisson ratio [\]
B = (E*h^3)/(12*(1-nu^2));     % bending stiffness

% frequencies of the pairs of strings
f_F4 = 349.23;    %[Hz]
f_A4 = 440.00;
f_C5 = 523.25;
f_E5 = 659.25;
f_G5 = 783.99;
notes = [f_F4, f_A4, f_C5, f_E5, f_G5];
notes_names = ['F4', 'A4', 'C5', 'E5', 'G5'];

%% Soundboard characterization

space_res = 100;             % resolution of x and y axis
x_axis = linspace(0,L_x,space_res);
y_axis = linspace(0,L_y,space_res);
[x_grid,y_grid] = meshgrid(x_axis,y_axis);

f = linspace(10,10000,10000);       % frequency axis [Hz]
omega = f*(2*pi);                   % [rad/s]

% Modal approach
max_nm = 16;          % max value to give to m and n
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
surf(x_grid,y_grid,Phi_nm{1,3}) % plot a mode shape


%% Input impedance computation of single point as a function of frequency

%considered point
x_index = 60;
y_index = 50;

%damping matrix with reyleigh damping
alpha = 0.01;
beta = 0.00001;
%damping matrix with isotropic loss factor
eta = 0.012;

Y = zeros(1, length(f));

for i = 1:max_nm
    for j = 1:max_nm
        m_mod = rho*h*L_x*L_y/4;
        k_mod = m_mod*omega_nm(i,j)^2;
        c_mod = alpha*m_mod + beta*k_mod;
        Y_num = 1i*omega*Phi_nm{i,j}(x_index,y_index)^2;
        %Y_den = -omega.^2*m_mod + 1i*omega*c_mod + k_mod;
        Y_den = -omega.^2*m_mod + (1+1i*eta)*k_mod;
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
f_value = notes;
f_index = zeros(1, length(notes));
for i = 1:length(notes)
    f_index(i) = find(abs(f-f_value(i))==min(abs(f-f_value(i))),1);
end

Y = zeros(length(x_axis), length(y_axis));
Z_notes = cell(1,length(notes));

for ii = 1:length(notes)
    for i = 1:max_nm
        for j = 1:max_nm
            m_mod = rho*h*L_x*L_y/4;
            k_mod = m_mod*omega_nm(i,j)^2;
            c_mod = alpha*m_mod + beta*k_mod;
            Y_num = 1i*omega(f_index(ii))*Phi_nm{i,j}'.*Phi_nm{i,j};
            %Y_den = -omega(f_index)^2*m_mod + 1i*omega(f_index)*c_mod + k_mod;
            Y_den = -omega(f_index(ii))^2*m_mod + (1+1i*eta)*k_mod;
            Y = Y + Y_num./Y_den;
        end
    end
    Z_notes{ii} = 1./Y;
end
%% surf plot
note = 1;  % 1==F2, 2==A4, 3==C5, 4==E5, 5==G5
Z = Z_notes{note};

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
note = 5;  % 1==F2, 2==A4, 3==C5, 4==E5, 5==G5
Z = Z_notes{note};

figure()
pcolor(x_grid,y_grid, db(abs(Z)))
xlabel('$x\,[m]$', 'interpreter', 'latex', 'fontsize', 17)
ylabel('$y\,[m]$', 'interpreter', 'latex', 'fontsize', 17)
title('$|Z_{INPUT}|\,[dB]$', 'interpreter', 'latex', 'fontsize', 20)
colorbar


%% Bridge design

x_F4 = 0.71; y_F4 = 0.482;
x_A4 = 0.77; y_A4 = 0.392;
x_C5 = 0.81; y_C5 = 0.350;
x_E5 = 0.86; y_E5 = 0.300;
x_G5 = 0.93; y_G5 = 0.232;
x_notes_coord = [x_F4, x_A4, x_C5, x_E5, x_G5];
y_notes_coord = [y_F4, y_A4, y_C5, y_E5, y_G5];

x_notes_index = zeros(1,5);
y_notes_index = zeros(1,5);

for i = 1:5
    x_notes_index(i) = find(abs(x_axis-x_notes_coord(i))==min(abs(x_axis-x_notes_coord(i))));
    y_notes_index(i) = find(abs(y_axis-y_notes_coord(i))==min(abs(y_axis-y_notes_coord(i))));
end

%%
for i = 1:length(notes)
    figure()
    hold on
    pcolor(x_grid,y_grid, db(abs(Z_notes{i})))
    plot(x_axis(x_notes_index(i)), y_axis(y_notes_index(i)), 'ro', 'linewidth', 2)
    hold off
    xlabel('$x\,[m]$', 'interpreter', 'latex', 'fontsize', 17)
    ylabel('$y\,[m]$', 'interpreter', 'latex', 'fontsize', 17)
    title('$|Z_{INPUT}|\,[dB]$', 'interpreter', 'latex', 'fontsize', 20)
    colorbar
end

%%

figure()
hold on
plot(x_axis,zeros(1,length(x_axis)), 'k', 'linewidth', 3)
plot(x_axis,L_y*ones(1,length(x_axis)), 'k', 'linewidth', 3)
plot(zeros(1,length(y_axis)),y_axis, 'k', 'linewidth', 3)
plot(L_x*ones(1,length(y_axis)),y_axis, 'k', 'linewidth', 3)
for i = 1:length(notes)
    if i<5
        line([x_axis(x_notes_index(i)), x_axis(x_notes_index(i+1))],[y_axis(y_notes_index(i)), y_axis(y_notes_index(i+1))], 'color', 'k', 'linewidth', 3)
    end
    plot(x_axis(x_notes_index(i)), y_axis(y_notes_index(i)), 'ro', 'linewidth', 2)
    text(x_axis(x_notes_index(i)), y_axis(y_notes_index(i))+0.03, strcat(notes_names(2*i-1),notes_names(2*i)), 'fontsize', 15)
end
hold off
xlim([-0.1, 1.5])
ylim([-0.1, 1.1])




