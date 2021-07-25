%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 5.1 Exercises 3
% Author: Jiahua Zhao
% Last Modified: 10/06/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;

%1) Initializing parameters
% Discretizing temporal domain
T = 100; % Total time, s
dt = 0.1; % Timestep, s
tlen = 5; % source duration, s
t = 0:dt:(T-dt);
nt = numel(t);
% Spatial domain
D = 100e+3; % Total depth, m
E = linspace(1e+10,1e+11,200); % Young’s module, Pa
Rho = linspace(1e+3,4e+3,200); % Rock density, kg/m^3
C = 1; % CFL Number

%2) Initializing receiver location and synthetic seismogram
recv = 1; % Receiver
s_seismogram=zeros(1,nt); % Synthetic seismogram of receiver

%3) Applying velocity-stress FD method method and plotting
c = zeros(1,numel(E));
for i = 1:numel(E)
    v1(i) = sqrt(5e+10 / Rho(i));
    dx1(i) = v1(i) .* dt / C;
    xnum1(i) = floor(D / dx1(i) + 1);
    % Initializing a displacement field and stresses
    uc1 = zeros(1,xnum1(i));  ub1 = uc1; ua1 = ub1; sc1 = zeros(1,xnum1(i));
    c1(i) = FD_staggered(ua1,ub1,uc1,sc1,5e+10,Rho(i),c(i),tlen,s_seismogram,recv,D,dx1(i),dt,nt,xnum1(i));
    
    v2(i) = sqrt(E(i) / 3e+3);
    dx2(i) = v2(i) .* dt / C;
    xnum2(i) = floor(D / dx2(i) + 1);
    % Initializing a displacement field and stresses
    uc2 = zeros(1,xnum2(i));  ub2 = uc2; ua2 = ub2; sc2 = zeros(1,xnum2(i));
    c2(i) = FD_staggered(ua2,ub2,uc2,sc2,E(i),3e+3,c(i),tlen,s_seismogram,recv,D,dx2(i),dt,nt,xnum2(i));
    
    for j = 1:numel(Rho)
        v3(i) = sqrt(E(i) / Rho(j));
        dx3(i) = v3(i) .* dt / C;
        xnum3(i) = floor(D / dx3(i) + 1);
        % Initializing a displacement field and stresses
        uc3 = zeros(1,xnum3(i));  ub3 = uc3; ua3 = ub3; sc3 = zeros(1,xnum3(i));
        c3(i,j) = FD_staggered(ua3,ub3,uc3,sc3,E(i),Rho(j),c(i),tlen,s_seismogram,recv,D,dx3(i),dt,nt,xnum3(i));
    end
end

% Plotting
subplot(2,2,1);
plot(Rho./1000,c1);
xlabel('Density (kg/m^3)');
ylabel('Velocity (m/s)');
title('E = 50 GPa');
grid on;

subplot(2,2,2); plot(E./1e+9,c2);
xlabel('Young’s module (GPa)');
ylabel('Velocity (m/s)');
title('Rho = 3 kg/m^3');
grid on;

subplot(2,2,[3,4]); pcolor(Rho./1000,E./1e+9,c3);
xlabel('Density (kg/m^3)');
ylabel('Young’s module (GPa)');
shading interp;
colorbar; colormap(jet);
ylabel(colorbar,'Velocity (m/s)');
hold on;
[cs, h]=contour(Rho./1000,E./1e+9,c3,[2000,3000,4000,5000,6000,7000,8000],'w--');
clabel(cs, h,'LabelSpacing',100, 'FontSize', 10, 'Color', 'k');
hold off;