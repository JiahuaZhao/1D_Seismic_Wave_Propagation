%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 5.1 Exercises 4
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
% Discretizing Spatial domain
D = 100e+3; % Total depth, m
Emax = 1e+11; % Max Young’s module, Pa
Rhomax = 3.3e+3; % Max rock density, kg/m^3
vmax = sqrt(Emax / Rhomax); % Max theoretical wave velocity, m/s
C = 1; % CFL number
dx = vmax * dt / C; % Space interval, m
xnum = floor(D / dx + 1); % Grid points
x = 0:dx:D;
% Distribution of elastic parameters
layer = 40; % Interface
E = [5e+10*ones(1,round(xnum/100*layer)),1e+11*ones(1,round(xnum-xnum/100*layer))]; % Young’s module, Pa
Rho = [2.5e+3*ones(1,round(xnum/100*layer)), 3.3e+3*ones(1,round(xnum-xnum/100*layer))]; % Rock density, kg/m^3
v = sqrt(E./Rho); % m/s

% Plotting
figure(1)
subplot(1,3,1)
plot(E./1e+9,x./1000,'r','LineWidth',2)
axis([10, 200, 0, 100])
set(gca,'YDir','reverse')
ylabel('Depth (km)')
xlabel('Young’s module (GPa)')
title('Initial Elastic Model')
subplot(1,3,2)
plot(Rho,x./1000,'b','LineWidth',2)
axis([2e+3, 4e+3, 0, 100])
set(gca,'YDir','reverse')
ylabel('Depth (km)')
xlabel('Density (kg/m^3)')
title('Initial Density Model')
subplot(1,3,3)
plot(v,x./1000,'g','LineWidth',2)
axis([4e+3, 6e+3, 0, 100])
set(gca,'YDir','reverse')
ylabel('Depth (km)')
xlabel('Velocity (m/s)')
title('Initial Velocity Model')

%2) Initializing source - receiver location and synthetic seismograms 
source = int16(xnum * 0.5); % Equake. source
recv = 1; % Receiver
s_seismogram=zeros(1,nt); % Synthetic seismogram of receiver

%3) Initializing a displacement field and stresses
uc = zeros(1,xnum); % Current dis. field, u(t)
ub = uc; % u(t-1)
ua = ub; % u(t+1)
sc = zeros(1,xnum); % Current stress, s(t)

%4) Applying velocity-stress FD method method and plotting
% n = 0;
for j = 1:nt
%     n = n + 1;
    % Stress FD
    for i = 2:xnum-1
        sc(i) = (E(i+1)+E(i))/2 * (uc(i+1) - uc(i)) / dx;
        sc(i-1) = (E(i)+E(i-1))/2 * (uc(i) - uc(i-1)) / dx;
    end
    % The upper boundary condition is stress-free
    sc(1) = 0.0;
    
    % Displacement FD
    for i = 2:xnum-1
        ua(i) = dt^2 * (sc(i) - sc(i-1)) / (Rho(i) * dx) + 2 * uc(i) - ub(i);
    end
    % The lower boundary condition is motion-free
    ua(xnum) = 0.0;
    ua(1) = ua(2);
    
    % Applying the source function for 5 seconds
    tc = j * dt;
    if (tc <= tlen)
        ua(source) = 1e-3 * (sin(pi * tc / 5)).^2;
    end
    
    % Updating the displacement field and record the seismograms
    ub = uc;
    uc = ua;
    s_seismogram(1,j) = ua(recv);
    
    % Plotting
    figure(2)
    plot(ua.*1000,x./1000)
    hold on
    text(-4,10,'Layer I');
    text(-4,15,'E = 5e+10 Pa');
    text(-4,20,'Rho = 2.5e+3 kg/m^3');
    plot([-5,5],[layer,layer],'--k');
    text(-4,50,'Layer II');
    text(-4,55,'E = 1e+11 Pa');
    text(-4,60,'Rho = 3.3e+3 kg/m^3');
    axis([-5, 5, 0, 100]);
    set(gca,'YDir','reverse')
    xlabel('Amplitude (x10^{-3})');
    ylabel('Depth (km)')
    title(['Wavefield at ',num2str(tc),'s',' /',num2str(T),'s'])
    hold off
    
%     % Making a GIF figure
%     f = getframe(gcf);
%     imind = frame2im(f);
%     [imind,cm] = rgb2ind(imind,256);
%     if n == 1
%         imwrite(imind,cm,'1D_staggered_layer.gif','gif','LoopCount',inf,'DelayTime',0.01);
%     else
%         imwrite(imind,cm,'1D_staggered_layer.gif','gif','WriteMode','append','DelayTime',0.01);
%     end
end

% Plotting synthetic seismogram
figure(3)
plot(t,s_seismogram(1,:).*1000);
xlabel('Time (s)');
ylabel('Amplitude (x10^{-3})');
ylim([-5 5])
title('The Synthetic Seismogram at the Surface')

fprintf('When the seismic waves propagate in media with different elastic parameters, \nthe phenomenon of wave impedance will appear. \nWe can find from the dynamic image that when the seismic wave propagates upwards, \nthere will be a small fluctuation generated at the interface and propagated downward, \nwhich will eventually affect the record of the synthetic seismogram. \n');