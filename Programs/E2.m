%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 5.1 Exercises 2
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
E = 5e+10; % Young’s module, Pa
Rho = 3e+3; % Rock density, kg/m^3
v = sqrt(E / Rho); % Theoretical wave velocity, m/s
C = 1; % CFL number
dx = v * dt / C; % Space interval, m
xnum = floor(D / dx + 1); % Grid points
x = 0:dx:D;

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
%n = 0;
for j = 1:nt
%     n = n + 1;
    % Stress FD
    for i = 2:xnum-1
        sc(i) = E * (uc(i+1) - uc(i)) / dx;
        sc(i-1) = E * (uc(i) - uc(i-1)) / dx;
    end
    % The upper boundary condition is stress-free
    sc(1) = 0.0;
    
    % Displacement FD
    for i = 2:xnum-1
        ua(i) = dt^2 * (sc(i) - sc(i-1)) / (Rho * dx) + 2 * uc(i) - ub(i);
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
    figure(1);
    plot(ua.*1000,x./1000);
    axis([-2, 2, 0, 100]);
    set(gca,'YDir','reverse');
    xlabel('Amplitude (x1000)');
    ylabel('Depth / km');
    title(['Wavefield at ',num2str(tc),'s',' /',num2str(T),'s']);
    
%     % Making a GIF figure
%     f = getframe(gcf);
%     imind = frame2im(f);
%     [imind,cm] = rgb2ind(imind,256);
%     if n == 1
%         imwrite(imind,cm,'1D_staggered.gif','gif','LoopCount',inf,'DelayTime',0.01);
%     else
%         imwrite(imind,cm,'1D_staggered.gif','gif','WriteMode','append','DelayTime',0.01);
%     end
end

% Plotting synthetic seismogram
figure(2);
plot(t,s_seismogram(1,:).*1000);
xlabel('Time / s');
ylabel('Amplitude');
ylim([-2 2]);
title('The Synthetic Seismogram at the Surface')

% % Calculating arrivel time
% for k = 1:numel(s_seismogram)
%     if s_seismogram(1,k)*1000 > 0.01
%         st = s_seismogram(1,k)*1000;
%         at = k * dt;
%         break;
%     end
% end
% 
% % Plotting picking
% figure(3);
% plot(t,s_seismogram(1,:).*1000,'g'); hold on;
% plot([at,at],[-2,2],'--'); hold on;
% plot([0,T],[st,st],'k');
% xlabel('Time / s');
% ylabel('Amplitude');
% ylim([-2 2]); hold off;
% 
% c = double(D / 2 / at);
% fprintf('The velocity of the propagating wave is %d m/s. \n',c);