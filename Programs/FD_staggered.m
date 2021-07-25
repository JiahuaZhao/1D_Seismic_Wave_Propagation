function [c] = FD_staggered(ua,ub,uc,sc,E,Ro,c,tlen,s_seismogram,recv,D,dx,dt,nt,xnum)

for j = 1:nt
    % Stress FD
    for i = 2:xnum-1
        sc(i) = E * (uc(i+1) - uc(i)) / dx;
        sc(i-1) = E * (uc(i) - uc(i-1)) / dx;
    end
    % The upper boundary condition is stress-free
    sc(1) = 0.0;
    
    % Displacement FD
    for i = 2:xnum-1
        ua(i) = dt^2 * (sc(i) - sc(i-1)) / (Ro * dx) + 2 * uc(i) - ub(i);
    end
    % The lower boundary condition is motion-free
    ua(xnum) = 0.0;
    ua(1) = ua(2);
    
    % Applying the source function for 5 seconds
    tc = j * dt;
    if (tc <= tlen)
        ua(int16(xnum/2)) = 1e-3 * (sin(pi * tc / 5)).^2;
    end
    % Updating the displacement field and record the seismograms
    ub = uc;
    uc = ua;
    s_seismogram(1,j) = ua(recv);
end

% Calculating arrivel time
for k = 1:numel(s_seismogram)
    if s_seismogram(1,k)*1000 > 0.1
        at = k * dt;
        c = double(D / 2 / at);
        break;
    end
end
end