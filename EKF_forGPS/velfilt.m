function [vx_his, vy_his, vx_gps, vy_gps] = velfilt(accfile, gpsfile)
    %% Import file
%     [biasx, biasy, biasz,gpsbias, station] = biascalc('accelerometer_still.csv','gps_still.csv');
    acc = csvread('accelerometer.csv',1);
%     acc = csvread(accfile,1);
    acc = acc(:,[2,4]).';
%     acc = [acc(:,2); acc(:,4)];
%     acc = [acc(:,2).'-biasx; acc(:,4).'-biasz].*9.81/station;
%     gps = csvread(gpsfile,1);
    gps = csvread('gps.csv',1);
    gps = gps(:,6).';%-gpsbias;
    dt = 0.05;
    t = 0:dt:dt*size(acc,2);
    %% low pass filter
    f0=0.25; %cut-off frequency
    w0=2*pi*f0;
    N=512;
    Fs=20; % sampling frequency
    [NUMs,DENs]=butter(2,w0,'s'); %Butterworth order 2.
    [NUMdp,DENdp] = bilinear(NUMs,DENs,Fs,f0) ;%with prewarping
    
    acc = filtfilt(NUMdp,DENdp,acc.').';

   % acclong = filtfilt(NUMdp,DENdp,acclong);
    %% Filtering
    vx_his = [0];
    vy_his = [0];
    vx_gps = [0];
    vy_gps = [0];
    posx = [0];
    posy = [0];
%     px = [0];
%     py = [0];
    P = eye(2);
    Q = 0;
    R = 0.1;
    for i = 1:size(acc,2);
        [V, P] = prediction([vx_his(i);vy_his(i)], [acc(1,i);-acc(2,i)], P, Q); % prediction
        
        if (mod(i,20)==0)
            [V, P, K] = measUpdate(V, P, R, gps(i/20)); % update, rate 20 iteration 
            vx_gps(i/20+1) = V(1);
            vy_gps(i/20+1) = V(2);
        end
        vx_his(i+1) = V(1);
        vy_his(i+1) = V(2);
        posx(i+1) = posx(i) + vx_his(i)*dt;
        posy(i+1) = posy(i) + vy_his(i)*dt;
        
    end
    %% Plotting
    figure(1);
    subplot(2,1,1);
    plot(t,vx_his)
    title('Vx');
    xlabel('time (sec)'); ylabel('Velx (m/s)');
    subplot(2,1,2);
    plot(t,vy_his)
    title('Vy');
    xlabel('time (sec)'); ylabel('Vely (m/s)');
    
    figure(2);
    subplot(2,1,1);
    plot(t,posx)
    title('Posx');
    xlabel('time (sec)'); ylabel('Posx (m)');
    subplot(2,1,2);
    plot(t,posy)
    title('Posy');
    xlabel('time (sec)'); ylabel('Posy (m)');
    
    figure(3);
    subplot(2,1,1);
    plot(vx_gps)
    title('Vx');
    xlabel('time (sec)'); ylabel('Velx (m/s)');
    subplot(2,1,2);
    plot(vy_gps)
    title('Vy');
    xlabel('time (sec)'); ylabel('Vely (m/s)');

    
    
end