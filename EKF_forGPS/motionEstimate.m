function posBellat = motionEstimate(accFile, gpsFile)
    % accFile = fullfile('Input','accelerometer_x5m.csv'); gpsFile = fullfile('Input','gps_x5m.csv');
    % posBel = motionEstimate(accFile,gpsFile);
    % 
    % Input
    %   accFile: path of the accelerometer file
    %   gpsFile: path of the gps file
    % Output
    %   posBel: array of the displacement
    %   plot: posBel plot
    
    %% data import
    [biasx, biasy, biasz,gpsbias, station] = biascalc('accelerometer_still.csv','gps_still.csv');
    acc = csvread('accelerometer.csv',1);
    gps = csvread('gps.csv',1);
    [a,b,c,d] = velfilt('accelerometer.csv', 'gps.csv');
    acclat =  acc(:,2)*9.81/station - biasx; % x-dir accelerometer
    acclong = (acc(:,4)*9.81/station) - biasz; % z-dir accelerometer
    gpsSpeed = vertcat(0,gps(:,6));%-gpsbias; % speed
%     gpsSpeed = gpsSpeed.*sin(theta);
    accT = 0.05; % accelerometer sample period
    gpsT = 1; % gps sample period

    %% low pass filter
%     [bl,al] = butter(12,0.1,'low');
    f0=0.25; %cut-off frequency
    w0=2*pi*f0;
    N=512;
    Fs=20; % sampling frequency
    [NUMs,DENs]=butter(2,w0,'s'); %Butterworth order 2.
    [NUMdp,DENdp] = bilinear(NUMs,DENs,Fs,f0) ;%with prewarping
    
    acclat = filtfilt(NUMdp,DENdp,acclat);
%     g = 0.9 * g + 0.1 * v
    acclong = filtfilt(NUMdp,DENdp,acclong);
        figure(1);
        subplot(3,1,1);
        plot(acclat);
        title('lat acc');
        subplot(3,1,2);
        plot(acclong);
        title('long acc');
        subplot(3,1,3);
        plot(gpsSpeed);
        title('gps speed');
    
%     theta = vertcat(0,atan2(acclat,acclong));

    %% main
    posBellat = 0; % position believe
    velBellat = 0; % velocity believe
    velBellong = 0;
    posBellong = 0;

    accVar = 1.0; % acceleromter variance
    GPSVar = 0.1; % GPS variance
    for ii = 1:size(acc,1)
        velBellat = vertcat(velBellat, velBellat(end) + acclat(ii)*accT);
        velBellong = vertcat(velBellong, velBellong(end) + acclong(ii)*accT);
        theta = atan2(acclat(ii),acclong(ii))
         if(mod(ii,20) == 0)
            gpsvel = c(ii/20);
            gpsvely = d(ii/20);%*sin(theta);
            combinedSpeedlat = (velBellat(end)/accVar + gpsvely/GPSVar)/(1/accVar + 1/GPSVar); % combine the gaussian 
            velBellat(end) = combinedSpeedlat;
            combinedSpeedlong = (velBellong(end)/accVar + gpsvel/GPSVar)/(1/accVar + 1/GPSVar); % combine the gaussian 
            velBellong(end) = combinedSpeedlong;
         end
        posBellat= vertcat(posBellat, posBellat(end)+velBellat(end)*accT);
        posBellong = vertcat(posBellong, posBellong(end)+velBellong(end)*accT);
    end
    figure(2);
    subplot(3,1,1);
    t = 0:accT:accT*size(acc,1);

    plot(t,posBellat)
    title('posBellat predict');
    xlabel('time (sec)'); ylabel('displacement (m)');
    subplot(3,1,2);
    plot(t,posBellong)
    title('posBellong predict');
    xlabel('time (sec)'); ylabel('displacement (m)');
    subplot(3,1,3);
    plot(t,theta)
    title('theta');
    xlabel('time (sec)'); ylabel('thetat (rad)');
    
end
  