function [v, P] = prediction(v, u, P, Q)
    % predict the states
        dt = 0.05; %time step for accelerometer
        v = v + u.*dt;
    % predict the covariance
        P = P + Q;         
end