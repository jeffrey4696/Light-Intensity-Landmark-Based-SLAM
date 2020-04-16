function [v, P, K] = measUpdate(v, P, R, z)
            [H, HT, z_hat] = mJacobian(v); % measurment Jacobian
            Qz = H*P*HT + R; % measurement covariance, dim: 3*3
            if det(Qz) < 1e-6
                QzInv = pinv(Qz);
            else
                QzInv = inv(Qz);
            end
            
            % Kalman gain
            K = P*HT*QzInv; 
            
            % update states and covariance
            v = v + K*(z-z_hat);
            P = (eye(2) - K*H)*P;
            
end