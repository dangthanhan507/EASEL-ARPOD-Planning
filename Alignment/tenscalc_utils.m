classdef tenscalc_utils
    methods (Static)
        function output = FullRotation(theta, us)
            %{
                theta = [roll_1, roll_2, ..., roll_n]
                        [pitch_1, pitch_2, ..., pitch_n]
                        [yaw_1, yaw_2, ..., yaw_n]
            %}

            %fully expanded rotation matrix but we vectorize this
            R_11 = cos(theta(3,:)).*cos(theta(2,:));
            R_12 = -sin(theta(3,:)).*cos(theta(1,:)) + cos(theta(3,:)).*sin(theta(2,:)).*sin(theta(1,:));
            R_13 = sin(theta(3,:)).*sin(theta(1,:)) + cos(theta(3,:)).*sin(theta(2,:)).*cos(theta(1,:));

            R_21 = sin(theta(3,:)).*cos(theta(2,:));
            R_22 = cos(theta(3,:)).*cos(theta(1,:)) + sin(theta(3,:)).*sin(theta(2,:)).*sin(theta(1,:));
            R_23 = -cos(theta(3,:)).*sin(theta(1,:)) + sin(theta(3,:)).*sin(theta(2,:)).*cos(theta(1,:));

            R_31 = -sin(theta(2,:));
            R_32 = cos(theta(2,:)).*sin(theta(1,:));
            R_33 = cos(theta(2,:)).*cos(theta(1,:));
            
            %size(us) = 3 x N
            output = [sum(R_11.*us(1,:) + R_12.*us(2,:) + R_13.*us(3,:), 1),... % N size
                      sum(R_21.*us(1,:) + R_22.*us(2,:) + R_23.*us(3,:), 1),...
                      sum(R_31.*us(1,:) + R_32.*us(2,:) + R_33.*us(3,:), 1)].';
        end
        function z_t = SensingVectorized(states)
            [statedim, N] = size(states);
            xs = states(1,:);
            ys = states(2,:);
            zs = states(3,:);

            % rho = norm(states(1:3,:),2); %1xN
            rho = sqrt(power(states(1,:),2) + power(states(2,:),2) + power(states(3,:),2));
            
            z_norm = (zs ./ rho);
            e1 = atan(ys ./ xs);
            e2 = atan( z_norm ./ sqrt((1+z_norm).*(1-z_norm)));
            z_t = Tzeros(3,N);
            z_t(1,:) = e1;
            z_t(2,:) = e2;
            z_t(3,:) = rho;
        end

        function dynamicConstraints = mheDynamicsTranslational(A, B, Tx0, Tx, Tatt0, Tatt, uk, d)
            xk = [Tx0,Tx(:,1:end-1)];
            attk = [Tatt0, Tatt(:,1:end-1)];

            u = uk + d;
            u = tenscalc_utils.FullRotation(attk(1:3,:), u(1:3,:) - u(4:6,:));
            dynamicConstraints = (Tx == A*xk + B*u);
        end
        function Jcost = mheObjectiveTranslational(mheQ, mheR, Tx, meas, disturbance)
            g_of_x = tenscalc_utils.SensingVectorized(Tx(:,1:backwardHorizon));
            Jcost = mheQ*norm2(disturbance) + mheR*norm2(meas - g_of_x);
        end

        function dynamicConstraints = mheDynamicsAttitude(A, B, Tatt0, Tatt, uatt)
            attk = [Tatt0, Tatt(:,1:end-1)];
            dynamicConstraints = (Tatt == A*attk + B*uatt);
        end
        function Jcost = mheObjectiveAttitude(mheQ, mheR, Tatt, meas)
            Jcost = mheR*norm2(meas - Tatt);
        end

        %%%%%%%% TRANSLATION MPC SETUP %%%%%%%%
        function dynamicConstraints = mpcDynamicsTranslational(A, B, Tx0, Tx, Tatt0, Tatt, uforward)
            xk = [Tx0, Tx(:,1:end-1)];  
            attk = [Tatt0, Tatt(:,1:end-1)];

            u = tenscalc_utils.FullRotation(attk(1:3,:), uforward(1:3,:) - uforward(4:6,:));
            dynamicConstraints = (Tx == A*xk + B*u);
        end
        function Jcost = mpcObjectiveTranslational(mpcQ, mpcR, Tx, uforward)
            Jcost = mpcQ*norm2(Tx) + mpcR*norm2(uforward);
        end

        %%%%%%%% ATTITUDINAL MPC SETUP %%%%%%%%
        function dynamicConstraints = mpcDynamicsAttitude(A, B, Tatt0, Tatt, uAtt)
            attk = [Tatt0, Tatt(:,1:end-1)];

            dynamicConstraints = (Tatt == A*attk + B*uAtt);
        end
        function Jcost = mpcObjectiveAttitude(mpcQ, mpcR, Tatt, uAtt)
            Jcost = mpcQ*norm2(Tatt) + mpcR*norm2(uAtt);
        end


        %%%%%%%% TRANSLATION MPCMHE SETUP %%%%%%%%
        function dynamicConstraints = mpcmheDynamicsTranslational(A, B, Tx0, Tx, Tatt0, Tatt, uback, uforward, D, disturbance)
            xk = [Tx0, Tx(:,1:end-1)];  
            attk = [Tatt0, Tatt(:,1:end-1)];
            uk = [uback, uforward];

            % uk = D*uk + disturbance;
            uk = D*uk;
            u = tenscalc_utils.FullRotation(attk(1:3,:), uk(1:3,:) - uk(4:6,:));
            % u = uk(1:3,:) - uk(4:6,:);
            dynamicConstraints = (Tx == A*xk + B*u);
        end

        function Jcost = mpcmheObjectiveTranslational(mpcQ, mpcR, mheQ, mheR, Tx, uForward, meas, disturbance, TD)
            [meas_dim, backwardHorizon] = size(meas);
            [control_dim, forwardHorizon] = size(uForward);

            g_of_x = tenscalc_utils.SensingVectorized(Tx(:,1:backwardHorizon));

            Jcost = mpcQ*norm2( Tx(:,backwardHorizon+1:forwardHorizon) ) + mpcR*norm2(uForward);
            Jcost = Jcost - mheQ*norm2(TD) - mheQ*norm2(disturbance) - mheR*norm2(meas - g_of_x );
            % Jcost = Jcost - mheQ*norm2(TD) - mheQ*norm2(disturbance);
        end


        %%%%%%%% ATTITUDINAL MPCMHE SETUP %%%%%%%%
        function Jcost = mpcmheObjectiveAttitude(mpcQ, mpcR, mheQ, mheR, Tx, uForward, sensor_disturbance)
            [meas_dim, backwardHorizon] = size(sensor_disturbance);
            [control_dim, forwardHorizon] = size(uForward);

            Jcost = mpcQ*norm2( Tx(:,backwardHorizon+1:forwardHorizon) ) + mpcR*norm2(uForward);
            Jcost = Jcost - mheR*norm2(sensor_disturbance); % only sensor noise established
        end
        function dynamicConstraints = mpcmheDynamicsAttitude(Aatt, Batt, Tx0, Tx, uback, uforward)
            xk = [Tx0, Tx(:,1:end-1)];
            uk = [uback, uforward];

            dynamicConstraints = (Tx == Aatt*xk + Batt*uk);
        end
        function sensorConstraints = mpcmheSensorConstraintsAttitude(Tx, sensorDisturb, measurements)
            [measdim, backwardHorizon] = size(measurements);
            C = eye(6);

            sensorConstraints = (measurements == C*Tx(:,1:backwardHorizon) + sensorDisturb);
        end
    end
end