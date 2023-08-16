classdef MPCMHE_Tcalcutils
    methods (Static)
        function dynamicConstraints = createLTIDynamicConstraints(A, B, Tx0, Tx, uback, uforward, disturbance, disturbType)
            %{
                Description:
                ------------

                Parameters:
                -----------

                Returns:
                --------
            %}

            xk = [Tx0, Tx(:,1:end-1)];
            uk = [uback, uforward];

            [control_dim, n] = size(uk);

            %parse disturbance
            %NOTE: pass disturbType = 0, if you don't want to use disturbance
            if disturbType == 0
                %leave empty
            elseif disturbType == 1
                %time-varying additive
                for i = 1:n
                    uk(:,i) = uk(:,i) + disturbance(:,i);
                end
            elseif disturbType == 2
                %time-invariant additive
                for i = 1:n
                    uk(:,i) = uk(:,i) + disturbance;
                end
            elseif disturbType == 3
                %time-invariant matrix transform
                for i = 1:n
                    uk(:,i) = disturbance*uk(:,i);
                end
            else
                %time-invariant matrix transform + additive
                dMatrix = disturbance(1:control_dim, 1:control_dim);
                dAdd    = disturbance(1:control_dim, end);
                for i = 1:n
                    uk(:,i) = dMatrix*uk(:,i) + dAdd;
                end
            end

            dynamicConstraints = (Tx == A*xk + B*uk);
        end
        function measurementConstraints = createLinearSensorConstraint(C, Tx, sensorDisturb, measurements, backwardHorizon)
            %{
                Description:
                ------------

                Parameters:
                -----------

                Returns:
                --------
            %}
            measurementConstraints = (measurements == C*Tx(:,1:backwardHorizon) + sensorDisturb);
        end
        function Jcost = objectiveMPCMHE(mpcQ, mpcR, mheQ, mheR, Tx, control, disturbance, sensor_disturbance)
            %{
                Description:
                ------------

                Parameters:
                -----------

                Returns:
                --------
            %}
            [meas_dim, backwardHorizon] = size(sensor_disturbance);
            [control_dim, forwardHorizon] = size(control);

            Jcost = mpcQ*norm2(Tx(:,backwardHorizon+1:forwardHorizon)) + mpcR*norm2(control);
            Jcost = Jcost - mheQ*norm2(disturbance) - mheR*norm2(sensor_disturbance);
        end
    end
end