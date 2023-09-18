classdef MPCMHE_6dofutils
    methods (Static)
        function [Tx0, Tx, Td, TD, Tuback, Tyback, Tvback, Tuforward] = createOfficialTenscalcVariables(disturbType, statedim, controldim, meas_dim, backHorizon, forwardHorizon)
            Tvariable x0a   [statedim,1];
            Tvariable x     [statedim, backHorizon + forwardHorizon];
            Tvariable uback [controldim, backHorizon];
            Tvariable vback [meas_dim, backHorizon];
            Tvariable ypast [meas_dim, backHorizon];
            Tvariable u     [controldim, forwardHorizon];

            Tvariable d [controldim, backHorizon + forwardHorizon];
            Tvariable D [controldim, controldim];

            Tx0 = x0a;
            Tx = x;
            Td = d;
            Tuback = uback;
            Tyback = ypast;
            Tvback = vback;
            Tuforward = u;
            TD = D;
        end

        function dynamicConstraints = createFixedBodyLTIDynamics(A, B, Tx0, Tx, uback, uforward, disturbance, disturbType)

            xk = [Tx0, Tx(:,1:end-1)];
            uk = [uback, uforward];

            [control_dim, n] = size(uk);
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

            % 6dof vs normal spacecraft
            %{
                u is 6x1
                [ 1 -1  0  0  0  0] [ux+]
                [ 0  0  1 -1  0  0] [ux-]    [ux]
                [ 0  0  0  0  1 -1] [uy+]  = [uy]
                                    [uy-]    [uz]
                                    [uz+]
                                    [uz-]
            %}
            dynamicConstraints = (Tx == A*xk + B*(uk(1:3,:) - uk(4:6,:)));
        end
        function opt_properties = setupOptimizationCells(costFunction, dynamicConstraints, sensorConstraints, dMax, uMax, vMax, Tx0, Tx, Td, TD, Tuback, Tyback, Tvback, Tuforward)
            %only difference is we set uMin to zero
            opt_properties = MPCMHE_properties;
            [control_dim, forwardT] = size(Tuforward);
            [meas_dim, backwardT] = size(Tyback);
            opt_properties = opt_properties.init(costFunction,  {Tuforward}, ...
                                                                {Td, Tvback,Tx0}, ...
                                                                {Tx},...
                                                                {Tuforward <= uMax*Tones(size(Tuforward)), Tuforward >= -1e-5*Tones(size(Tuforward))},...
                                                                {Td<=dMax*Tones(size(Td)), Td>=-dMax*Tones(size(Td)),...
                                                                    Tvback<=vMax*Tones(size(Tvback)), Tvback>=-vMax*Tones(size(Tvback)),...
                                                                    sensorConstraints},...
                                                                {dynamicConstraints},...
                                                                {Tuforward,Td,Tx0,Tx,Tvback},...
                                                                {Tuback,Tyback});
            opt_properties = opt_properties.addTD(TD);
        end
        function opt_properties = setupUnconstrainedOptimizationCells(costFunction, dynamicConstraints, dMax, uMax, vMax, Tx0, Tx, Td, TD, Tuback, Tyback, Tvback, Tuforward)
            opt_properties = MPCMHE_properties;
            opt_properties = opt_properties.init(costFunction,  {Tuforward}, ...
                                                                {Td, Tvback,Tx0}, ...
                                                                {Tx},...
                                                                    {Tuforward <= uMax*Tones(size(Tuforward)), Tuforward >= -1e-5*Tones(size(Tuforward))},...
                                                                {Td<=dMax*Tones(size(Td)), Td>=-dMax*Tones(size(Td)),...
                                                                    Tvback<=vMax*Tones(size(Tvback)), Tvback>=-vMax*Tones(size(Tvback))},...
                                                                {dynamicConstraints},...
                                                                {Tuforward,Td,Tx0,Tx,Tvback},...
                                                                {Tuback,Tyback});
            opt_properties = opt_properties.addTD(TD);
        end        

        %======== Adding Rotation =============
        function rotmat = Rotation(roll, pitch, yaw)
            Rx = [1 0 0; 0 cos(roll) -sin(roll); 0 sin(roll) cos(roll)];
            Ry = [cos(pitch) 0 sin(pitch); 0 1 0; -sin(pitch) 0 cos(pitch)];
            Rz = [cos(yaw) -sin(yaw) 0; sin(yaw) cos(yaw) 0; 0 0 1];
            rotmat = Rz*Ry*Rx;
        end
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
        function dynamicConstraints = createOfficialDynamics(A, B, Tx0, Tx, Tatt0, Tatt, uback, uforward, D, disturbance, disturbType)
            xk = [Tx0, Tx(:,1:end-1)];  
            attk = [Tatt0, Tatt(:,1:end-1)];
            uk = [uback, uforward];

            [control_dim, n] = size(uk);

            uk = D*uk + disturbance;
            u = MPCMHE_6dofutils.FullRotation(attk(1:3,:),uk(1:3,:) - uk(4:6,:));
            dynamicConstraints = (Tx == A*xk + B*u);
        end
        function dynamicConstraints = createCoupledLTIDynamics(A, B, Tx0, Tx, Tatt0, Tatt, uback, uforward, disturbance, disturbType)
            xk = [Tx0, Tx(:,1:end-1)];  
            attk = [Tatt0, Tatt(:,1:end-1)];
            uk = [uback, uforward];
            
            % 6dof vs normal spacecraft
            %{
                u is 6x1
                [ 1 -1  0  0  0  0] [ux+]
                [ 0  0  1 -1  0  0] [ux-]    [ux]
                [ 0  0  0  0  1 -1] [uy+]  = [uy]
                                    [uy-]    [uz]
                                    [uz+]
                                    [uz-]
            %}

            [control_dim, n] = size(uk);
            u = Tzeros(3, n);

            for i = 1:n
                % R = MPCMHE_6dofutils.Rotation(attk(1,i), attk(2,i), attk(3,i));
                if disturbType == 0
                elseif disturbType == 1
                    % uk(:,i) = (uk(:,i) + disturbance(:,i));
                elseif disturbType == 2
                    uk(:,i) = (uk(:,i) + disturbance);

                elseif disturbType == 3
                    uk(:,i) = disturbance*uk(:,i);
                else
                    dMatrix = disturbance(1:control_dim, 1:control_dim);
                    dAdd    = disturbance(1:control_dim, end);
                    uk(:,i) = dMatrix*uk(:,i) + dAdd;
                end
            end
            uk = uk + disturbance;
            u = MPCMHE_6dofutils.FullRotation(attk(1:3,:),uk(1:3,:) - uk(4:6,:));
            dynamicConstraints = (Tx == A*xk + B*u);

        end
        function Jcost = objectiveUnconstrainedMPCMHENonlinearMeas(A, B, mpcQ, mpcR, mheQ, mheR, Tx, uBack, uForward, meas, disturbance, TD, sensor_disturbance, disturbType)
            [meas_dim, backwardHorizon] = size(sensor_disturbance);
            [control_dim, forwardHorizon] = size(uForward);
            
            Jcost = mpcQ*norm2( Tx(:,backwardHorizon+1:forwardHorizon) ) + mpcR*norm2(uForward);
            Jcost = Jcost - mheQ*norm2(TD) - mheQ*norm2(disturbance) - mheR*norm2(meas - MPCMHE_Tcalcutils.applyNonlinearSensor(Tx(:,1:backwardHorizon)) );
        end
        function [window,att_window, D0] = mpcmhe_solve(window, att_window, opt, mu0, maxIter, saveIter, use_attitude)
            [status,iter,time] = solve(opt,mu0,int32(maxIter),int32(saveIter));

            disp("Status:")
            disp(status)
            disp("Iterations:")
            disp(iter)
            
            [Jcost, mpcUs, mheDs, x0s, xs, vs, D0, attmpcUs, attmheDs, attx0s, attxs, attvs] = getOutputs(opt);
            [measdim,backwardT] = size(vs);
            attmheXs = attxs(:,1:backwardT);
            attmpcXs = attxs(:,backwardT+1:end);

            att_window.window_controlDisturbances = attmheDs;%TODO: FIX THIS
            att_window.window_measError = attvs;
            att_window.window_mhestates = attmheXs;
            att_window.window_mpcstates = attmpcXs;
            att_window.window_mpccontrols = attmpcUs; %TODO: FIX THIS

            mheXs = xs(:,1:backwardT);
            mpcXs = xs(:,backwardT+1:end);

            disp("Jcost:")
            disp(Jcost)

            window.window_controlDisturbances = mheDs;
            window.window_measError = vs;
            window.window_mhestates = mheXs;
            window.window_mpcstates = mpcXs;
            window.window_mpccontrols = mpcUs;

            disp("MPC Future")
            disp(window.window_mpcstates)
            disp("MPC control horizon")
            disp(window.window_mpccontrols)
        end
        %======================================
    end
end