classdef MPCMHE_Tcalcutils
    methods (Static)
        function [Tx0, Tx, Td, Tuback, Tyback, Tvback, Tuforward] = createTenscalcVariables(disturbType, statedim, controldim, meas_dim, backHorizon, forwardHorizon)
            Tvariable x0a   [statedim,1];
            Tvariable x     [statedim, backHorizon + forwardHorizon];
            Tvariable uback [controldim, backHorizon];
            Tvariable vback [statedim, backHorizon];
            Tvariable ypast [meas_dim, backHorizon];
            Tvariable u     [controldim, forwardHorizon];

            if disturbType == 0
                Tvariable d     [controldim, 1];
            elseif disturbType == 1
                Tvariable d     [controldim, backHorizon + forwardHorizon];
            elseif disturbType == 2
                Tvariable d     [controldim, 1];
            elseif disturbType == 3
                Tvariable d     [controldim, controldim];
            else
                Tvariable d     [controldim, controldim+1]; %transform + additive concatenated
            end
            Tx0 = x0a;
            Tx = x;
            Td = d;
            Tuback = uback;
            Tyback = ypast;
            Tvback = vback;
            Tuforward = u;
        end
        
        function [Tx0, Tx, Td, Tuback, Tyback, Tvback, Tuforward] = createTenscalcVariablesAttitude(disturbType, statedim, controldim, backHorizon, forwardHorizon)
            %the Tvariable names MATTER SINCE THEY ARE GLOBAL
            Tvariable att0a   [statedim,1];
            Tvariable att     [statedim, backHorizon + forwardHorizon];
            Tvariable attuback [controldim, backHorizon];
            Tvariable attvback [statedim, backHorizon];
            Tvariable attypast [statedim, backHorizon];
            Tvariable attu     [controldim, forwardHorizon];

            if disturbType == 0
                Tvariable attd     [controldim, 1];
            elseif disturbType == 1
                Tvariable attd     [controldim, backHorizon + forwardHorizon];
            end
            Tx0 = att0a;
            Tx = att;
            Td = attd;
            Tuback = attuback;
            Tyback = attypast;
            Tvback = attvback;
            Tuforward = attu;
        end
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
        function dynamicConstraints = createLTIDynamicConstraintsMPC(A,B,Tx, uforward, disturbance, disturbType)
            [control_dim, forwardT] = size(uforward);
            [state_dim, backforT]   = size(Tx);
            backwardT = backforT - forwardT;
            xk = Tx(:,backwardT:end-1);
            uk = uforward;

            n = forwardT;
            if disturbType == 0
            elseif disturbType == 1
                %time-varying additive
                for i = 1:n
                    uk(:,i) = uk(:,i) + disturbance(:,i+backwardT-1);
                end
            elseif disturbType == 2
                %time-invariant additive
                for i = 1:n
                    uk(:,i) = uk(:,i) + disturbance;
                end
            elseif disturbType == 3
                %time-invariant matrix transform
                for i = 1:n
                    uk(:,i) = disturbance*uk(:,i+backwardT-1);
                end
            else
                %time-invariant matrix transform + additive
                dMatrix = disturbance(1:control_dim, 1:control_dim);
                dAdd    = disturbance(1:control_dim, end);
                for i = 1:n
                    uk(:,i) = dMatrix*uk(:,i) + dAdd;
                end
            end

            dynamicConstraints = (Tx(:,backwardT+1:end) == A*xk + B*uk);
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

            Jcost = mpcQ*norm2( Tx(:,backwardHorizon+1:forwardHorizon) ) + mpcR*norm2(control);
            Jcost = Jcost - mheQ*norm2(disturbance) - mheR*norm2(sensor_disturbance);
        end
        function Jcost = objectiveUnconstrainedMPCMHE(A, B, C, mpcQ, mpcR, mheQ, mheR, Tx, uBack, uForward, meas, disturbance, sensor_disturbance, disturbType)
            [meas_dim, backwardHorizon] = size(sensor_disturbance);
            [control_dim, forwardHorizon] = size(uForward);

            uk = uBack;
            n = backwardHorizon-1;
            if disturbType == 0
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
            xk = [Tx0, Tx(:,1:backwardHorizon-1)];
            Jcost = mpcQ*norm2( Tx(:,backwardHorizon+1:forwardHorizon) ) + mpcR*norm2(uForward);
            Jcost = Jcost - mheQ*norm2(Tx(:,1:backwardHorizon) - (A*xk + B*uk)) - mheR*norm2(meas - (C*Tx(:,1:backwardHorizon) + sensor_disturbance));
        end

        function opt_properties = setupOptimizationCells(costFunction, dynamicConstraints, sensorConstraints, dMax, uMax, vMax, Tx0, Tx, Td, Tuback, Tyback, Tvback, Tuforward)
            opt_properties = MPCMHE_properties;
            opt_properties = opt_properties.init(costFunction,  {Tuforward}, ...
                                                                {Td, Tvback,Tx0}, ...
                                                                {Tx},...
                                                                {Tuforward <= uMax*Tones(size(Tuforward)), Tuforward >= -uMax*Tones(size(Tuforward))},...
                                                                {Td<=dMax*Tones(size(Td)), Td>=-dMax*Tones(size(Td)),...
                                                                    Tvback<=vMax*Tones(size(Tvback)), Tvback>=-vMax*Tones(size(Tvback)),...
                                                                    sensorConstraints},...
                                                                {dynamicConstraints},...
                                                                {Tuforward,Td,Tx0,Tx,Tvback},...
                                                                {Tuback,Tyback});
        end
        function setupOptimizationVarsTrans(opt, window)

            %setting translation stuff
            setP_uback(opt, window.window_mhecontrols);
            setP_ypast(opt, window.window_measurements);
            setV_u(opt, window.window_mpccontrols);
            setV_d(opt, window.window_controlDisturbances);
            setV_x0a(opt, window.x0);
            setV_x(opt, [window.window_mhestates, window.window_mpcstates]);
            setV_vback(opt, window.window_measError);
        end
        function setupOptimizationVarsAtt(opt, window)
            setP_attuback(opt, window.window_mhecontrols);
            setP_attypast(opt, window.window_measurements);
            setV_attu(opt, window.window_mpccontrols);
            setV_attd(opt, window.window_controlDisturbances);
            setV_att0a(opt, window.x0);
            setV_att(opt, [window.window_mhestates, window.window_mpcstates]);
            setV_attvback(opt, window.window_measError);
        end
        function [window,att_window] = mpcmhe_solve(window, att_window, opt, mu0, maxIter, saveIter, use_attitude)
            [status,iter,time] = solve(opt,mu0,int32(maxIter),int32(saveIter));

            if use_attitude
                [Jcost, mpcUs, mheDs, x0s, xs, vs, attmpcUs, attmheDs, attx0s, attxs, attvs] = getOutputs(opt);
                [measdim,backwardT] = size(vs);
                attmheXs = attxs(:,1:backwardT);
                attmpcXs = attxs(:,backwardT+1:end);

                att_window.window_controlDisturbances = attmheDs;%TODO: FIX THIS
                att_window.window_measError = attvs;
                att_window.window_mhestates = attmheXs;
                att_window.window_mpcstates = attmpcXs;
                att_window.window_mpccontrols = attmpcUs; %TODO: FIX THIS

            else
                [Jcost, mpcUs, mheDs, x0s, xs, vs] = getOutputs(opt); %get results
                [measdim,backwardT] = size(vs);
            end
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
        
        %{
            arctan function for MHE
        %}
        function z_t = Sensing(state)
            x = state(1,:);
            y = state(2,:);
            z = state(3,:);
        
            norm = sqrt(x.*x+y.*y+z.*z);
            e1 = atan(y ./ x);
            z_norm = z ./ norm;
            e2 = atan( z_norm ./ sqrt((1+z_norm)*(1-z_norm)));
            e3 = norm;
            z_t = [e1;e2;e3];
        end
        function gX = applyNonlinearSensor(Tstates)
            [state_dim, num] = size(Tstates);
            gX = Tzeros(size(Tstates));
            for i = 1:num
                gX(:,i) = MPCMHE_Tcalcutils.Sensing(Tstates(:,i));
            end
            %return sensor version of everything
        end
        function measurementConstraints = createNonlinearSensorConstraint(Tx, sensorDisturb, measurements, backwardHorizon)
            measurementConstraints = (measurements == MPCMHE_Tcalcutils.applyNonlinearSensor(Tx(:,1:backwardHorizon)) + sensorDisturb);
        end
        function Jcost = objectiveUnconstrainedMPCMHENonlinearMeas(A, B, mpcQ, mpcR, mheQ, mheR, Tx, uBack, uForward, meas, disturbance, sensor_disturbance, disturbType)
            [meas_dim, backwardHorizon] = size(sensor_disturbance);
            [control_dim, forwardHorizon] = size(uForward);

            uk = uBack;
            n = backwardHorizon-1;
            if disturbType == 0
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
            xk = [Tx0, Tx(:,1:backwardHorizon-1)];
            Jcost = mpcQ*norm2( Tx(:,backwardHorizon+1:forwardHorizon) ) + mpcR*norm2(uForward);
            Jcost = Jcost - mheQ*norm2(Tx(:,1:backwardHorizon) - (A*xk + B*uk)) - mheR*norm2(meas - ( MPCMHE_Tcalcutils.applyNonlinearSensor(Tx(:,1:backwardHorizon)) + sensor_disturbance));
        end
    end
end