classdef MPCMHE
    properties
        %Tenscalc variables for translational dynamics
        Tx0
        Tx
        Td
        Tuback
        Tyback
        Tvback
        Tuforward

        %Tenscalc variables for attitudinal dynamics
        attTx0
        attTx
        attTd
        attTuback
        attTyback
        attTvback
        attTuforward

        %Tenscalc constraint booleans
        dynamicConstraints
        sensorConstraints

        attdynamicConstraints
        attsensorConstraints

        %Tenscalc cost function
        costFunction
        attcostFunction

        %Variable Limits
        dMax
        vMax
        uMax

        %Optimizer wrapper
        opt_properties

        %parameters of mpcmhe
        backwardT
        forwardT

        %choice of mpcmhe design
        disturbType

        %opt
        opt

        %window data structures
        window
        att_window
    end

    methods
        function obj = init(obj, disturbType, backwardT, forwardT, state_dim, meas_dim, control_dim, att_dim)
            [obj.Tx0, obj.Tx, obj.Td, obj.Tuback, obj.Tyback, obj.Tvback, obj.Tuforward] = MPCMHE_Tcalcutils.createTenscalcVariables(disturbType, state_dim, control_dim, meas_dim, backwardT, forwardT);
            [obj.attTx0, obj.attTx, obj.attTd, obj.attTuback, obj.attTyback, obj.attTvback, obj.attTuforward] = MPCMHE_Tcalcutils.createTenscalcVariablesAttitude(0,att_dim, att_dim/2, backwardT, forwardT);

            obj.backwardT = backwardT;
            obj.forwardT = forwardT;
            obj.disturbType = disturbType;
        end

        %%%% TENSCALC SETUP OPTIMIZATION CODE BEGIN %%%%
        function obj = setupUnconstrainedFixedBody(obj, A, B, C, mpcQ, mpcR, mheQ, mheR)
            %runs objective and constraint at once
            obj.dynamicConstraints = MPCMHE_Tcalcutils.createLTIDynamicConstraintsMPC(A,B,obj.Tx, obj.Tuforward, obj.Td, obj.disturbType);
            obj.costFunction       = MPCMHE_Tcalcutils.objectiveUnconstrainedMPCMHE(A,B,C, mpcQ, mpcR, mheQ, mheR, obj.Tx, obj.Tuback, obj.Tuforward, obj.Tyback, obj.Td, obj.disturbType);
        end
        function obj = setupLTILinearFixedBody(obj, A, B, C)
            obj.dynamicConstraints = MPCMHE_Tcalcutils.createLTIDynamicConstraints(A,B,obj.Tx0, obj.Tx, obj.Tuback, obj.Tuforward, obj.Td, obj.disturbType);
            obj.sensorConstraints = MPCMHE_Tcalcutils.createLinearSensorConstraint(C, obj.Tx, obj.Tvback, obj.Tyback, obj.backwardT);
        end
        function obj = setupSimpleObjective(obj, mpcQ, mpcR, mheQ, mheR)
            obj.costFunction = MPCMHE_Tcalcutils.objectiveMPCMHE(mpcQ,mpcR,mheQ,mheR,obj.Tx, obj.Tuforward, obj.Td, obj.Tvback);
        end
        function obj = setupVariableLimits(obj, dMax, uMax, vMax)
            obj.dMax = dMax;
            obj.uMax = uMax;
            obj.vMax = vMax;
        end
        function obj = createOptimization(obj)
            obj.opt_properties = MPCMHE_Tcalcutils.setupOptimizationCells(obj.costFunction,obj.dynamicConstraints, obj.sensorConstraints,...
                                                        obj.dMax, obj.uMax, obj.vMax,...
                                                        obj.Tx0, obj.Tx, obj.Td, obj.Tuback, obj.Tyback, obj.Tvback,...
                                                        obj.Tuforward);
        end
        function obj = setupAttitudeConstraints(obj, A, B, C)
            obj.attdynamicConstraints = MPCMHE_Tcalcutils.createLTIDynamicConstraints(A,B,obj.attTx0, obj.attTx, obj.attTuback, obj.attTuforward, obj.attTd, 0);
            obj.attsensorConstraints = MPCMHE_Tcalcutils.createLinearSensorConstraint(C, obj.attTx, obj.attTvback, obj.attTyback, obj.backwardT);
        end
        function obj = setupAttitudeCost(obj, mpcQ, mpcR, mheQ, mheR)
            obj.attcostFunction = MPCMHE_Tcalcutils.objectiveMPCMHE(mpcQ,mpcR,mheQ,mheR,obj.attTx, obj.attTuforward, obj.attTd, obj.attTvback);
        end
        function obj = addAttitudeOptimization(obj)
            other_properties = MPCMHE_Tcalcutils.setupOptimizationCells(obj.attcostFunction, obj.attdynamicConstraints, obj.attsensorConstraints,...
                                                                        obj.dMax, obj.uMax, obj.vMax,...
                                                                        obj.attTx0, obj.attTx, obj.attTd, obj.attTuback, obj.attTyback, obj.attTvback,...
                                                                        obj.attTuforward);
            obj.opt_properties = obj.opt_properties.concat(other_properties);
        end
        function obj = setupOptimizationVar(obj)
            obj.opt = obj.opt_properties.setupOptimization();
        end
        %%%%  TENSCALC SETUP OPTIMIZATION CODE END  %%%%

        %%%%  Window Code BEGIN %%%%
        function obj = setupWindows(obj, x0, window_states, window_measurements, window_controls)
            obj.window = MPCMHE_windows;
            obj.window = obj.window.init(x0, window_states, window_measurements, window_controls, obj.forwardT, obj.disturbType);
        end
        function obj = setupAttWindows(obj, x0, window_states, window_measurements, window_controls)
            obj.att_window = MPCMHE_windows;
            obj.att_window = obj.att_window.init(x0, window_states, window_measurements, window_controls, obj.forwardT, 0);
        end
        function obj = shift(obj, meas, control, att_meas, att_control)
            obj.window = obj.window.shiftWindows(meas,control);
            obj.att_window = obj.att_window.shiftWindows(att_meas, att_control);
        end
        %%%%  Window Code END   %%%%

        function obj = optimize(obj)
            %mutating void functions
            MPCMHE_Tcalcutils.setupOptimizationVarsTrans(obj.opt, obj.window);
            MPCMHE_Tcalcutils.setupOptimizationVarsAtt(obj.opt, obj.att_window);

            [obj.window, obj.att_window] = MPCMHE_Tcalcutils.mpcmhe_solve(obj.window, obj.att_window, obj.opt, 1, 1000, -1, 1);
        end
        function [mheXs, mheDs, mheVs, mpcXs, mpcUs] = getOptimizeResult(obj)
            mheXs = obj.window.window_mhestates;

            mheDs = obj.window.window_controlDisturbances;
            
            mheVs = obj.window.window_measError;
            
            mpcUs = obj.window.window_mpccontrols;
            
            mpcXs = obj.window.window_mpcstates;
            
        end
        function [mheXs, mheDs, mheVs, mpcXs, mpcUs] = getOptimizeResultAttitude(obj)
            mheXs = obj.att_window.window_mhestates;
            mheDs = obj.att_window.window_controlDisturbances;
            mheVs = obj.att_window.window_measError;
            mpcUs = obj.att_window.window_mpccontrols;
            mpcXs = obj.att_window.window_mpcstates;
        end
        function u = getMPCControl(obj)
            u = obj.window.window_mpccontrols(:,1);
        end
        function u = getMPCAttControl(obj)
            u = obj.att_window.window_mpccontrols(:,1);
        end
        function x = getMHEState(obj)
            x = obj.window.window_mhestates(:,end);
        end
        function x = getMHEAttState(obj)
            x = obj.att_window.window_mhestates(:,end);
        end
        
    end
end