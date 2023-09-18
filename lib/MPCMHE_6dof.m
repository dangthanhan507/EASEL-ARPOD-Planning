classdef MPCMHE_6dof < MPCMHE
    properties
        % uMat = [1, -1,  0,  0,  0,  0;
        %         0,  0,  1, -1,  0,  0;
        %         0,  0,  0,  0,  1, -1];
        uMat = [1,  0,  0, -1,  0,  0;
                0,  1,  0,  0, -1,  0;
                0,  0,  1,  0,  0, -1];

        TD
        D0
        official
    end
    methods 
        %MATLAB subclassing is just like inheritance
        %all methods are inherited from MPCMHE
        function obj = init(obj, disturbType, backwardT, forwardT, state_dim, meas_dim, control_dim, att_dim, official)
            if official == true
                [obj.Tx0, obj.Tx, obj.Td, obj.TD, obj.Tuback, obj.Tyback, obj.Tvback, obj.Tuforward] = MPCMHE_6dofutils.createOfficialTenscalcVariables(disturbType, state_dim, control_dim, meas_dim, backwardT, forwardT);
            else
                [obj.Tx0, obj.Tx, obj.Td, obj.Tuback, obj.Tyback, obj.Tvback, obj.Tuforward] = MPCMHE_Tcalcutils.createTenscalcVariables(disturbType, state_dim, control_dim, meas_dim, backwardT, forwardT);
            end
            
            [obj.attTx0, obj.attTx, obj.attTd, obj.attTuback, obj.attTyback, obj.attTvback, obj.attTuforward] = MPCMHE_Tcalcutils.createTenscalcVariablesAttitude(0,att_dim, att_dim/2, backwardT, forwardT);

            obj.backwardT = backwardT;
            obj.forwardT = forwardT;
            obj.disturbType = disturbType;
            obj.official = official;

            obj.D0 = eye(control_dim);
        end
        function obj = setupLTIUncoupled(obj, A, B, C, mpcQ, mpcR, mheQ, mheR)
            obj.dynamicConstraints = MPCMHE_6dofutils.createFixedBodyLTIDynamics(A,B,obj.Tx0, obj.Tx, obj.Tuback, obj.Tuforward, obj.Td, obj.disturbType);
            obj.sensorConstraints  = MPCMHE_Tcalcutils.createLinearSensorConstraint(C, obj.Tx, obj.Tvback, obj.Tyback, obj.backwardT);
            obj.costFunction = MPCMHE_Tcalcutils.objectiveMPCMHE(mpcQ,mpcR,mheQ,mheR,obj.Tx, obj.Tuforward, obj.Td, obj.Tvback);
            obj.opt_properties = MPCMHE_6dofutils.setupOptimizationCells(obj.costFunction,obj.dynamicConstraints, obj.sensorConstraints,...
                                                        obj.dMax, obj.uMax, obj.vMax,...
                                                        obj.Tx0, obj.Tx, obj.Td, obj.TD, obj.Tuback, obj.Tyback, obj.Tvback,...
                                                        obj.Tuforward);
        end
        function obj = setupCoupledUnconstrainedNonlinear(obj, A, B, mpcQ, mpcR, mheQ, mheR)
            obj.dynamicConstraints = MPCMHE_6dofutils.createOfficialDynamics(A,B,obj.Tx0, obj.Tx, obj.attTx0, obj.attTx, obj.Tuback, obj.Tuforward, obj.TD, obj.Td, obj.disturbType);
            obj.costFunction       = MPCMHE_6dofutils.objectiveUnconstrainedMPCMHENonlinearMeas(A, B, mpcQ, mpcR, mheQ, mheR, obj.Tx, obj.Tuback, obj.Tuforward, obj.Tyback, obj.Td, obj.TD, obj.Tvback, obj.disturbType);
            obj.opt_properties = MPCMHE_6dofutils.setupUnconstrainedOptimizationCells(obj.costFunction,obj.dynamicConstraints,...
                                                        obj.dMax, obj.uMax, obj.vMax,...
                                                        obj.Tx0, obj.Tx, obj.Td, obj.TD, obj.Tuback, obj.Tyback, obj.Tvback,...
                                                        obj.Tuforward);
        end
        function obj = setupUncoupledUnconstrainedNonlinear(obj, A, B, mpcQ, mpcR, mheQ, mheR)
            obj.dynamicConstraints = MPCMHE_6dofutils.createFixedBodyLTIDynamics(A,B,obj.Tx0, obj.Tx, obj.Tuback, obj.Tuforward, obj.Td, obj.disturbType);
            obj.costFunction       = MPCMHE_6dofutils.objectiveUnconstrainedMPCMHENonlinearMeas(A, B, mpcQ, mpcR, mheQ, mheR, obj.Tx, obj.Tuback, obj.Tuforward, obj.Tyback, obj.Td, obj.TD, obj.Tvback, obj.disturbType);
            obj.opt_properties = MPCMHE_6dofutils.setupUnconstrainedOptimizationCells(obj.costFunction,obj.dynamicConstraints,...
                                                        obj.dMax, obj.uMax, obj.vMax,...
                                                        obj.Tx0, obj.Tx, obj.Td, obj.TD, obj.Tuback, obj.Tyback, obj.Tvback,...
                                                        obj.Tuforward);
        end
        function obj = setupLTICoupled(obj, A, B, C, mpcQ, mpcR, mheQ, mheR)
            obj.dynamicConstraints = MPCMHE_6dofutils.createCoupledLTIDynamics(A,B,obj.Tx0, obj.Tx, obj.attTx0, obj.attTx, obj.Tuback, obj.Tuforward, obj.Td, obj.disturbType);
            obj.sensorConstraints =  MPCMHE_Tcalcutils.createLinearSensorConstraint(C, obj.Tx, obj.Tvback, obj.Tyback, obj.backwardT);
            obj.costFunction = MPCMHE_Tcalcutils.objectiveMPCMHE(mpcQ,mpcR,mheQ,mheR,obj.Tx, obj.Tuforward, obj.Td, obj.Tvback);
            obj.opt_properties = MPCMHE_6dofutils.setupOptimizationCells(obj.costFunction,obj.dynamicConstraints, obj.sensorConstraints,...
                                                        obj.dMax, obj.uMax, obj.vMax,...
                                                        obj.Tx0, obj.Tx, obj.Td, obj.TD, obj.Tuback, obj.Tyback, obj.Tvback,...
                                                        obj.Tuforward);
        end
        function obj = optimize(obj)
            %mutating void functions
            MPCMHE_Tcalcutils.setupOptimizationVarsTrans(obj.opt, obj.window);
            MPCMHE_Tcalcutils.setupOptimizationVarsAtt(obj.opt, obj.att_window);
            
            
            maxIter = 10;
            if obj.official == true
                setV_D(obj.opt, obj.D0);
                [obj.window, obj.att_window, D1] = MPCMHE_6dofutils.mpcmhe_solve(obj.window, obj.att_window, obj.opt, 1, maxIter, -1, 1);
                obj.D0 = D1;
            else
                [obj.window, obj.att_window] = MPCMHE_Tcalcutils.mpcmhe_solve(obj.window, obj.att_window, obj.opt, 1, maxIter, -1, 1);
            end

        end
    end
end