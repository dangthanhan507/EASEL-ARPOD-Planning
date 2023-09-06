classdef MPCMHE_6dof < MPCMHE
    properties
        uMat = [1, -1,  0,  0,  0,  0;
                0,  0,  1, -1,  0,  0;
                0,  0,  0,  0,  1, -1];
    end
    methods 
        %MATLAB subclassing is just like inheritance
        %all methods are inherited from MPCMHE
        function obj = setupLTIUncoupled(obj, A, B, C, mpcQ, mpcR, mheQ, mheR)
            obj.dynamicConstraints = MPCMHE_6dofutils.createFixedBodyLTIDynamics(A,B,obj.Tx0, obj.Tx, obj.Tuback, obj.Tuforward, obj.Td, obj.disturbType);
            obj.sensorConstraints  = MPCMHE_Tcalcutils.createLinearSensorConstraint(C, obj.Tx, obj.Tvback, obj.Tyback, obj.backwardT);
            obj.costFunction = MPCMHE_Tcalcutils.objectiveMPCMHE(mpcQ,mpcR,mheQ,mheR,obj.Tx, obj.Tuforward, obj.Td, obj.Tvback);
            obj.opt_properties = MPCMHE_6dofutils.setupOptimizationCells(obj.costFunction,obj.dynamicConstraints, obj.sensorConstraints,...
                                                        obj.dMax, obj.uMax, obj.vMax,...
                                                        obj.Tx0, obj.Tx, obj.Td, obj.Tuback, obj.Tyback, obj.Tvback,...
                                                        obj.Tuforward);
        end
        function obj = setupUncoupledUnconstrainedNonlinear(obj, A, B, mpcQ, mpcR, mheQ, mheR)
            obj.dynamicConstraints = MPCMHE_6dofutils.createFixedBodyLTIDynamics(A,B,obj.Tx0, obj.Tx, obj.Tuback, obj.Tuforward, obj.Td, obj.disturbType);
            obj.costFunction       = MPCMHE_Tcalcutils.objectiveUnconstrainedMPCMHENonlinearMeas(A, B, mpcQ, mpcR, mheQ, mheR, obj.Tx, obj.Tuback, obj.Tuforward, obj.Tyback, obj.Td, obj.Tvback, obj.disturbType);
            obj.opt_properties = MPCMHE_6dofutils.setupUnconstrainedOptimizationCells(obj.costFunction,obj.dynamicConstraints,...
                                                        obj.dMax, obj.uMax, obj.vMax,...
                                                        obj.Tx0, obj.Tx, obj.Td, obj.Tuback, obj.Tyback, obj.Tvback,...
                                                        obj.Tuforward);
        end
        function obj = setupLTICoupled(obj, A, B, C, mpcQ, mpcR, mheQ, mheR)
            obj.dynamicConstraints = MPCMHE_6dofutils.createCoupledLTIDynamics(A,B,obj.Tx0, obj.Tx, obj.attTx0, obj.attTx, obj.Tuback, obj.Tuforward, obj.Td, obj.disturbType);
            obj.sensorConstraints =  MPCMHE_Tcalcutils.createLinearSensorConstraint(C, obj.Tx, obj.Tvback, obj.Tyback, obj.backwardT);
            obj.costFunction = MPCMHE_Tcalcutils.objectiveMPCMHE(mpcQ,mpcR,mheQ,mheR,obj.Tx, obj.Tuforward, obj.Td, obj.Tvback);
            obj.opt_properties = MPCMHE_6dofutils.setupOptimizationCells(obj.costFunction,obj.dynamicConstraints, obj.sensorConstraints,...
                                                        obj.dMax, obj.uMax, obj.vMax,...
                                                        obj.Tx0, obj.Tx, obj.Td, obj.Tuback, obj.Tyback, obj.Tvback,...
                                                        obj.Tuforward);
        end
    end
end