classdef MPCMHE_thruster
    properties
        %TensCalc Specific parameters
        version
        codeType
        compilerFlags

        % window of stuff to consider
        window_mhecontrols
        window_mpccontrols

        window_measurements

        window_mhestates
        window_mpcstates

        window_controlDisturbances
        window_measError

        window_phase

        % pertinent information regarding algorithm
        x0
        backwardT
        forwardT

        % tens calc optimizer
        opt

        %T variables from tenscalc
        Tx0
        Tx
        Td
        Tuback
        Tyback
        Tvback
        Tuforward

        % misalignment variables
        % not there yet
        A_x %matrix bias
        b_x %additive bias
    end
    methods (Static)
        %padding solution for MHE
        function modified_sense = senseModify(measurement)
            [dim,n] = size(measurement);
            modified_sense = measurement;
            if dim == 2
                modified_sense(3) = 0;
            end
        end
        function modified_control = controlModify(control, A, b)
            modified_control = A*control + b;
        end
        function [s,b] = estimateBias(d_ks, u_ks)
        end
    end
    methods
        %setup tenscalc global params
        function obj = setupTensCalc(obj)
            % obj.version = @cmex2equilibriumLatentCS;
            obj.version = @class2equilibriumLatentCS;
            obj.codeType = 'C';
            obj.compilerFlags = '-O1';
        end
        
        %initialize class
        function obj = init(obj, x0_, backwardT, forwardT, tstep)
            [dim,n] = size(x0_); 
            control_dim = dim/2;
            meas_dim = 3;

            obj.x0 = x0_;
            obj.backwardT = backwardT;
            obj.forwardT = forwardT;

            Tvariable x0a [dim,1];
            Tvariable x [dim, backwardT + forwardT];
            Tvariable d [control_dim, backwardT + forwardT];
            Tvariable uback [control_dim, backwardT];
            Tvariable yback [meas_dim, backwardT];
            Tvariable vback [meas_dim, backwardT];
            Tvariable uforward [control_dim, forwardT];

            obj.Tx0 = x0a;
            obj.Tx = x;
            obj.Td = d;
            obj.Tuback = uback;
            obj.Tyback = yback;
            obj.Tuforward = uforward;
            obj.Tvback = vback;
        end
        function obj = initWindows(obj, window_mhecontrols, window_mpccontrols, window_measurements, window_mhestates, window_mpcstates, window_phases)
            obj.window_mpccontrols = window_mpccontrols;
            obj.window_mhecontrols = window_mhecontrols;
            obj.window_measurements = window_measurements;
            obj.window_mhestates = window_mhestates;
            obj.window_mpcstates = window_mpcstates;

            obj.window_measError = 100*ones(3,obj.backwardT);
            obj.window_controlDisturbances = 100*ones(3,obj.backwardT);

            obj.window_phase = window_phases;
        end

        function obj = windowShift(obj, meas, control, tstep, phase)
            [A,B] = ARPOD_Dynamics.linearHCWDynamics(tstep);
            [dim,num] = size(obj.window_measurements);


            obj.x0 = obj.window_mhestates(:,1);
            obj.window_mhestates = [obj.window_mhestates(:,2:end), A*obj.window_mhestates(:,end) + B*control];
            obj.window_mpcstates = [obj.window_mpcstates(:,2:end), A*obj.window_mpcstates(:,end)];

            obj.window_measurements = [obj.window_measurements(:,2:end), MPCMHE_thruster.senseModify(meas)];
            obj.window_mhecontrols = [obj.window_mhecontrols(:,2:end), control];

            obj.window_mpccontrols = [obj.window_mpccontrols(:,2:end), zeros(3,1)];

            obj.window_controlDisturbances = [obj.window_controlDisturbances(:,2:end), 100*ones(6,1)];

            obj.window_phase = [obj.window_phase(2:end-1), phase];

            state_meas = ARPOD_Sensor.measure(obj.window_mhestates(:,end));
            target_meas = MPCMHE_thruster.senseModify(meas);
            if target_meas(3) == 0
                state_meas(3) = 0;
            end
            error = target_meas - state_meas;
            obj.window_measError = [obj.window_measError(:,2:end), error];
        end        

        %setting up constraints
        function xDynamics = setupDynamics(obj, tstep)
            x_k = [obj.Tx0, obj.Tx(:,1:end-1)];
            u = [obj.Tuback, obj.Tuforward];
            [A,B] = ARPOD_Dynamics.linearHCWDynamics(tstep);
            xDynamics = ( obj.Tx == A*x_k + B*( u + obj.Td) );
        end
        function yConstr = setupMeasConstr(obj)
            x_k = obj.Tx(:,1:obj.backwardT);

            state_meas = [];
            
            for i = 1:obj.backwardT
                state_meas = [ state_meas, ARPOD_Sensor.Tmeasure(x_k(:,i)) ];
                
                
                if obj.window_phase(i) == 1
                    state_meas(:,i) = diag([1,1,0])*state_meas(:,i);
                end
            end
            
            disp("ALO")
            yConstr = (obj.Tyback == state_meas + obj.Tvback);
        end

        %setting up objective
        function J = objectiveF(obj,mpcQ,mpcR,mheQ,mheR)
            J = 10*norm2(obj.Tx(:,obj.backwardT+1:obj.forwardT));
            J = J + 10*norm2(obj.Tuforward);
            J = J - 10*norm2(obj.Td) - 100*norm2(obj.Tvback);

        end

        %setting up optimizer
        function obj = tensCalcOpt(obj, mpcQ, mpcR, mheQ, mheR, tstep, vMax, dMax, uMax)
            %setup objective
            J = obj.objectiveF(mpcQ,mpcR,mheQ,mheR);

            %setup constraints
            xDynamics = obj.setupDynamics(tstep);
            measConstr = obj.setupMeasConstr();

            maxIter = 200;
            allowSave = true;
            %add all of this into TensCalc
            %%%%
            %%%% PROBLEM IS IN THE MEASUREMENT CONSTRAINT RIGHT NOW
            %%%%
            %%%% 
            classname=obj.version(...
                'classname','tmpC_target_chaser_main',...
                'method','primalDual',...
                'P1objective',J,...
                'P2objective',-J,...
                'P1optimizationVariables',{obj.Tuforward},...
                'P2optimizationVariables',{obj.Td, obj.Tvback, obj.Tx0},...
                'latentVariables',{obj.Tx},...
                'P1constraints',{obj.Tuforward<=uMax*Tones(size(obj.Tuforward)),...
                                 obj.Tuforward>=-uMax*Tones(size(obj.Tuforward))},...
                'P2constraints',{obj.Td<=dMax*Tones(size(obj.Td)),...
                                 obj.Td>=-dMax*Tones(size(obj.Td)), ...
                                 obj.Tvback <= vMax*Tones(size(obj.Tvback)),...
                                 obj.Tvback >= -vMax*Tones(size(obj.Tvback)),...
                                 measConstr},...
                'latentConstraints',{xDynamics},...
                'outputExpressions',{J,obj.Tuforward,obj.Td,obj.Tx0,obj.Tx, obj.Tvback},...
                'parameters',{obj.Tuback,obj.Tyback},...
                'muFactorAggressive',.5,...
                'muFactorConservative',.99,...
                'desiredDualityGap',0.1,...
                'maxIter',maxIter,...
                'addEye2Hessian',0*1e-9,...
                'gradTolerance',1e-3,...
                'skipAffine',true,...
                'delta',2,...
                'codeType',obj.codeType,...
                'compilerOptimization',obj.compilerFlags,...
                'allowSave',allowSave,...
                'profiling',false,...
                'solverVerboseLevel',1,... 
                'verboseLevel',1); 
            obj.opt = feval(classname);
        end
        function [mheXs, mheDs, vs, mpcXs, mpcUs]  = optimize(obj)
            %uses same names as defined in the init

            %define parameters
            setP_uback(obj.opt,obj.window_mhecontrols);
            setP_yback(obj.opt,obj.window_measurements);

            %define variables
            setV_uforward(obj.opt, obj.window_mpccontrols);
            setV_d(obj.opt, obj.window_controlDisturbances);
            setV_x0a(obj.opt, obj.x0);

            combined_x = [obj.window_mhestates, obj.window_mpcstates];
            setV_x(obj.opt, combined_x);

            mu0 = 1;
            maxIter = 200;
            saveIter = -1;

            solve(obj.opt,mu0,int32(maxIter), int32(saveIter));
            [Jcost, mpcUs, mheDs, x0s, xs, vs] = getOutputs(obj.opt);

            mheXs = xs(:,1:obj.backwardT);
            mpcXs = xs(:,obj.backwardT+1:end);
        end
        function obj = estimate(obj, meas, control, tstep)
            obj.windowShift(meas,control,tstep);
            obj.optimize();

            %estimate
        end
    end
end