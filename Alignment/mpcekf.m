classdef mpcekf
    properties
        maxIter = 20;
        compileFlag = '-Ofast';
        optimizer = @cmex2optimizeCS;
        % optimizer = @class2optimizeCS;
        opt

        window_mpcXs
        window_mpcXsAtt
        window_mpcUs
        window_mpcUsAtt

        ekf_hcw
        ekf_att

        forwardT

        x0
        att0
    end
    methods 
        function obj = init(obj, forwardHorizon, x0, att0)
            obj.forwardT = forwardHorizon;
            
            obj.window_mpcUs = zeros(6,obj.forwardT);
            obj.window_mpcXs = zeros(6,obj.forwardT);

            obj.window_mpcUsAtt = zeros(3,obj.forwardT);
            obj.window_mpcUs    = zeros(6,obj.forwardT);

            obj.x0 = x0;
            obj.att0 = att0;
        end

        % for general use
        function [control, attcontrol, state, stateAtt] = outputResults(obj)
            control = obj.window_mpcUs(:,1);
            attcontrol = obj.window_mpcUs(:,1);

            state = obj.x0;
            stateAtt = obj.att0;
        end
        function obj = estimate_and_control(obj, control, controlAtt, meas, measAtt)
            obj.ekf_hcw.estimateEKF(control, meas);
            obj.ekf_att.estimateKF(controlAtt, measAtt, eye(6));

            hcw_state = obj.ekf_hcw.outputResult();
            att_state = obj.ekf_att.outputResult();

            %setup mpc
            obj.x0 = hcw_state;
            obj.att0 = att_state;

            %% SHIFT WINDOWS
            obj.window_mpcXs = [obj.window_mpcXs(:,2:end), obj.window_mpcXs(:,end)];
            [dim,n] = size(obj.window_mpcUs);
            obj.window_mpcUs = [obj.window_mpcUs(:,2:end), zeros(dim,1)];

            obj.window_mpcXsAtt = [obj.window_mpcXsAtt(:,2:end), obj.window_mpcXsAtt(:,end)];
            [dim, n] = size(obj.window_mpcUsAtt);
            obj.window_mpcUsAtt = [obj.window_mpcUsAtt(:,2:end), zeros(dim,1)];

            setP_x0(obj.opt, obj.x0);
            setV_x(obj.opt,  obj.window_mpcXs);
            setV_u(obj.opt,  obj.window_mpcUs);
            
            setP_att0(obj.opt, obj.att0);
            setV_att(obj.opt, obj.window_mpcXsAtt);
            setV_u(obj.opt, obj.window_mpcUsAtt);

            mu0 = 1;
            saveIter = -1;
            addEyeHessian = [1e-5, 1e-5];
            [status, iter, time] = solve(obj.opt, mu0, int32(obj.maxIter), int32(saveIter), addEyeHessian);

            disp("Status:")
            disp(status)
            disp("Iterations:")
            disp(iter)

            [Jtotal, mpcXs, mpcAttXs, mpcUs, mpcAttUs] = getOutputs(obj.opt);

            disp("Jcost:")
            disp(Jtotal)

            obj.window_mpcXs = mpcXs;
            obj.window_mpcUs = mpcUs;
            obj.window_mpcXsAtt = mpcAttXs;
            obj.window_mpcUsAtt = mpcAttUs;
        end



        function obj = setupEKF(obj, cov0, covatt0, Q, R, Qatt, Ratt, tstep)
            [A_HCW, B_HCW] = utils.linearHCWDynamics(tstep);
            [Aatt,  Batt]  = utils.attitudeDiscrete(tstep);

            obj.ekf_hcw = EKF;
            obj.ekf_hcw = obj.ekf_hcw.init(obj.x0, cov0, A_HCW, B_HCW, Q, R);

            obj.ekf_att = EKF;
            obj.ekf_att = obj.ekf_att.init(obj.att0, covatt0, Aatt, Batt, Qatt, Ratt);
        end
        function obj = setupOptimizationCode(obj, tstep)
            statedim = 6;
            controldim = 6;
            Tvariable x0 [statedim,1];
            Tvariable x  [statedim,   obj.forwardT];
            Tvariable u  [controldim, obj.forwardT];

            statedim = 6;
            controldim = 3;
            Tvariable att0 [statedim, 1];
            Tvariable att  [statedim, obj.forwardT];
            Tvariable attu [controldim, obj.forwardT];

            [A_HCW, B_HCW] = utils.linearHCWDynamics(tstep);
            [Aatt,  Batt]  = utils.attitudeDiscrete(tstep);

            dynamicsConstraintsTranslational = tenscalc_utils.mpcDynamicsTranslational(A_HCW,B_HCW,x0,x,att0,att,u);
            dynamicsConstraintsAttitude      = tenscalc_utils.mpcDynamicsAttitude(Aatt, Batt, att0, att, attu);

            mpcQ = 1e3;
            mpcR = 1e1;
            Jcost = tenscalc_utils.mpcObjectiveTranslational(mpcQ,mpcR,x,u);

            mpcQ = 1e3;
            mpcR = 1e1;
            attJcost = tenscalc_utils.mpcObjectiveAttitude(mpcQ,mpcR,att,attu);

            Jtotal = Jcost + attJcost;

            uMax = 0.1;
            attuMax = 0.3;

            classname = obj.optimizer(...
                'classname', 'tmpC',...
                'optimizationVariables', {x,u,att,u},...
                'objective',Jtotal,...
                'constraints',{...
                               dynamicsConstraintsTranslational,...
                               dynamicsConstraintsAttitude,...
                               u>=-1e-5*Tones(size(u)),...
                               u<=uMax*Tones(size(u)),...
                               attu <= attuMax*Tones(size(attu)),...
                               attu >= -attuMax*Tones(size(attu))...
                               },...
                'parameters',{x0, att0},...
                'outputExpressions',{Jtotal, x, att, u, attu},...
                'smallerNewtonMatrix',true,...
                'addEye2Hessian',true,...
                'muFactorAggressive',0.99,...
                'muFactorConvservative',0.99,...
                'desiredDualityGap',0.1,...
                'maxIter',obj.maxIter,...
                'gradTolerance',1e-3,...
                'skipAffine',true,...
                'delta',2,...
                'compilerOptimization',obj.compileFlag,...
                'allowSave',true,...
                'profiling',false,...
                'solverVerboseLevel',1,...
                'verboseLevel',1);
            
            obj.opt = feval(classname);
        end
    end
end