classdef mpc_mhe
    properties
        maxIter = 20;
        compileFlag = '-Ofast';
        optimizer = @cmex2optimizeCS;
        % optimizer = @class2optimizeCS;
        opt

        opt_est

        window_mpcXs
        window_mpcXsAtt
        window_mpcUs
        window_mpcUsAtt


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
            obj.window_mpcXsAtt = zeros(6,obj.forwardT);

            obj.x0 = x0;
            obj.att0 = att0;
        end

        % for general use
        function [control, attcontrol, state, stateAtt] = outputResults(obj)
            control = obj.window_mpcUs(:,1);
            attcontrol = obj.window_mpcUsAtt(:,1);

            state = obj.x0;
            stateAtt = obj.att0;
        end
        function obj = estimate_and_control(obj, control, controlAtt, meas, measAtt)
            % P = [eye(3), -eye(3)];

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
            setV_attu(obj.opt, obj.window_mpcUsAtt);

            mu0 = 1;
            saveIter = -1;
            addEyeHessian = [1e-5; 1e-5];
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

        function obj = setupOptimizationEstimation(obj, tstep)
            statedim = 6;
            measdim = 3;
            controldim = 6;


            Tvariable x0past [statedim,1];
            Tvariable xpast  [statedim,   obj.backwardT];
            Tvariable upast  [controldim, obj.backwardT];
            Tvariable ypast  [measdim, obj.backwardT];
            Tvariable dpast  [controldim, obj.backwardT];

            statedim = 6;
            measdim = 6;
            controldim = 3;
            Tvariable att0past [statedim, 1];
            Tvariable attpast  [statedim, obj.backwardT];
            Tvariable attupast [statedim, obj.backwardT];
            Tvariable attypast [statedim, obj.backwardT];


            [A_HCW, B_HCW] = utils.linearHCWDynamics(tstep);
            [Aatt,  Batt]  = utils.attitudeDiscrete(tstep);

            dynamicsConstraintsTranslational = tenscalc_utils.mheDynamicsTranslational(A_HCW, B_HCW, x0past, xpast, att0past, attpast, upast, dpast);
            dynamicsConstraintsAttitude      = tenscalc_utils.mheDynamicsAttitude(Aatt, Batt, att0past, attpast, attupast);

            mheQ = 1e1;
            mheR = 1e1;
            Jcost = tenscalc_utils.mheObjectiveTranslational(mheQ,mheR,xpast,ypast,dpast);

            mheQ = 1e1;
            mheR = 1e1;
            attJcost = tenscalc_utils.mheObjectiveAttitude(mheQ, mheR, attpast, attypast);

            Jtotal = Jcost + attJcost;

            dMax = 0.05;

            classname = obj.optimizer(...
                'classname','tmpC_est',...
                'optimizationVariables', {xpast, attpast},...
                'objective',Jtotal,...
                'constraints',{
                                dynamicsConstraintsTranslational,...
                                dynamicsConstraintsAttitude,...
                                dpast >= dMax*Tones(size(dpast)),...
                                dpast <= -dMax*Tones(size(dpast)),...
                              },...
                'parameters', {x0past, att0past},...
                'outputExpressions', {Jtotal, xpast, attpast},...
                'smallerNewtonMatrix',false,...
                'addEye2Hessian',true,...
                'muFactorAggressive',0.99,...
                'muFactorConservative',0.99,...
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
            obj.opt_est = feval(classname);
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
                'optimizationVariables', {x,u,att,attu},...
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
                'smallerNewtonMatrix',false,...
                'addEye2Hessian',true,...
                'muFactorAggressive',0.99,...
                'muFactorConservative',0.99,...
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