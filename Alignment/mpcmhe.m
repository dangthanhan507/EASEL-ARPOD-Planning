classdef mpcmhe
    properties
        maxIter = 40;
        compileFlag = '-Ofast';
        % optimizer = @cmex2equilibriumLatentCS;
        optimizer = @class2equilibriumLatentCS;
        opt


        %windows (translational)
        window_mheXs % mhe states
        window_mpcXs % mpc states
        window_mheUs % mhe control past
        window_mpcUs % mpc control future
        window_mheYs % mhe measurement past

        Dmat         % Time-Invariant disturbance
        window_ds % Time-varying additive disturbance

        %windows (attitudinal)
        window_mheXsAtt % mhe states
        window_mpcXsAtt % mpc states
        window_mheUsAtt % mhe control past
        window_mpcUsAtt % mpc control future
        window_mheYsAtt % mhe measurement past

        window_mheVsAtt % mhe sensor noise past

        time_ctr = 1;
        backwardT
        forwardT
        x0
        att0
    end
    methods
        function obj = init(obj, backwardHorizon, forwardHorizon, x0, att0)

            obj.backwardT = backwardHorizon;
            obj.forwardT  = forwardHorizon;

            %windows (translational)
            obj.window_mheXs = zeros(6, backwardHorizon); % mhe states
            obj.window_mpcXs = zeros(6, forwardHorizon); % mpc states
            obj.window_mheUs = zeros(6, backwardHorizon); % mhe control past
            obj.window_mpcUs = zeros(6, forwardHorizon); % mpc control future
            obj.window_mheYs = zeros(3, backwardHorizon); % mhe measurement past

            obj.Dmat = eye(6);         % Time-Invariant disturbance
            obj.window_ds = zeros(6, backwardHorizon + forwardHorizon); % Time-varying additive disturbance

            %windows (attitudinal)
            obj.window_mheXsAtt = zeros(6, backwardHorizon); % mhe states
            obj.window_mpcXsAtt = zeros(6, forwardHorizon); % mpc states
            obj.window_mheUsAtt = zeros(3, backwardHorizon); % mhe control past
            obj.window_mpcUsAtt = zeros(3, forwardHorizon); % mpc control future
            obj.window_mheYsAtt = zeros(6, backwardHorizon); % mhe measurement past

            obj.window_mheVsAtt = zeros(6, backwardHorizon); % mhe sensor noise past 


            obj.x0 = x0;
            obj.att0 = att0;
        end

        % for general use
        %
        function [control, attcontrol, state, stateAtt] = outputResults(obj)
            control = obj.window_mpcUs(:,1);
            attcontrol = obj.window_mpcUsAtt(:,1);
            state = obj.window_mheXs(:,end);
            stateAtt = obj.window_mheXsAtt(:,end);
        end
        function obj = estimate_and_control(obj, control, controlAtt, meas, measAtt)

            
            if obj.time_ctr <= obj.backwardT
                %only for first-time use when we don't have all the measurements yet
                obj.window_mheYs(:,obj.time_ctr) = meas;
                obj.window_mheUs(:,obj.time_ctr) = control;

                obj.window_mheYsAtt(:,obj.time_ctr) = measAtt;
                obj.window_mheUsAtt(:,obj.time_ctr) = controlAtt;

                obj.window_mpcUs = zeros(6, obj.forwardT);
                obj.window_mpcUsAtt = zeros(3, obj.forwardT);

                obj.time_ctr = obj.time_ctr + 1;
            else
                % shift all the windows

                %%%%% TRANSLATIONAL SHIFT %%%%%%%%
                obj.x0 = obj.window_mheXs(:,1);
                obj.window_mheXs = [obj.window_mheXs(:,2:end), obj.window_mpcXs(:,1)];

                obj.window_mpcXs = [obj.window_mpcXs(:,2:end), obj.window_mpcXs(:,end)];

                [dim,n] = size(obj.window_ds);
                obj.window_ds = [obj.window_ds(:,2:end), zeros(dim,1)];

                [dim,n] = size(obj.window_mpcUs);
                obj.window_mheUs = [obj.window_mheUs(:,2:end), control];
                obj.window_mpcUs = [obj.window_mpcUs(:,2:end), zeros(dim,1)];

                obj.window_mheYs = [obj.window_mheYs(:,2:end), meas];

                %%%%% ATTITUDINAL SHIFT %%%%%%%%
                obj.att0 = obj.window_mheXsAtt(:,1);
                obj.window_mheXsAtt = [obj.window_mheXsAtt(:,2:end), obj.window_mpcXsAtt(:,1)];
                obj.window_mpcXsAtt = [obj.window_mpcXsAtt(:,2:end), obj.window_mpcXsAtt(:,end)];

                [dim, n] = size(obj.window_mpcUsAtt);
                obj.window_mheUsAtt = [obj.window_mheUsAtt(:,2:end), controlAtt];
                obj.window_mpcUsAtt = [obj.window_mpcUsAtt(:,2:end), zeros(dim,1)];
                obj.window_mheYsAtt = [obj.window_mheYsAtt(:,2:end), measAtt];
            end

            if obj.time_ctr > obj.backwardT
                %main code
                % disp(obj.window_mheYs);
                % disp(obj.window_mheYsAtt);
                
                % translational variable setting
                setP_uback(obj.opt, obj.window_mheUs);
                setP_ypast(obj.opt, obj.window_mheYs);
                setP_x0(obj.opt, obj.x0);
                setV_u(obj.opt, obj.window_mpcUs);
                setV_D(obj.opt, obj.Dmat);
                setV_d(obj.opt, obj.window_ds);
                setV_x(obj.opt, [obj.window_mheXs, obj.window_mpcXs]);

                % attitudinal variable setting
                setP_attuback(obj.opt, obj.window_mheUsAtt);
                setP_attypast(obj.opt, obj.window_mheYsAtt);
                setP_att0(obj.opt, obj.att0);
                setV_attu(obj.opt, obj.window_mpcUsAtt);
                setV_att(obj.opt, [obj.window_mheXsAtt, obj.window_mpcXsAtt]);
                setV_attvback(obj.opt, obj.window_mheVsAtt);

                mu0 = 1;
                saveIter = -1;
                addEyeHessian = [5e1;5e1];
                [status, iter, time] = solve(obj.opt, mu0, int32(obj.maxIter), int32(saveIter), addEyeHessian);

                disp("Status:")
                disp(status)
                disp("Iterations:")
                disp(iter)

                [Jtotal, Xs, XsAtt, mpcUs, mpcUsAtt, D, d, mheVsAtt] = getOutputs(obj.opt);

                disp("Jcost:")
                disp(Jtotal)                

                mheXs = Xs(:,1:obj.backwardT);
                mpcXs = Xs(:,obj.backwardT+1:end);

                obj.window_mpcXs = mpcXs;
                obj.window_mheXs = mheXs;
                obj.window_mpcUs = mpcUs;
                obj.Dmat = D;
                obj.window_ds = d;

                mheXsAtt = XsAtt(:,1:obj.backwardT);
                mpcXsAtt = XsAtt(:,obj.backwardT+1:end);

                obj.window_mpcXsAtt = mpcXsAtt;
                obj.window_mheXsAtt = mheXsAtt;
                obj.window_mpcUsAtt = mpcUsAtt;
                obj.window_mheVsAtt = mheVsAtt;
            end
        end



        % behind the scenes work (optimization)

        %function to create C code for optimization (setup before use)
        function obj = setupOptimizationCode(obj, tstep)
            backHorizon = obj.backwardT;
            forwardHorizon = obj.forwardT;

            statedim = 6;
            controldim = 6;
            meas_dim = 3;

            %translation tcalc variables
            Tvariable x0    [statedim,1];
            Tvariable x     [statedim, backHorizon + forwardHorizon];
            Tvariable uback [controldim, backHorizon];
            Tvariable ypast [meas_dim, backHorizon];
            Tvariable u     [controldim, forwardHorizon];
            Tvariable D [controldim, controldim];
            Tvariable d [controldim, backHorizon + forwardHorizon];

            statedim = 6;
            controldim = 3;
            meas_dim = 6;

            %attitude tcalc variables
            Tvariable att0    [statedim,1];
            Tvariable att     [statedim, backHorizon + forwardHorizon];
            Tvariable attuback [controldim, backHorizon];
            Tvariable attvback [statedim, backHorizon];
            Tvariable attypast [statedim, backHorizon];
            Tvariable attu     [controldim, forwardHorizon];

            [A_HCW, B_HCW] = utils.linearHCWDynamics(tstep);
            [Aatt,  Batt]  = utils.attitudeDiscrete(tstep);


            % translational setup
            dynamicsConstraintsTranslational = tenscalc_utils.mpcmheDynamicsTranslational(A_HCW, B_HCW,...
                                                                             x0, x, att0, att, uback, u, D, d);

            mpcQ = 1e6;
            mpcR = 1e1;
            mheQ = 1e6;
            mheR = 1e6;
            Jcost = tenscalc_utils.mpcmheObjectiveTranslational(mpcQ,mpcR,mheQ,mheR,x,u,ypast,d,D);


            % attitudinal setup
            mpcQ = 1e1;
            mpcR = 1e1;
            mheQ = 1e1;
            mheR = 1e1;
            dynamicsConstraintsAttitude = tenscalc_utils.mpcmheDynamicsAttitude(Aatt, Batt, att0, att, attuback, attu);
            sensorConstraintsAttitude   = tenscalc_utils.mpcmheSensorConstraintsAttitude(att, attvback, attypast);
            attJcost = tenscalc_utils.mpcmheObjectiveAttitude(mpcQ, mpcR, mheQ, mheR, att, attu, attvback);

            Jtotal = Jcost + attJcost;

            uMax = 0.1;
            attuMax = 0.3;
            DMax = 1.4;
            dMax = 0.1;
            attvMax = 0.8;

            classname = obj.optimizer(...
                'classname','tmpC',...
                'P1objective',Jtotal,...
                'P2objective',-Jtotal,...
                'P1optimizationVariables',{u, attu},...
                'P2optimizationVariables',{D, d, attvback},...
                'latentVariables',{x,att},...
                'P1constraints',{...
                                 u >= -1e-3*Tones(size(u)),...
                                 u <= uMax*Tones(size(u)),...
                                 attu <= attuMax*Tones(size(attu)),...
                                 attu >= -attuMax*Tones(size(attu))...
                                 },...
                'P2constraints',{...
                                 sensorConstraintsAttitude,...
                                 D <= DMax*Tones(size(D)),...
                                 D >= -DMax*Tones(size(D)),...
                                 d >= -dMax*Tones(size(d)),...
                                 d <= dMax*Tones(size(d)),...
                                 attvback >= -attvMax*Tones(size(attvback)),...
                                 attvback <= attvMax*Tones(size(attvback))...
                                 },...
                'latentConstraints',{dynamicsConstraintsTranslational, dynamicsConstraintsAttitude},...
                'outputExpressions',{Jtotal, x, att, u, attu, D, d, attvback},...
                'parameters',{uback, attuback, ypast, attypast, x0, att0},...
                'smallerNewtonMatrix',false,...
                'addEye2Hessian',true,...
                'muFactorAggressive',0.99,...
                'muFactorConservative',0.99,...
                'desiredDualityGap',1.0,...
                'maxIter',obj.maxIter,...
                'gradTolerance',1.0,...
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