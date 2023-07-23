classdef MpcMheThruster
    %{
        MPCMHE problem working with the ARPOD problem.
        A full abstraction that should allow for modifications


    %}
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

        hcwA % discrete HCW dynamics matrix for state
        hcwB % discrete HCW dynamics matrix for control input

        % misalignment variables
        % not there yet
        % TODO: lets get there
        A_x %matrix bias
        b_x %additive bias
    end
    methods
        function obj = init(obj, x0_, backwardT, forwardT, tstep, meas_dim, control_dim)
            %{
                Description:
                ------------
                initializes MPCMHE object.
                setups up the parameters for the MPCMHE that are invariant every problem such as
                    -> backward time horizon
                    -> forward time horizon
                    -> seconds per discrete timestep

                In short, this should be used as the CONSTRUCTOR of the class. That means call this first 
                before calling anything else.
                
                Parameters
                -----------
                @params backwardT: backward facing time window (scalar)
                @params forwardT:  forward facing time window (scalar)
                @params tstep:     number of seconds per timestep (scalar)
            %}

            obj.version = @class2equilibriumLatentCS;
            obj.codeType = 'C';
            obj.compilerFlags = '-O1';

            [dim,n] = size(x0_);

            obj.x0 = x0_;
            obj.backwardT = backwardT;
            obj.forwardT  = forwardT;

            %{
                TensCalc specific variables:
                ----------------------------
                    -> lists all variables that will be used during optimization.
                    -> x0a (intiial state)
                    -> x (combined states from mpc and mhe to optimize)
                    -> d (noise) added to control dim
                    -> uback (past control inputs)
                    -> vback (mhe requires you to model noise as a variable.... this is v)
            %}
            Tvariable x0a   [dim,1];
            Tvariable x     [dim, backwardT + forwardT];
            Tvariable d     [control_dim, backwardT + forwardT];
            Tvariable uback [control_dim, backwardT];
            Tvariable vback [meas_dim, backwardT];
            Tvariable u     [control_dim, forwardT];
            Tvariable ypast [meas_dim, backwardT];

            %set class members to be equal to Tvariables to allow access throughout class instance.
            obj.Tx0 = x0a;
            obj.Tx = x;
            obj.Td = d;
            obj.Tuback = uback;
            obj.Tuforward = u;
            obj.Tvback = vback;
            obj.Tyback = ypast;


            % define dynamics utilized in the controller
            [obj.hcwA,obj.hcwB] = ARPOD_Dynamics.linearHCWDynamics(tstep);
        end
        function xDynamics = setupDynamicConstraints(obj)
            %{
                Description:
                -------------
                The MpcMhe problem is a huge optimization procedure. In the formulation, we are 
                concerned with the dynamical constraints
                
                x_k+1 = A*x_k + B*u_k

                This is where we will define this in a way we can interface with Tenscalc

                You can consider this a helper function for the optimize() function.

                Parameters:
                ------------
                N/A

                Returns:
                --------
                returns dynamic contraint instance for tenscalc
            %}

            %the dynamic constraints occur for everything but last element
            xk = [obj.Tx0, obj.Tx(:,1:end-1)];
            u = [obj.Tuback, obj.Tuforward];

            %return
            xDynamics = (obj.Tx == obj.hcwA*xk + obj.hcwB*(u+obj.Td));
        end
        function yConstr = setupSensorConstr(obj)
            %{
                Description:
                ------------
                Sets up the measurement constraints specifically for the MHE portion of the
                estimation process.

                Normally measurement constraints are in this form
                yk = g(xk) + vk

                yk is sensor measurement
                xk is the state
                g(*) is the function to convert the latent state into a sensor measurement
                vk is the error between g() and yk....

                Parameters:
                -----------
                N/A

                Returns:
                ---------
                returns constraint instance for TensCalc for measurements


                RESEARCH DETAILS: (IMPORTANT)
                ==============================
                    -> we intend to use an indentity matrix for g
                    -> yk = Ix + vk
                    -> this is because we tried using the ARPOD sensor measurements
                        and it didn't work out too well. (infeasible optimization)
                    
                    -> if there is interest in using ARPOD code talk to An. He has done this before
                        and can show you how to reproduce the code.
            %}

            xk = obj.Tx(:,1:obj.backwardT);
            [dim,n] = size(obj.x0);

            C = eye(dim); % measurement function

            yConstr = (obj.Tyback == C*xk + obj.Tvback);
        end
        %setting up objective
        function J = objectiveF(obj)
            %{
                Description:
                ------------
                Creates an objective function to minimize/maximize with TensCalc

                Cost is formulated as

                J = (x^T @ Q @ x) + (u^T @ R @ u)


                Parameters:
                -----------
                N/A
                TODO: might be good to add Q,R matrices for this

                Returns:
                ----------
                returns objective function instance for TensCalc format
            %}
            J = 10*norm2(obj.Tx(:,obj.backwardT+1:obj.forwardT));
            J = J + 10*norm2(obj.Tuforward);
            J = J - 10*norm2(obj.Td) - 100*norm2(obj.Tvback);
        end
        function obj = setupOptimize(obj, vMax, dMax, uMax, maxIter)
            %{
                Description:
                ------------
                Sets up optimization problem using objective function, constraints, and
                current horizon of state/measurements/control inputs

                Parameters:
                ------------
                @params vMax: the max absolute value that the vk variable is allowed to be at.
                                this means -vMax <= vk <= vMax
                @params dMax: the max absolute value that the dk variable is allowed to be at.
                                this means -dMax <= dk <= dMax
                @params uMax: the max absolute value that the uk variable is allowed to be at.
                                this means -uMax <= uk <= uMax

                @params maxIter: the max iterations the MPC-MHE is allowed to do.

                Returns:
                -------
                class instance with full optimization setup.
            %}
            J = obj.objectiveF();

            xDynamics = obj.setupDynamicConstraints();
            yConstr = obj.setupSensorConstr();

            allowSave = true;

            %{
                creating optimization problem for MPC-MHE

            %}
            classname=obj.version(...
                'classname','tmpC_target_chaser_main',...
                'P1objective',J,...
                'P2objective',-J,...
                'P1optimizationVariables',{obj.Tuforward},...
                'P2optimizationVariables',{obj.Td, obj.Tvback, obj.Tx0},...
                'latentVariables',{obj.Tx},...
                'P1constraints',{obj.Tuforward<=uMax*Tones(size(obj.Tuforward)),...
                                 obj.Tuforward>=-uMax*Tones(size(obj.Tuforward))},...
                'P2constraints',{obj.Td<=dMax*Tones(size(obj.Td)),...
                                 obj.Td>=-dMax*Tones(size(obj.Td)),...
                                 obj.Tvback <= vMax*Tones(size(obj.Tvback)),...
                                 obj.Tvback >= -vMax*Tones(size(obj.Tvback)),...
                                 yConstr},...
                'latentConstraints',{xDynamics},...
                'outputExpressions',{J,obj.Tuforward,obj.Td,obj.Tx0,obj.Tx, obj.Tvback},...
                'parameters',{obj.Tuback,obj.Tyback},...
                'muFactorAggressive',.5,...
                'muFactorConservative',.99,...
                'desiredDualityGap',0.1,...
                'maxIter',maxIter,...
                'gradTolerance',1e-3,...
                'skipAffine',true,...
                'delta',2,...
                'compilerOptimization',obj.compilerFlags,...
                'allowSave',allowSave,...
                'profiling',false,...
                'solverVerboseLevel',1,... 
                'verboseLevel',1); 
            obj.opt = feval(classname); %saves this work into the obj.opt variable
            %now you only need to call an equivalent of solve() to run this
        end
        function obj = setupWindows(obj, window_states, window_measurements, window_controls)
            %{

            %}
            [control_dim, n] = size(window_controls);
            [state_dim, n] = size(window_states);
            [meas_dim, n] = size(window_measurements);


            %setup mhe windows 
            obj.window_mhestates = window_states;
            obj.window_mhecontrols = window_controls;
            obj.window_measurements = window_measurements;

            %setup mhe disturbance window
            obj.window_controlDisturbances = zeros(control_dim, obj.forwardT + obj.backwardT);
            obj.window_measError = zeros(meas_dim, obj.backwardT);

            %setup mpc windows
            obj.window_mpccontrols = zeros(control_dim, obj.forwardT);
            obj.window_mpcstates = zeros(state_dim, obj.forwardT);
 
        end
        function obj = shiftWindows(obj, measurement, controls)
            %{
                TODO: write this
            %}
            error("TODO: implement this!");
        end
        function [mheXs, mheDs, vs, mpcXs, mpcUs] = optimize(obj)
            %{
                Description:
                ------------
                Performs optimization using the problem that has been setup in this class
                using TensCalc. Afterward, it returns all of the results of the optimization including
                whether it succeeded.

                Parameters:
                -----------
                N/A

                Returns:
                --------
                @return mheXs: states estimated by MHE [N, backwardT]
                @return mheDs: control input error estimated by MHE [N, backwardT]
                @return vs   : measurement error estimated by MHE [C, backwardT]
                @return mpcXs: states for optimal path created by MPC [N, forwardT]
                @return mpcUs: control inputs used for optimal path created by MPC [M, forwardT]

                N: dimension of state
                M: dimension of control input
                C: dimension of sensor measurement
            %}

            %define parameters used in TensCalc optimizer
            setP_uback(obj.opt, obj.window_mhecontrols);
            setP_ypast(obj.opt, obj.window_measurements);
            
            %define variables to optimize
            %NOTE: this can be used for warm-starting....
            setV_u(obj.opt, obj.window_mpccontrols);
            setV_d(obj.opt, obj.window_controlDisturbances);
            setV_x0a(obj.opt, obj.x0);
            setV_x(obj.opt, [obj.window_mhestates, obj.window_mpcstates]);
            setV_vback(obj.opt, obj.window_measError);

            mu0 = 1;
            maxIter = 1;
            saveIter = -1;

            [status, iter, time] = solve(obj.opt, mu0, int32(maxIter), int32(saveIter));
            [Jcost, mpcUs, mheDs, x0s, xs, vs] = getOutputs(obj.opt); %get results

            mheXs = xs(:,1:obj.backwardT);
            mpcXs = xs(:,obj.backwardT+1:end);
        end
    end
end