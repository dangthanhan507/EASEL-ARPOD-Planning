% ============================================================
% MpcMheThruster.m
% 
% This code centralizes in the utilization of MPC-MHE for
% the thruster misalignment problem. MHE performs state
% estimation while MPC performs control optimization.
% This is a research implementation using TensCalc as its
% base library to perform min-max optimization.
%
%
% Author: An Dang, David Copp
%
% ============================================================

classdef MpcMheThruster
    properties
        %TensCalc Specific parameters
        version
        codeType
        compilerFlags

        % window of stuff to consider

        %mhe windows
        window_mhecontrols %parameter in the estimation
        window_measurements%parameter in the estimation
        window_mhestates   %decision variable in optimization
        window_measError   %decision variable in optimization
        window_controlDisturbances % decision variable in optimization

        %mpc windows
        window_mpccontrols %decision variable in optimization
        window_mpcstates   %decision variable in optimiation


        % pertinent information regarding algorithm
        x0 %this is a variable we set to not change during optimization of mpc-mhe
        backwardT
        forwardT

        % tens calc optimizer
        opt

        %Tenscalc variables
        Tx0
        Tx
        Td
        Tuback
        Tyback
        Tvback
        Tuforward

        Ak % invariant discrete dynamics matrix for state
        Bk % invariant discrete dynamics matrix for control input

        % misalignment variables
        % TODO: not there yet... lets get there
        A_x %matrix bias
        b_x %additive bias
    end
    methods
        function obj = init(obj, x0_, backwardT, forwardT, A, B, meas_dim, control_dim)
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
                @params backwardT:   backward facing time window (scalar)
                @params forwardT:    forward facing time window (scalar)
                @params A        :   propagation dynamics matrix for the state
                @params B        :   propagation dynamics matrix for the control input
                @params meas_dim:    scalar representing vector length of sensor measurements
                @params control_dim: scalar representing vector length of control inputs
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
            obj.Ak = A;
            obj.Bk = B;
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
            xDynamics = (obj.Tx == obj.Ak*xk + obj.Bk*(u+obj.Td));
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
                    
                    -> if there is interest in using ARPOD Sensor code talk to An. He has 
                        done this before and can show you how to reproduce the code.
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
        function obj = shiftWindows(obj, measurement, control)
            %{
                Description:
                ------------
                Everything in the MPC-MHE class relies on the maintaining of
                its windows data structures. 

                We opt to shift the window data structures everytime we receive
                a measurement and control input.

                Parameters:
                -----------
                @params measurement: (M,1) vector representing incoming sensor measurement
                @params control    : (D,1) vector representing recent control input used

                M: dimension of sensor measurement.
                D: dimension of control input.

                Returns:
                --------
                class instance with its window data structures shifted.


                Access Notes:
                --------------
                Other than setupWindows()...,
                Ideally this will be the only function called before we call 
                setupOptimize() and optimize().
            %}
            

            obj.x0 = obj.window_mhestates(:,1);

            %shift and use obj.window_mpcstates as a "predicted variable" for warm-start
            obj.window_mhestates = [obj.window_mhestates(:,2:end), obj.window_mpcstates(:,1)];
            %propagate zero control input thrust as warm-start for last element
            obj.window_mpcstates = [obj.window_mpcstates(:,2:end), obj.Ak*obj.window_mpcstates(:,end)];

            %assume last element no control disturbance for warm-start
            [dim,n] = size(obj.window_controlDisturbances);
            obj.window_controlDisturbances = [obj.window_controlDisturbances(:,2:end), zeros(dim,1)];

            %assume no measError on last element for warm-start
            [dim,n] = size(obj.window_measError);
            obj.window_measError = [obj.window_measError(:,2:end), zeros(dim,1)];

            %push mhecontrols with new control input
            %assume no new thrust at tail end of mpc
            [dim,n] = size(obj.window_mpccontrols);
            obj.window_mhecontrols = [obj.window_mhecontrols(:,2:end), control];
            obj.window_mpccontrols = [obj.window_mpccontrols(:,2:end), zeros(dim,1)];

            %push new measurement into mhe measurements
            obj.window_measurements = [obj.window_measurements(:,2:end), measurement];
        end
        function obj = optimize(obj)
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
                -----------
                Saves results of optimization into class instance.
                Look at getOptimizeResult() for full list of variables that optimization gets.
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
            maxIter = 100;
            saveIter = -1;

            %TODO: parse the status, iter, and time for anything useful to know (throw errors if necessary)
            [status, iter, time] = solve(obj.opt, mu0, int32(maxIter), int32(saveIter));
            [Jcost, mpcUs, mheDs, x0s, xs, vs] = getOutputs(obj.opt); %get results

            mheXs = xs(:,1:obj.backwardT);
            mpcXs = xs(:,obj.backwardT+1:end);


            %save everything into mpcmhe
            obj.window_controlDisturbances = mheDs;
            obj.window_measError   = vs;
            obj.window_mhestates   = mheXs;
            obj.window_mpcstates   = mpcXs;
            obj.window_mpccontrols = mpcUs;
        end
        function [mheXs, mheDs, mheVs, mpcXs, mpcUs] = getOptimizeResult(obj)
            %{
                Description:
                ------------
                Returns everything optimize() gets you. This function exists
                as a design detail where we want to keep everything inside of th
                class object. Then we can access the internals of the class through getter
                methods. This keeps the code concise and abstracted.

                Parameters:
                -----------
                N/A

                Returns:
                --------
                @return mheXs: states estimated by MHE [N, backwardT]
                @return mheDs: control input error estimated by MHE [N, backwardT]
                @return mheVs: measurement error estimated by MHE [C, backwardT]
                @return mpcXs: states for optimal path created by MPC [N, forwardT]
                @return mpcUs: control inputs used for optimal path created by MPC [M, forwardT]

                N: dimension of state
                M: dimension of control input
                C: dimension of sensor measurement

                Access Notes:
                -------------
                Only call this after calling optimize(). Otherwise, you will not get
                results stored in the window objects and what this returns will pretty
                much be useless.
            %}
            mheXs = obj.window_mhestates;
            mheDs = obj.window_controlDisturbances;
            mheVs = obj.window_measError;
            mpcXs = obj.window_mpcstates;
            mpcUs = obj.window_mpccontrols;
        end
        function u = getMPCControl(obj)
            %{
                Description:
                ------------
                Getter method to access MPC optimized contorl input

                Parameters:
                -----------
                N/A

                Access Notes:
                -------------
                Only call this after calling optimize(). Otherwise, you will not get
                results stored in the window objects and what this returns will pretty
                much be useless.
            %}
            u = obj.window_mpccontrols(:,1);
        end
        %{
            TODO: implement a runController(meas,control) that outputs current estimated state and controls

            for benchmarking, it's goiing to be useful to perform some type of docktyping/polymorphism
            every class we want to benchmark MUST follow this format which will allow for easier plots.
        %}
    end
end