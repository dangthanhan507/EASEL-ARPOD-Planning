classdef EKF
    %{
        This is the discrete HCW EKF
        Prediction is based on discrete matrices Ax+Bu
        Measurement is based on using jacobian with nonlinear functions.
        x_t+1 = Ax_t + Bu_t
        c_t = h(x_t)
        Example Using EKF:
        ------------------
            Call state_estimator = ChaserEKF
            Q = eye(6)
            R = eye(3)
            state0 = zeros(6,1);
            cov0 = zeros(6,6);
            u = ....;
            tstep = 1;
            state_estimator.initEKF(state0, cov0);
            state_estimator.prediction(u, Q, tstep)
            state_estimator.sensor_correct(z_t, measCov, phase)
            
    %}
    properties
        %{
            These are properties of EKF tracked over time.
            The only ones that matter.
        %}
        A
        B
        state = [0;0;0;0;0;0];
        cov = eye(6);

        Q % system Covariance
        R % measure Covariance
    end
    methods 
        function obj = init(obj, state0, cov0, A, B, systemCov, measCov)
            %{
                Parameters:
                ------------
                    state0: 6x1 vector that defines state of spacecraft.
                    cov0: 6x6 matrix that defines current covariance
                Description:
                ------------
                    Simply initializing object for state estimation.
            %}
            obj.state = state0;
            obj.cov = cov0;

            obj.A = A;
            obj.B = B;
            obj.Q = systemCov;
            obj.R = measCov;
        end
        function obj = prediction(obj, u)
            %{
                Parameters:
                ------------
                    u: 3x1 vector control input of thrusters (in accel)
                    systemCov: 6x6 matrix system covariance of EKF.
                    
                Description:
                ------------
                    Implementation of the prediction of spacecraft state
                    given the previous state. This is a rough estimate that
                    uses linear dynamics of the system. Namely, the
                    Hill-Clohessy-Wiltshire equations for relative space
                    motion. This assumes discrete control input over time.
            %}
            obj.state = obj.A*obj.state + obj.B*u;
            obj.cov = obj.A*obj.cov*transpose(obj.A) + obj.Q;
        end
        function obj = sensor_correctEKF(obj, z_t)
            %{
                Parameters:
                -----------
                    z_t: 3x1 or 2x1 measurement vector of EKF.
                    measCov: 3x3 or 2x2 covariance matrix of meas.
                    phase: the phase that spacecraft is in
                Description:
                ------------
                    NOTE: requires prediction to have been run beforehand.
                    Runs the sensor_correct that corrects prediction to a
                    much more accurate state by conditioning on the
                    observed measurement variable in the HMM formulation of
                    EKF.
            %}

            %take jacobian of the measurements based on the state
            H = utils.measureJacobian(obj.state);
            %calculate riccati K gain
            K_gain = obj.cov*transpose(H)*inv(H*obj.cov*transpose(H)+obj.R);

            %convert predicted state into a measurement to compare with the
            %real measurement, then correct covariance and state.
            measure = utils.measure(obj.state);
            obj.state = obj.state + K_gain*(z_t-measure);
            obj.cov = (eye(6) - K_gain*H)*obj.cov;
        end
        function obj = sensor_correctKF(obj, z_t, C)
            H = C;
            measure = C*obj.state;
            K_gain = obj.cov*transpose(H)*inv(H*obj.cov*transpose(H) + obj.R);

            obj.state = obj.state + K_gain*(z_t-measure);
            obj.cov = (eye(6) - K_gain*H)*obj.cov;
        end
        function obj = estimateEKF(obj, u, z_t)
            %{
                Combines the use of both the predict and correct step
                For the sake of generalizing state estimators.
                NOTE: I hope ducktyping works :)
            %}
            obj = obj.prediction(u);
            obj = obj.sensor_correctEKF(z_t);
        end
        function obj = estimateKF(obj, u, z_t, C)
            obj = obj.prediction(u);
            obj = obj.sensor_correctKF(z_t, C);
        end

        function state = outputResult(obj)
            state = obj.state;
        end
    end
end