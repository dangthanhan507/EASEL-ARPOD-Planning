function benchmark = mpcmhe_main(mode)
    delete("tmpC*") %delete any temporary files created
    close all
    clc
    %fixes randomness for reproducibility of any plots
    rng(1);


    %Keep 2d = false
    %Keep attitude = true
    use2D = false;
    useNonlinear = true;
    useAttitude = true;
    mpc_horizon = 10; %NOTE:MUST STAY IMBALANCED. Without this imbalance, the cost will get completely cancelled out and no optimization
    mhe_horizon = 5;
    total_time = 30+mhe_horizon; %total time in seconds + setup time
    tstep = 1;        %each time step is 1 second

    state_dim = 6;
    meas_dim = 6;
    control_dim = 3;
    att_dim = 6;
    att_control_dim = 3;
    
    mpcmhe = MPCMHE;
    mpcmhe = mpcmhe.init(mode, mhe_horizon, mpc_horizon, state_dim, meas_dim, control_dim, att_dim);
    [A,B] = ARPOD_Dynamics.linearHCWDynamics(tstep, ARPOD_Mission.mu, ARPOD_Mission.a, use2D);
    [Aatt,Batt] = ARPOD_Dynamics.attitudeLVLH(tstep, use2D);
    mpcmhe = mpcmhe.setupLTILinearFixedBody(A,B,eye(state_dim));
    mpcmhe = mpcmhe.setupAttitudeConstraints(Aatt, Batt, eye(att_dim));
    
    if mode == 1
        disturbance = [0.03;0.03;0.03];
        disturbance_fn = @(u) disturbance + u;

        mpcmhe = mpcmhe.setupSimpleObjective(1e3,10,1,1);
        mpcmhe = mpcmhe.setupAttitudeCost(1e3,10,1,1);
        mpcmhe = mpcmhe.setupVariableLimits(1.3,0.1,0.2);
    elseif mode == 2
        disturbance = [0.03;0.03;0.03];
        disturbance_fn = @(u) disturbance + u;

        mpcmhe = mpcmhe.setupSimpleObjective(1e3,10,1,1);
        mpcmhe = mpcmhe.setupAttitudeCost(1e3,10,1,1);
        mpcmhe = mpcmhe.setupVariableLimits(1.3,0.1,0.2);
    elseif mode == 3
        disturbance = eye(3);

        % a = 0;
        % b = pi/3;
        % c = -pi/3;
        % Rx = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)];
        % Ry = [cos(b) 0 sin(b); 0 1 0; -sin(b) 0 cos(b)];
        % Rz = [cos(c) -sin(c) 0; sin(c) cos(c) 0; 0 0 1];
        % disturbance = Rz*Ry*Rx;

        % disturbance = [0,1,0;
        %             0,0,1;
        %             1,0,0];
        disturbance_fn = @(u) disturbance*u;

        mpcmhe = mpcmhe.setupSimpleObjective(1e3,10,1,1);
        mpcmhe = mpcmhe.setupAttitudeCost(1e3,10,1,1);
        mpcmhe = mpcmhe.setupVariableLimits(1.3,0.1,0.2);
    else
        %Time invariant matrix transform and additive disturbance on control vector
        disturbanceD = eye(3);
        disturbanced = zeros(3,1);
        disturbance_fn = @(u) disturbanceD*u + disturbanced;
        mpcmhe = mpcmhe.setupSimpleObjective(1e3,10,1,1);
        mpcmhe = mpcmhe.setupAttitudeCost(1e3,10,1,1);
        mpcmhe = mpcmhe.setupVariableLimits(1.3,0.1,0.2);
    end

    %TODO: remove 2d from code.
    mpcmhe = mpcmhe.createOptimization();
    mpcmhe = mpcmhe.addAttitudeOptimization();
    mpcmhe = mpcmhe.setupOptimizationVar();


    benchmark = ThrusterBenchmark;
    benchmark = benchmark.init(useNonlinear, mhe_horizon, mpc_horizon, total_time, tstep, use2D, useAttitude);
    %{
        The noiseQ, noiseR, and traj0 dimensions have to agree with the use2D and useAttitude
    %}

    %======================= NOISE ===========================
    % 2d no att
    % noiseQ = @() [0;0;0;0];
    % noiseR = @() [0;0;0;0];
    % traj0 = [1;1;0.001;0.001];

    % 3d w/ att
    noiseQ = @() [0;0;0;0;0;0];
    noiseR = @() [0;0;0;0;0;0];
    traj0 = [1;1;1;1;1;1;0;0;0;0;0;0];
    %======================= NOISE END ========================

    benchmark = benchmark.runBenchmark(traj0, noiseQ, noiseR, disturbance_fn, mpcmhe);

    disp("Actual transform: ")
    disp(disturbance)
    delete("tmpC*") %delete any temporary files created
end