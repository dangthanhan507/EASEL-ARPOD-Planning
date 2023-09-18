function benchmark = mpcmhe_6dof_main()
    % delete("tmpC*") %delete any temporary files created
    close all
    clc
    %fixes randomness for reproducibility of any plots
    rng(1);


    mode = 1;
    %Keep 2d = false
    %Keep attitude = true
    use2D = false;
    useNonlinear = true;
    useAttitude = true;
    mpc_horizon = 10; %NOTE:MUST STAY IMBALANCED. Without this imbalance, the cost will get completely cancelled out and no optimization
    mhe_horizon = 5;
    total_time = 10+mhe_horizon; %total time in seconds + setup time
    tstep = 1;        %each time step is 1 second

    state_dim = 6;
    control_dim = 6;
    att_dim = 6;
    
    if useNonlinear == true
        meas_dim = 3;
    else
        meas_dim = 6;
    end

    mpcmhe = MPCMHE_6dof;
    mpcmhe = mpcmhe.init(mode, mhe_horizon, mpc_horizon, state_dim, meas_dim, control_dim, att_dim);
    [A,B] = ARPOD_Dynamics.linearHCWDynamics(tstep, ARPOD_Mission.mu, ARPOD_Mission.a, use2D);
    [Aatt,Batt] = ARPOD_Dynamics.attitudeLVLH(tstep, use2D);
    
    
    % d_disturbance = 0.01*ones(6,1);
    d_disturbance = zeros(6,1);
    D_disturbance = eye(3);
    bigD = blkdiag(D_disturbance,D_disturbance);
    disturbance_fn = @(u) bigD* (u + d_disturbance);

    mpcmhe = mpcmhe.setupVariableLimits(0.05,0.1,0.3);
    mpcmhe = mpcmhe.setupUncoupledUnconstrainedNonlinear(A,B,1e10,1e5,1e3,1e3);
    mpcmhe = mpcmhe.setupAttitudeConstraints(Aatt, Batt, eye(att_dim));
    mpcmhe = mpcmhe.setupAttitudeCost(1e10,1,1,1);

    mpcmhe = mpcmhe.addAttitudeOptimization();
    mpcmhe = mpcmhe.setupOptimizationVar();


    benchmark = ThrusterBenchmark;
    benchmark = benchmark.init(useNonlinear, mhe_horizon, mpc_horizon, total_time, tstep, true, useAttitude);
    %{
        The noiseQ, noiseR, and traj0 dimensions have to agree with the use2D and useAttitude
    %}
    %======================= NOISE ===========================

    % 3d w/ att
    noiseQ = @() [0;0;0;0;0;0];
    % noiseR = @() [0;0;0;0;0;0];
    if useNonlinear == true
        noiseR = @() [0;0;0];
    else
        noiseR = @() [0;0;0;0;0;0];
    end
    

    % traj0 = [1;1;1;1;1;1;0;0;0;0;0;0];
    traj0 = [1;1;1;0;0;0;0;0;0;0;0;0];
    %======================= NOISE END ========================

    disp("passed compiliation")
    benchmark = benchmark.runBenchmark(traj0, noiseQ, noiseR, disturbance_fn, mpcmhe);

    % disp("Actual transform: ")
    % delete("tmpC*") %delete any temporary files created
end