function nothing = mpcmhe_main(mode)
    delete("tmpC*") %delete any temporary files created
    close all
    clc
    %fixes randomness for reproducibility of any plots
    rng(1);


    %Keep 2d = false
    %Keep attitude = true
    use2D = false;
    useNonlinear = true;
    useAttitude = false;
    mpc_horizon = 20; %NOTE:MUST STAY IMBALANCED. Without this imbalance, the cost will get completely cancelled out and no optimization
    mhe_horizon = 10;
    total_time = 100; %total time in seconds
    tstep = 1;        %each time step is 1 second

    if mode == 1
        %time varying additive disturbance
        disturbance = [0.03;0.03;0.03];

        disturbance_fn = @(u) disturbance + u;
        mpcmheType = MpcMheThruster;

        mpcmheType = mpcmheType.setupCostGains(1e3,10,1,1);
        mpcmheType = mpcmheType.setupOptimizeGains(0.05,0.05,0.1,1000);
    elseif mode == 2
        %time invariant additive disturbance
        disturbance = [0.03;0.03;0.03];
        disturbance_fn = @(u) disturbance + u;
        mpcmheType = MpcMheNaive;
        mpcmheType = mpcmheType.setupCostGains(1e3,10,1,1);
        mpcmheType = mpcmheType.setupOptimizeGains(0.05,0.05,0.1,1000);
    elseif mode == 3
        %Time invariant matrix transform disturbance

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
        mpcmheType = MpcMheTransform;
        mpcmheType = mpcmheType.setupCostGains(1e3,10,1,1);
        mpcmheType = mpcmheType.setupOptimizeGains(0.07,1.3,0.1,1000);
    else
        %Time invariant matrix transform and additive disturbance on control vector
    end

    benchmark = ThrusterBenchmark;
    benchmark = benchmark.init(useNonlinear, mhe_horizon, mpc_horizon, total_time, tstep, use2D, useAttitude);
    %{
        The noiseQ, noiseR, and traj0 dimensions have to agree with the use2D and useAttitude
    %}

    % 2d no att
    % noiseQ = @() [0;0;0;0];
    % noiseR = @() [0;0;0;0];
    % traj0 = [1;1;0.001;0.001];

    % 3d w/ att
    % noiseQ = @() [0;0;0;0;0;0];
    % noiseR = @() [0;0;0;0;0;0];
    % traj0 = [1;1;1;1;1;1;0;0;0;0;0;0];

    %3d no att
    noiseQ = @() [0;0;0;0;0;0];
    noiseR = @() [0;0;0;0;0;0];
    traj0 = [10;10;10;0;0;0];

    benchmark = benchmark.runBenchmark(traj0, noiseQ, noiseR, disturbance_fn,mpcmheType);

    disp("Actual transform: ")
    disp(disturbance)
    delete("tmpC*") %delete any temporary files created
end