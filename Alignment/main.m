function benchmark = main(algo_type)
    clc
    close all
    rng(1);

    delete("tmpC*")
    if exist("@tmpC/", 'dir')
        rmdir("@tmpC", 's')
    end
    time_step = 1; % 1 second
    sim_time  = 40;


    x0   = [10;10;10];
    att0 = [1;1;1];
    xdot0 = [0;0;0];
    attdot0 = [0;0;0];

    process_noise = [0,0,0,0,0,0];
    noiseQ        = @() mvnrnd([0;0;0;0;0;0], [0,0,0,0,0,0] + process_noise).';

    sensor_noise = [1e-3,1e-3,1e-5];
    noiseR        = @() mvnrnd([0;0;0], [0,0,0] + sensor_noise);


    d = [0;0;0;0;0;0];
    D = eye(6);
    % disturbance_fn = @(u) D*u + d;
    disturbance_fn = @(u) u;
    % apply_fn = @(u) u(1:3,:) - u(4:6,:);
    apply_fn = @(u) u;
    
    forwardHorizon  = 10;
    backwardHorizon = 5;

    xTraj = [x0;xdot0];
    attTraj = [att0;attdot0];

    if algo_type == 1
        algo = mpcmhe;
        algo = algo.init(backwardHorizon, forwardHorizon, xTraj, attTraj);
        algo = algo.setupOptimizationCode(time_step);
    else
        algo = mpcekf;
        algo = algo.init;
    end

    disp("Starting Simulation")
    % control    = zeros(6,1);
    control = zeros(3,1);
    controlAtt = zeros(3,1);
    

    ctr = 1;
    window_mheXsWarmUp = zeros(6, backwardHorizon);

    for i = time_step:time_step:sim_time

        %%% DO NOT TOUCH SIMULATION CODE %%%
        %%%%% SIMULATION ONLY %%%%%
        % xTraj = utils.nonlinearMotionSolver(xTraj, utils.rotRPY(attTraj)*apply_fn(disturbance_fn(control)), time_step);
        xTraj = utils.nonlinearMotionSolver(xTraj, apply_fn(disturbance_fn(control)), time_step);
        attTraj = utils.attitudeSolver(attTraj, controlAtt, time_step);

        % meas    = utils.measure(xTraj);
        meas = xTraj;
        measAtt = attTraj;
        %%%%% SIMULATION ONLY END %%%%%

        if ctr <= backwardHorizon
            window_mheXsWarmUp(:,ctr) = xTraj;
            if ctr == backwardHorizon
                algo.window_mheXs = window_mheXsWarmUp;
            end
        end
        ctr = ctr + 1;

        algo = algo.estimate_and_control(control, controlAtt, meas, measAtt);
        [control, controlAtt, state, stateAtt] = algo.outputResults();

        if i >= backwardHorizon
            disp("State Estimated")
            disp(state)
            disp("True State")
            disp(xTraj)
            disp("Meas Estimated")
            disp(utils.measure(state))
            disp("True Meas")
            disp(meas)
            disp("Control")
            disp(control)
        end
    end
    disp("Finished")
    disp("Final Traj")
    disp(xTraj)
    disp("Final Att")
    disp(attTraj)

    disp("Estimated Traj")
    disp(state)


    benchmark = algo;
end