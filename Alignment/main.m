function benchmark = main()
    clc
    close all
    rng(1);

    delete("tmpC*")
    if exist("@tmpC/", 'dir')
        rmdir("@tmpC", 's')
    end
    time_step = 1; % 1 second
    sim_time  = 40;


    x0   = [1;1;1];
    att0 = [1;1;1];
    xdot0 = [0;0;0];
    attdot0 = [0;0;0];

    noiseQ    = [0;0;0;0;0;0];
    noiseQAtt = [0;0;0;0;0;0];

    noiseR    = [0;0;0];
    noiseRAtt = [0;0;0;0;0;0];


    d = [0;0;0;0;0;0];
    D = eye(6);
    disturbance_fn = @(u) D*u + d;
    apply_fn = @(u) u(1:3,:) - u(4:6,:);
    
    backwardHorizon = 15;
    forwardHorizon  = 20;

    xTraj = [x0;xdot0];
    attTraj = [att0;attdot0];

    algo = mpcmhe;
    algo = algo.init(backwardHorizon, forwardHorizon, xTraj, attTraj);
    algo = algo.setupOptimizationCode(time_step);

    disp("Starting Simulation")
    control    = zeros(6,1);
    controlAtt = zeros(3,1);


    for i = time_step:time_step:sim_time

        %%% DO NOT TOUCH SIMULATION CODE %%%
        %%%%% SIMULATION ONLY %%%%%
        xTraj = utils.nonlinearMotionSolver(xTraj, utils.rotRPY(attTraj)*apply_fn(disturbance_fn(control)), time_step);
        % xTraj = utils.nonlinearMotionSolver(xTraj, apply_fn(disturbance_fn(control)), time_step);
        attTraj = utils.attitudeSolver(attTraj, controlAtt, time_step);

        meas    = utils.measure(xTraj);
        measAtt = attTraj;
        %%%%% SIMULATION ONLY END %%%%%

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
        end
    end
    disp("Finished")
    disp("Final Traj")
    disp(xTraj)
    disp("Final Att")
    disp(attTraj)

    disp("Estimated Traj")
    disp(state)


    benchmark = 1;
end