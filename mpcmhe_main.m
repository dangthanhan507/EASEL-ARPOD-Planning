close all
clear
clc
rng(1);

Mission = ARPOD_Mission;

traj = [-0.2;-0.2;0.2;0.001;0.001;0.001];
total_time = 100;
tstep = 1; % update every second
Mission.initMission(traj, tstep);

mpc_horizon = 30;
mhe_horizon = 10;

mpc_Q = [1,0,0,0,0,0;
        0,1,0,0,0,0;
        0,0,1,0,0,0;
        0,0,0,100,0,0;
        0,0,0,0,100,0;
        0,0,0,0,0,100];
mpc_R = 100*eye(3);

mhe_Q = 1e20*eye(6);
mhe_R = 1e3*eye(3);

stats = MPCMHE_stats;

%setup code before reaching horizon
window_states = zeros(6,mhe_horizon);
window_measurements = zeros(3,mhe_horizon);
window_controls = zeros(3,mhe_horizon);
window_phases = zeros(mhe_horizon);


Mission = Mission.initMission(traj, tstep);
phase = Mission.phase;

mpcmhe = MPCMHE_thruster;
mpcmhe = mpcmhe.setupTensCalc();
mpcmhe = mpcmhe.init(traj, mhe_horizon, mpc_horizon, tstep);
%mpcmhe = mpcmhe.tensCalcOpt(mpc_Q,mpc_R,mhe_Q,mhe_R,tstep, 1, 1, ARPOD_Mission.ubar);

for i = tstep:tstep:total_time
    if phase == 1
        noiseQ = @() [0;0;0;0;0;0];
        noiseR = @() mvnrnd([0;0;0], [0.001, 0.001, 0.01]).';

        ekfQ = 1e-20*eye(6);
        ekfR = diag([0.001,0.001,0.01]);
    elseif phase == 2
        noiseQ = @() [0;0;0;0;0;0];
        noiseR = @() 10*mvnrnd([0;0;0], [0.001, 0.001, 0.01]).';

        ekfQ = 1e-20*eye(6);
        ekfR = diag([0.001,0.001,0.01]);
    else
        noiseQ = @() [0;0;0;0;0;0];
        noiseR = @() mvnrnd([0;0;0], [0.001, 0.001, 1e-5]).';

        ekfQ = 1e-20*eye(6);
        ekfR = diag([0.001,0.001,1e-5]);
    end

    if i/tstep < mhe_horizon
        % setup the control input as 0
        u = [0;0;0];
    else

    end
    Mission = Mission.nextStep(u,noiseQ, noiseR);

    meas = Mission.sensor;
    traj = Mission.traj;
    phase = Mission.phase;

    if i/tstep <= mhe_horizon
        %setup measurements
        window_measurements(:,i/tstep) = MPCMHE_thruster.senseModify(meas);
        window_controls(:,i/tstep) = u;
        window_states(:,i/tstep) = chaserEKF.state;
        window_phases(i/tstep) = phase;

        if i/tstep == mhe_horizon
            mpcmhe = mpcmhe.initWindows(window_controls, 100*ones(3,mpc_horizon), window_measurements, window_states, 100*ones(6,mpc_horizon), window_phases);
            mpcmhe = mpcmhe.tensCalcOpt(mpc_Q,mpc_R,mhe_Q,mhe_R,tstep, 1, 1, ARPOD_Mission.ubar);
        end
    else
        error("oof");
        
    end
        
end