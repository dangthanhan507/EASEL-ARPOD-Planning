%{
    mpcmhe_main.m
    ===============
        -> mpcmhe algorithm runner that performs estimation and planning simultaneously
        -> TODO: fix this up
            -> remove arpod sensor (we can assume good state)
            -> get this to solve in closed loop format

    NOTES:
    ========
        -> stay within phase 2



    
    IMPORT FILES:
    =============
        -> mpcmhe_lib
%}
delete("tmpC*") %delete any temporary files created
close all
clear
clc

rng(1);

Mission = ARPOD_Mission;

traj = [1;1;-1;0.001;0.001;0.001];
total_time = 100;
tstep = 1; % update every second

mpc_horizon = 30;
mhe_horizon = 10;

%setup code before reaching horizon
window_states = zeros(6,mhe_horizon);
window_measurements = zeros(6,mhe_horizon);
window_controls = zeros(3,mhe_horizon);
window_phases = zeros(mhe_horizon);

Mission = Mission.initMission(traj, tstep);
phase = Mission.phase;

%noise functions...
noiseQ = @() [0;0;0;0;0;0];
noiseR = @() [0;0;0;0;0;0];


mpcmhe = MpcMheThruster;
mpcmhe = mpcmhe.init(traj,mhe_horizon,mpc_horizon, tstep, 6, 3);

u = [0.0;0.0;0.0];
for i = tstep:tstep:total_time

    %control input

    Mission = Mission.nextStepMod(u, noiseQ,noiseR);


    if i/tstep <= mhe_horizon
        %this is setup code when we still need to fill up data for the MHE horizon
        window_measurements(:,i/tstep) = Mission.sensor;
        window_states(:,i/tstep) = Mission.traj;
        window_controls(:,i/tstep) = u;

        if i/tstep == mhe_horizon
            %finally filled up mhe horizon with datapoints.
            %perform optimization and get our first real MPC output
            mpcmhe = mpcmhe.setupWindows(window_states, window_measurements, window_controls);
            mpcmhe = mpcmhe.setupOptimize(0.01,0.01,0.1,200);
            mpcmhe = mpcmhe.optimize();
            u = mpcmhe.getMPCControl();
        end
    else
        %needs window shift
        meas = Mission.sensor;
        control = u;
        mpcmhe = mpcmhe.shiftWindows(meas,control);
        mpcmhe = mpcmhe.optimize();
        u = mpcmhe.getMPCControl();
        disp(Mission.traj);
        disp(u);
    end
end