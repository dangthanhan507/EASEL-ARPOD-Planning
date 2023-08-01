% ==============================================================
% thruster_env.m
% 
% This file is interested in creating a specialized environment
% meant running realtime simulation and setting up the code
% running as a benchmark. It will save all of its data
% after running which we can use to save stats and benchmark.
% This should support running with attitude, 2d, and 3d.
% 
% Author: An Dang
%
% ==============================================================


classdef ThrusterBenchmark
    properties
        %members of ThrusterBenchmark
        init_traj

        %algorithm parameters
        mhe_horizon
        mpc_horizon
        
        %simulation parameters
        total_time
        tstep
        
        %simulation options
        use_nonlinear
        use_attitude
        use_2d


        %benchmark save data
        true_trajs
        estimated_trajs
        measurements
        control_vectors

        true_att_trajs
        est_att_trajs
        att_controls

        %TODO: add its usage after discussing with copp
        disturbances
    end
    methods 
        function obj = init(obj, useNonlinear, mhe_horizon, mpc_horizon, total_time, tstep, use2D, useAttitude)
            %{
                TODO: fill this docstring
                Description:
                ------------


                Parameters:
                -----------
                @params useNonlinear
                @params mhe_horizon
                @params mpc_horizon
                @params total_time
                @params tstep
                @params use2D
                @params useAttitude

                Returns:
                --------
                Class instance with the parameters we want to start benchmark.
            %}

            %set simulation options
            obj.use_nonlinear = useNonlinear;
            obj.use_attitude = useAttitude;
            obj.use_2d = use2D;

            %set simulation parameters
            obj.total_time = total_time;
            obj.tstep = tstep;

            %set algorithm parameters
            obj.mhe_horizon = mhe_horizon;
            obj.mpc_horizon = mpc_horizon;

        end
        function [posTraj, attTraj] = decomposeTraj(obj, traj)
            if obj.use_attitude
                if obj.use_2d
                    % [x,y,yaw]
                    % 1->2 pos, 3 att, 4->5 posdot, 6 attdot
                    posTraj = [traj(1:2,:);traj(4:5,:)];
                    attTraj = [traj(3,:);traj(6,:)];
                else
                    % [x,y,z,roll,pitch,yaw]
                    % 1->3 pos, 4->6 att, 7->9 posdot, 10->12 attdot
                    posTraj = [traj(1:3,:);traj(7:9,:)];
                    attTraj = [traj(4:6,:);traj(10:12,:)];
                end
            else
                posTraj = traj;
                attTraj = nargin;
            end
        end
        function obj = runBenchmark(obj, traj0, noiseQ, noiseR)
            %{
                TODO: fill this docstring
                Description:
                ------------


                Parameters:
                -----------
                @params traj0
                @params noiseQ
                @params noiseR


                Returns:
                --------

                NOTES:
                ---------
                traj0 should contain both position and attitude and velocity of both too
            %}


            % TODO: make this check work with attitude
            % if obj.use_2d
            %     [dim,n] = size(traj0);
            %     if dim ~= 4
            %         error("Non 2D Trajectory used");
            %     end
            % end
            [traj0,att0] = obj.decomposeTraj(traj0);

            %MISSION SETUP
            %=====================
            %hard-coded thruster type
            Mission = ARPOD_Mission;
            Mission = Mission.initMission(traj0, obj.tstep, obj.use_2d, 1);
            Mission = Mission.useAttitude(att0);

            [A,B] = ARPOD_Dynamics.linearHCWDynamics(obj.tstep, Mission.mu, Mission.a, obj.use_2d);
            [Aatt,Batt] = ARPOD_Dynamics.attitudeLVLH(obj.tstep, obj.use_2d);
            %=====================

            
            %MPC-MHE setup
            %=====================
            %plug in state dim for measurement dim (we assume full state measurements)
            mpcmhe = MpcMheThruster;

            if obj.use_2d
                control_dim = 2;
                attcontrol_dim = 1;
            else
                control_dim = 3;
                attcontrol_dim = 3;
            end
            [state_dim,n] = size(traj0);
            [attstate_dim,n] = size(att0);
            mpcmhe = mpcmhe.init(traj0, obj.mhe_horizon, obj.mpc_horizon, A, B, state_dim, control_dim);

            if obj.use_attitude
                mpcmhe = mpcmhe.initAtt(att0, Aatt, Batt, attcontrol_dim);
            end
            mpcmhe = mpcmhe.setupOptimize(0.01, 0.01, 0.1, 200);

            %setup MPC-MHE windows
            window_states = zeros(state_dim,obj.mhe_horizon);
            window_measurements = zeros(state_dim,obj.mhe_horizon);
            window_controls = zeros(control_dim,obj.mhe_horizon);


            window_attstates = zeros(attstate_dim, obj.mhe_horizon);
            window_attmeasurements = zeros(attstate_dim, obj.mhe_horizon);
            window_attcontrols = zeros(attcontrol_dim, obj.mhe_horizon);
            %=====================

            %Fill up MPC-MHE windows with measurements
            %==========================================
            u = zeros(control_dim,1);
            uatt = zeros(attcontrol_dim,1);

            noiseQatt = @() zeros(length(att0),1);
            noiseRatt = @() zeros(length(att0),1);
            for i = obj.tstep:obj.tstep:obj.mhe_horizon % [tstep -> obj.mhe_horizon] inclusive steps by tstep
                Mission = Mission.nextStep(u,noiseQ,noiseR,true);
                if obj.use_attitude
                    Mission = Mission.nextAttStep(uatt, noiseQatt, noiseRatt);
                end

                true_traj = Mission.traj;
                meas = Mission.sensor;
                window_measurements(:,i/obj.tstep) = meas;
                window_states(:,i/obj.tstep) = true_traj;
                window_controls(:,i/obj.tstep) = u;

                if obj.use_attitude
                    window_attstates(:,i/obj.tstep) = Mission.att;
                    window_attmeasurements(:,i/obj.tstep) = Mission.att_sensor;
                    window_attcontrols(:,i/obj.tstep) = uatt;
                end
            end
            mpcmhe = mpcmhe.setupWindows(window_states, window_measurements, window_controls);

            if obj.use_attitude
                mpcmhe = mpcmhe.setupAttWindows(window_attstates, window_attmeasurements, window_attcontrols);
            end
            mpcmhe = mpcmhe.optimize();

            u = mpcmhe.getMPCControl();

            uatt = mpcmhe.getMPCAttControl();
            %==========================================

            %Run MPC-MHE and benchmark now
            %==============================
            [n, num_steps] = size((obj.tstep+obj.mhe_horizon):obj.tstep:obj.total_time);

            obj.true_trajs      = zeros(state_dim, num_steps);
            obj.estimated_trajs = zeros(state_dim, num_steps);
            obj.measurements    = zeros(state_dim,  num_steps);
            obj.control_vectors = zeros(control_dim, num_steps);

            obj.true_att_trajs  = zeros(attstate_dim, num_steps);
            obj.est_att_trajs   = zeros(attstate_dim, num_steps);
            obj.att_controls    = zeros(attcontrol_dim, num_steps);
            for i = (obj.tstep+obj.mhe_horizon):obj.tstep:obj.total_time
                idx = (i - obj.mhe_horizon)/obj.tstep;

                Mission = Mission.nextStep(u,noiseQ,noiseR,true);
                if obj.use_attitude
                    Mission = Mission.nextAttStep(uatt,noiseQatt,noiseRatt);      

                    att_meas = Mission.att_sensor;
                    mpcmhe = mpcmhe.shiftAttitudeWindows(att_meas,uatt);
                end
                true_traj = Mission.traj;
                meas = Mission.sensor;

                mpcmhe = mpcmhe.shiftWindows(meas, u);
                mpcmhe = mpcmhe.optimize();
                u = mpcmhe.getMPCControl();
                est_traj = mpcmhe.getMHEState();

                if obj.use_attitude                    
                    uatt = mpcmhe.getMPCAttControl();
                    est_att_traj = mpcmhe.getMHEAttState();

                    obj.true_att_trajs(:,idx) = Mission.att;
                    obj.est_att_trajs(:,idx) = est_att_traj;

                    obj.att_controls(:,idx) = uatt;
                end



                
                obj.true_trajs(:,idx) = true_traj;
                obj.estimated_trajs(:,idx) = est_traj;
                obj.measurements(:,idx) = meas;
                obj.control_vectors(:,idx) = u;
            end
            %==============================

            %return obj
        end
    end
end