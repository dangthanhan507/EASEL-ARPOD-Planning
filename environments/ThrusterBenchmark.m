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


        %benchmark save data
        true_trajs
        estimated_trajs
        measurements
        control_vectors
        true_control

        true_att_trajs
        est_att_trajs
        att_controls
        att_measurements


        %save output of the controller per timestep
        mheXs
        mheDs
        mheVs
        mpcXs
        mpcUs

        mheXsAtt
        mheDsAtt
        mheVsAtt
        mpcXsAtt
        mpcUsAtt

        thruster6dof

    end
    methods 
        function obj = init(obj, useNonlinear, mhe_horizon, mpc_horizon, total_time, tstep, thruster6dof, useAttitude)
            %{
                Description:
                ------------


                Parameters:
                -----------
                @params useNonlinear
                @params mhe_horizon
                @params mpc_horizon
                @params total_time
                @params tstep
                @params useAttitude

                Returns:
                --------
                Class instance with the parameters we want to start benchmark.
            %}
            obj.thruster6dof = thruster6dof;

            %set simulation options
            obj.use_nonlinear = useNonlinear;
            obj.use_attitude = useAttitude;

            %set simulation parameters
            obj.total_time = total_time;
            obj.tstep = tstep;

            %set algorithm parameters
            obj.mhe_horizon = mhe_horizon;
            obj.mpc_horizon = mpc_horizon;

        end
        function [posTraj, attTraj] = decomposeTraj(obj, traj)
            if obj.use_attitude
                % [x,y,z,roll,pitch,yaw]
                % 1->3 pos, 4->6 att, 7->9 posdot, 10->12 attdot
                posTraj = [traj(1:3,:);traj(7:9,:)];
                attTraj = [traj(4:6,:);traj(10:12,:)];
            else
                posTraj = traj;
                attTraj = nargin;
            end
        end
        function obj = runBenchmark(obj, traj0, noiseQ, noiseR, disturbance_fn, mpcmhe)
            [traj0,att0] = obj.decomposeTraj(traj0);

            %MISSION SETUP
            %=====================
            %hard-coded thruster type
            Mission = ARPOD_Mission;
            Mission = Mission.initMission(traj0, obj.tstep, false, 1);
            Mission = Mission.useAttitude(att0);
            %=====================

            
            %MPC-MHE setup
            %=====================
            %plug in state dim for measurement dim (we assume full state measurements)
            if obj.thruster6dof
                control_dim = 6;
            else
                control_dim = 3;
            end
            attcontrol_dim = 3;
            [state_dim,n] = size(traj0);
            [attstate_dim,n] = size(att0);

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

            if obj.thruster6dof
                uMat = mpcmhe.uMat;
            else
                uMat = eye(3);
            end
            for i = obj.tstep:obj.tstep:obj.mhe_horizon % [tstep -> obj.mhe_horizon] inclusive steps by tstep
                Mission = Mission.nextStep(uMat*u,noiseQ,noiseR,true);
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
            mpcmhe = mpcmhe.setupWindows(traj0, window_states, window_measurements, window_controls);

            if obj.use_attitude
                mpcmhe = mpcmhe.setupAttWindows(att0, window_attstates, window_attmeasurements, window_attcontrols);
            end
            mpcmhe = mpcmhe.optimize();

            u = mpcmhe.getMPCControl();
            if obj.use_attitude
                uatt = mpcmhe.getMPCAttControl();
            end
            %==========================================

            %Run MPC-MHE and benchmark now
            %==============================
            [n, num_steps] = size((obj.tstep+obj.mhe_horizon):obj.tstep:obj.total_time);

            obj.true_trajs      = zeros(state_dim, num_steps);
            obj.estimated_trajs = zeros(state_dim, num_steps);
            obj.measurements    = zeros(state_dim,  num_steps);
            obj.control_vectors = zeros(control_dim, num_steps);
            obj.true_control    = zeros(control_dim, num_steps);

            obj.mheXs           = cell(num_steps);
            obj.mheDs           = cell(num_steps);
            obj.mheVs           = cell(num_steps);
            obj.mpcXs           = cell(num_steps);
            obj.mpcUs           = cell(num_steps);

            if obj.use_attitude
                obj.true_att_trajs  = zeros(attstate_dim, num_steps);
                obj.est_att_trajs   = zeros(attstate_dim, num_steps);
                obj.att_measurements        = zeros(attstate_dim, num_steps);
                obj.att_controls    = zeros(attcontrol_dim, num_steps);

                obj.mheXsAtt           = cell(num_steps);
                obj.mheDsAtt           = cell(num_steps);
                obj.mheVsAtt           = cell(num_steps);
                obj.mpcXsAtt           = cell(num_steps);
                obj.mpcUsAtt           = cell(num_steps);
            end
            for i = (obj.tstep+obj.mhe_horizon):obj.tstep:obj.total_time
                idx = (i - obj.mhe_horizon)/obj.tstep;

                Mission = Mission.nextStep(uMat*disturbance_fn(u),noiseQ,noiseR,true);
                Mission = Mission.nextAttStep(uatt,noiseQatt,noiseRatt);      

                att_meas = Mission.att_sensor;
                true_traj = Mission.traj;
                meas = Mission.sensor;

                disp(u);
                disp(meas);

                mpcmhe = mpcmhe.shift(meas, u, att_meas, uatt);
                mpcmhe = mpcmhe.optimize();
                u = mpcmhe.getMPCControl();
                est_traj = mpcmhe.getMHEState();


                [mheXs, mheDs, mheVs, mpcXs, mpcUs] = mpcmhe.getOptimizeResult();
                obj.mheXs{idx} = mheXs;
                obj.mheDs{idx} = mheDs;
                obj.mheVs{idx} = mheVs;
                obj.mpcXs{idx} = mpcXs;
                obj.mpcUs{idx} = mpcUs;
                if obj.use_attitude
                    [mheXs, mheDs, mheVs, mpcXs, mpcUs] = mpcmhe.getOptimizeResultAttitude();               
                    uatt = mpcmhe.getMPCAttControl();
                    est_att_traj = mpcmhe.getMHEAttState();

                    obj.true_att_trajs(:,idx) = Mission.att;
                    obj.est_att_trajs(:,idx) = est_att_traj;
                    obj.att_measurements(:,idx) = att_meas;
                    obj.att_controls(:,idx) = uatt;

                    obj.mheXsAtt{idx} = mheXs;
                    obj.mheDsAtt{idx} = mheDs;
                    obj.mheVsAtt{idx} = mheVs;
                    obj.mpcXsAtt{idx} = mpcXs;
                    obj.mpcUsAtt{idx} = mpcUs;
                end
                
                obj.true_trajs(:,idx) = true_traj;
                obj.estimated_trajs(:,idx) = est_traj;
                obj.measurements(:,idx) = meas;
                obj.control_vectors(:,idx) = u;
                obj.true_control(:,idx) = disturbance_fn(u);
            end
            %==============================

            %return obj
        end
    end
end