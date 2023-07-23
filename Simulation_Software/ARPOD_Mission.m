classdef ARPOD_Mission
    %TODO: clean this up -> AN
    properties (Constant)
        t_e = 14400; % eclipse time (in seconds)
        t_f = 43200; % total mission duration (in seconds)
        rho_r = 1; % maximum distance for range measurements (1 km)
        rho_d = 0.1; % docking phase initial radius (0.1 km)
        m_t = 2000; % mass of target (2000 kg)
        m_c = 500; % mass of chaser (500 kg)
        mu = 398600.4; %earth's gravitational constant (in km^3/s^2)
        a = 42164; % semi-major axis of GEO (42164 km)
        Vbar = 5 * 10.^(-5); % max closing velocity while docking (in km/s)
        theta = 60; % LOS Constraint angle (in degrees)
        c = [-1;0;0]; % LOS cone direction
        x_docked = [0;0;0;0;0;0]; % docked position in km, and km/s
        x_relocation = [0;20;0;0;0;0]; %relocation position in km, km/s
        x_partner = [0;30;0;0;0;0]; %partner position in km, km/s
        ubar = 0.1; % limit of the thrust allowed
    end
    methods (Static)
        function inLOS = isInsideLOS(traj)
            theta1 = ARPOD_Mission.theta * pi / 180; %docking angle in radians
            theta2 =  theta1;
            LOS_mtx = [ sin(theta1/2), cos(theta1/2), 0; 
                    sin(theta1/2), -cos(theta1/2), 0;
                    sin(theta2/2), 0, cos(theta2/2);
                    sin(theta2/2), 0, -cos(theta2/2)];
            xyz = [traj(1);traj(2);traj(3)];
            b = LOS_mtx*xyz;
            if b(1) <= 0 && b(2) <= 0 && b(3) <= 0 && b(4) <= 0
                inLOS = 1;
            else
                inLOS = 0;
            end
        end
        function inLOS = isInsideConeLOS(traj)
            theta = ARPOD_Mission.theta * pi / 180; %docking angle in radians
            rho_d = ARPOD_Mission.rho_d;
            c = rho_d * ARPOD_Mission.c;
            rho = [traj(1);traj(2);traj(3)];
            
            dot = (rho.' * c) / (norm(rho)*norm(c));
            if (dot(1) >= cos(theta/2))
                % within LOS
                inLOS = 1;
            else
                % not within LOS
                inLOS = 0;
            end
        end
        
        function phase = calculatePhase(traj)
            norm = traj(1:3,:);
            norm = sqrt(sum(norm.^2));
            if (norm > ARPOD_Mission.rho_r)
                % ARPOD phase 1: Rendezvous w/out range
                phase = 1;
            elseif ( (norm > ARPOD_Mission.rho_d) || (ARPOD_Mission.isInsideLOS(traj) == 0) ) 
                % ARPOD phase 2: Rendezvous with range
                phase = 2;
            else 
                %ARPOD phase 3: Docking
                phase = 3;
            end
        end

    end
    properties
        dynamics % using dynamics from ARPOD_Dynamics
        traj % current trajectory of mission
        sensor % current measurements from sensors
        phase % current phase of mission
        inLOS % boolean as whether it is in or not
        tstep
        los_type %type of los, linear or nonlinear
    end
    methods
        function obj = initMission(obj, traj, tstep, is2D, thruster_type)
            obj.traj = traj;
            obj.phase = ARPOD_Mission.calculatePhase(traj);
            obj.inLOS = ARPOD_Mission.isInsideLOS(traj);
            obj.tstep = tstep;

            obj.dynamics = ARPOD_Dynamics;
            obj.dynamics = obj.dynamics.initDynamics(traj, is2D, thruster_type);
        end
        function obj = nextStep(obj, control, system_noise, sensor_noise, useFullSense)
            obj.dynamics = obj.dynamics.nextStep(control, obj.tstep, ARPOD_Mission.a, ARPOD_Mission.mu, system_noise);
            obj.traj = obj.dynamics.currentTraj();

            obj.phase = ARPOD_Mission.calculatePhase(obj.traj);

            if useFullSense
                obj.sensor = ARPOD_Sensor.fullSense(obj.traj, sensor_noise);
            else
                obj.sensor = ARPOD_Sensor.sense(obj.traj, sensor_noise, obj.phase);
            end
            obj.inLOS = ARPOD_Mission.isInsideLOS(obj.traj);
        end
    end
end