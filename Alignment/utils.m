classdef utils
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
        %%%%%%%%%% DYNAMICS FOR MPCMHE %%%%%%%%%%
        function [A,B] = linearHCWDynamics(T)
            %{
                Description:
                ------------
                Create relative orbital dynamics discrete propagation matrices
                using HCW equations. (A,B)

                Parameters:
                -----------
                @params T    : (scalar) time step that A,B propagate by
                @params mu_GM: (scalar) gravitational constant used in ARPOD
                @params R    : (scalar) radial distance from spacecraft to earth
                @params is2D : (bool) whether to use 2d version or not.    

                Returns:
                --------
                returns (A,B) matrices to use in xk+1 = A*xk + B*uk
            %}
            mu_GM = utils.mu; % gravitational constant we define in mission.
            R = utils.a; % orbital radius or semi-major axis

            n = sqrt(mu_GM / (R.^3) );
            A = zeros(6,6);
            B = zeros(6,3);
            S = sin(n * T);
            C = cos(n * T);

            %ZOH on HCW ODE produces these matrices (A,B)
            A(1,:) = [4-3*C,0,0,S/n,2*(1-C)/n,0];
            A(2,:) = [6*(S-n*T),1,0,-2*(1-C)/n,(4*S-3*n*T)/n,0];
            A(3,:) = [0,0,C,0,0,S/n];
            A(4,:) = [3*n*S,0,0,C,2*S,0];
            A(5,:) = [-6*n*(1-C),0,0,-2*S,4*C-3,0];
            A(6,:) = [0,0,-n*S,0,0,C];

            B(1,:) = [(1-C)/(n*n),(2*n*T-2*S)/(n*n),0];
            B(2,:) = [-(2*n*T-2*S)/(n*n),(4*(1-C)/(n*n))-(3*T*T/2),0];
            B(3,:) = [0,0,(1-C)/(n*n)];
            B(4,:) = [S/n,2*(1-C)/n, 0];
            B(5,:) = [-2*(1-C)/n,(4*S/n) - (3*T),0];
            B(6,:) = [0,0,S/n];

        end

        function [Ak,Bk] = attitudeDiscrete(T)
            %{
                Description:
                ------------
                Create simple attitude dynamics for LVLH system. This represents relative
                attitude dynamics to the target/chief spacecraft.

                Parameters:
                -----------
                @params T    : (scalar) time step that A,B propagate by
                @params is2D : (bool) whether to use 2d version or not.    
                
                Returns:
                --------
                returns (A,B) matrices to use in xk+1 = A*xk + B*uk
            %}
            %define the dynamics xdot = Ax + Bu
            %roll,pitch,yaw,rolldot,pitchdot,yawdot

            % x (6,1)
            % u (3,1)
            A = [zeros(3,3),eye(3);zeros(3,6)];
            B = [zeros(3,3);eye(3)];
            statespace = ss(A,B,eye(6),zeros(6,3));
            out = c2d(statespace,T);

            Ak = out.A;
            Bk = out.B;
        end

        %%%%%%%%%% SIMULATION CODE %%%%%%%%%%
            %%%%%%%%% ODE SOLVE %%%%%%%%%
        %{
            Nonlinear Relative Motion ODE Solve for Simulation
        %}
        function dtraj = nonlinearMotion(t,traj0, u)
            x = traj0(1);
            y = traj0(2);
            z = traj0(3);
            
            
            ut = u; %function
            ux = ut(1); %indexing thrusters
            uy = ut(2);
            uz = ut(3);

            mu_GM = utils.mu; % gravitational constant we define in mission.
            R = utils.a; % orbital radius or semi-major axis

            n = sqrt(mu_GM / (R.^3)); %orbital velocity
        
            %distance formula on chaser orbital radius ^3
            %resembles gravitational formula but generalized for 3d.
            const = ((R+x).^2 + y.^2 + z.^2).^(3/2); 
            
            
            xdot = traj0(4);
            ydot = traj0(5);
            zdot = traj0(6);

            xdotdot = 2*n*ydot + n*n*(R+x) - mu_GM*(R+x)/const + ux;
            ydotdot = -2*n*xdot + n*n*y - mu_GM*y/const + uy;
            zdotdot = -mu_GM*z/const + uz;
        
            dtraj = [xdot;ydot;zdot;xdotdot;ydotdot;zdotdot];
        end
        function traj = nonlinearMotionSolver(traj0, u, tstep)
            motion = @(t,traj) utils.nonlinearMotion(t, traj, u);
            [ts, trajs] = ode45(motion, 0:tstep, traj0);
            traj = transpose(trajs(length(trajs),:));
        end

        %{
            Attitude Dynamics Solve for Simulation
        %}
        function dtraj = attitudeMotion(t, traj0, u)
            %{
                no natural dynamics occurring.
                acceleration is mainly due to actuation.
            %}
            ut = u;
            roll = traj0(1);
            pitch = traj0(2);
            yaw = traj0(3);

            rolldot = traj0(4);
            pitchdot = traj0(5);
            yawdot = traj0(6);

            u_roll = ut(1);
            u_pitch = ut(2);
            u_yaw = ut(3);

            rolldotdot = u_roll;
            pitchdotdot = u_pitch;
            yawdotdot = u_yaw;

            dtraj = [rolldot;pitchdot;yawdot;rolldotdot;pitchdotdot;yawdotdot];
        end
        function traj = attitudeSolver(traj0, u, tstep)
            motion = @(t,traj) utils.attitudeMotion(t,traj,u);
            [ts,trajs] = ode45(motion, 0:tstep, traj0);
            traj = transpose(trajs(length(trajs),:));
        end

            %%%%%%%%% SENSING %%%%%%%%%
        function z_t = measure(state)
            %{
                Assume we are in Phase 2
            %}
            x = state(1,:);
            y = state(2,:);
            z = state(3,:);


            norm = sqrt(x*x+y*y+z*z);
            e1 = atan(y/x);
            e2 = asin(z/norm);
            e3 = norm;
            z_t = [e1;e2;e3];
        end
        function jacobian = measureJacobian(state)
            x = state(1,:);
            y = state(2,:);
            z = state(3,:);
            
            jacobian = zeros(3,6);
            %dArctan
            partialX = -y/(x*x+y*y);
            partialY = x/(x*x+y*y);
            jacobian(1,1) = partialX;
            jacobian(1,2) = partialY;
            %rest of the partials are zero so it doesn't matter anyways.

            %dArcsin
            norm = sqrt( (x*x+y*y)/(x*x+y*y+z*z) ) * (x*x+y*y+z*z).^(3/2);
            partialX = -x*z/norm;
            partialY = -y*z/norm;
            partialZ = sqrt( (x*x+y*y)/(x*x+y*y+z*z) ) / sqrt(x*x+y*y+z*z);
            jacobian(2,1) = partialX;
            jacobian(2,2) = partialY;
            jacobian(2,3) = partialZ;

            %drho
            norm = sqrt(x*x+y*y+z*z);
            partialX = x/norm;
            partialY = y/norm;
            partialZ = z/norm;
            jacobian(3,1) = partialX;
            jacobian(3,2) = partialY;
            jacobian(3,3) = partialZ;

            %return jacobian
        end

        function rot = rotRPY(theta)
            R11 = cos(theta(3,:)).*cos(theta(2,:));
            R12 = -sin(theta(3,:)).*cos(theta(1,:)) + cos(theta(3,:)).*sin(theta(2,:)).*sin(theta(1,:));
            R13 = sin(theta(3,:)).*sin(theta(1,:)) + cos(theta(3,:)).*sin(theta(2,:)).*cos(theta(1,:));

            R21 = sin(theta(3,:)).*cos(theta(2,:));
            R22 = cos(theta(3,:)).*cos(theta(1,:)) + sin(theta(3,:)).*sin(theta(2,:)).*sin(theta(1,:));
            R23 = -cos(theta(3,:)).*sin(theta(1,:)) + sin(theta(3,:)).*sin(theta(2,:)).*cos(theta(1,:));

            R31 = -sin(theta(2,:));
            R32 = cos(theta(2,:)).*sin(theta(1,:));
            R33 = cos(theta(2,:)).*cos(theta(1,:));

            rot = [R11, R12, R13;...
                   R21, R22, R23;...
                   R31, R32, R33];
        end
    end
end