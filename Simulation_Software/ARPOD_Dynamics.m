% ============================================================
% ARPOD_Dynamics.m
% 
% This file specifically should contain all things dynamics
% in the relative orbital motion problem. 
%
% TODO: support 2d dynamics, attitude dynamics and choose
%       between all of them.
% 
% Author: An Dang
%
% ============================================================

classdef ARPOD_Dynamics
    % gravitational constant in km^2/s^2 from chance-constr MPC
    properties 
        traj
        att_traj
        is2D %determine whether simulation is 2d or 3d
        thruster_type %determine thruster type in simulation
        use_attitude = false;

    end
    methods (Static)
        function [A,B] = linearHCWDynamics(T, mu_GM, R, is2D)
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

            if is2D
                A = [A(1:2,1:2), A(1:2,4:5); A(4:5,1:2), A(4:5,4:5)];
                B = [B(1:2,1:2); B(4:5,1:2)];
            end
        end
        function [Ak,Bk] = attitudeLVLH(T, is2D)
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
            if is2D
                %yaw, yawdot

                % x (2,1)
                % u (1,1)
                A = [0,1;0,0];
                B = [0;1];
                statespace = ss(A,B,eye(2),zeros(2,1));
                out = c2d(statespace,T);
            else
                %roll,pitch,yaw,rolldot,pitchdot,yawdot

                % x (6,1)
                % u (3,1)
                A = [zeros(3,3),eye(3);zeros(3,6)];
                B = [zeros(3,3);eye(3)];
                statespace = ss(A,B,eye(6),zeros(6,3));
                out = c2d(statespace,T);
            end
            Ak = out.A;
            Bk = out.B;
        end
        function dtraj = attitudeMotion(t, traj0, u, is2D)
            ut = u(t);
            if is2D
                yaw = traj0(1);

                yawdot = traj0(2);

                u_yaw = ut(1);

                yawdotdot = u_yaw;

                dtraj = [yawdot;yawdotdot];
            else
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
        end
        function dtraj = nonlinearMotion(t,traj0, mu_GM, R, u, is2D)
            x = traj0(1);
            y = traj0(2);
            if ~is2D
                z = traj0(3);
            end
            
            
            ut = u(t); %function
            ux = ut(1); %indexing thrusters
            uy = ut(2);
            if ~is2D
                uz = ut(3);
            end

            n = sqrt(mu_GM / (R.^3)); %orbital velocity
        
            %distance formula on chaser orbital radius ^3
            %resembles gravitational formula but generalized for 3d.
            if ~is2D
                const = ((R+x).^2 + y.^2 + z.^2).^(3/2); 
            else
                const = ((R+x).^2 + y.^2).^(3/2); 
            end
            
            
            if ~is2D
                xdot = traj0(4);
                ydot = traj0(5);
                zdot = traj0(6);
            else
                xdot = traj0(3);
                ydot = traj0(4);
            end

            xdotdot = 2*n*ydot + n*n*(R+x) - mu_GM*(R+x)/const + ux;
            ydotdot = -2*n*xdot + n*n*y - mu_GM*y/const + uy;
            if ~is2D
                zdotdot = -mu_GM*z/const + uz;
            end
            

            if is2D
                dtraj = [xdot;ydot;xdotdot;ydotdot];
            else
                dtraj = [xdot;ydot;zdot;xdotdot;ydotdot;zdotdot];
            end
        end
        function traj = nonlinearMotionSolver(traj0, R, mu_GM, u, tstep, is2D)
            motion = @(t,traj) ARPOD_Dynamics.nonlinearMotion(t, traj, mu_GM, R, u, is2D);
            [ts, trajs] = ode45(motion, 0:tstep, traj0);
            traj = transpose(trajs(length(trajs),:));
        end
        function traj = attitudeSolver(traj0, u, tstep, is2D)
            motion = @(t,traj) ARPOD_Dynamics.attitudeMotion(t,traj,u,is2D);
            [ts,trajs] = ode45(motion, 0:tstep, traj0);
            traj = transpose(trajs(length(trajs),:));
        end
        function jacobianMat = motionJacobian(t,traj0, mu_GM, R,u)
            n = sqrt(mu_GM / (R.^3)); %orbital velocity
            x = traj0(1);
            y = traj0(2);
            z = traj0(3);

            jacobianMat = zeros(6,6);

            jacobianMat(1,4) = 1;
            jacobianMat(2,5) = 1;
            jacobianMat(3,6) = 1;

            norm = ( (R+x).^2 + y*y + z*z ).^(5/2);
            jacobianMat(4,1) = n*n - mu_GM*((-2*R*R)-4*R*x-2*x*x+y*y+z*z)/norm;
            jacobianMat(4,2) = mu_GM*(3*y*(R+x))/norm;
            jacobianMat(4,3) = mu_GM*(3*z*(R+x))/norm;
            jacobianMat(4,5) = 2*n;

            jacobianMat(5,1) = mu_GM*(3*y*(R+x))/norm;
            jacobianMat(5,2) = n*n - mu_GM*(R*R + 2*R*x + x*x - 2*y*y + z*z) / norm;
            jacobianMat(5,3) = mu_GM * 3*y*z/norm;
            jacobianMat(5,4) = -2*n;

            jacobianMat(6,1) = mu_GM*3*z*(R+x)/norm;
            jacobianMat(6,2) = mu_GM*3*y*z/norm;
            jacobianMat(6,3) = -mu_GM*(R*R+2*R*x+x*x+y*y-2*z*z)/norm;
        end
    end
    methods 
        function obj = initDynamics(obj, traj0, is2D, thruster_type)
            obj.is2D = is2D;
            obj.traj = traj0;
            obj.thruster_type = thruster_type;
        end
        function obj = initAttitude(obj, traj0)
            obj.att_traj = traj0;
            obj.use_attitude = true;
        end
        function obj = nextStep(obj, control_input, tstep, R, mu_GM, noise_model)
            %{
                Thruster Types:
                    Discrete-Time
                    Impulsive
            %}
            if obj.thruster_type == 1
                %discrete time
                ut = @(t) control_input;
            else
                %impulsive
                ut = @(t) [0;0;0];
                obj.traj(4:6) = obj.traj(4:6) + control_input;
            end
            obj.traj = ARPOD_Dynamics.nonlinearMotionSolver(obj.traj, R, mu_GM, ut, tstep, obj.is2D) + noise_model();
        end
        function obj = nextAttitudeStep(obj, control_input, tstep, noise_model)
            ut = @(t) control_input;
            obj.att_traj = ARPOD_Dynamics.attitudeSolver(obj.att_traj, ut, tstep, obj.is2D) + noise_model();
        end
        function traj = currentTraj(obj)
            traj = obj.traj;
        end
        function traj = currentAtt(obj)
            traj = obj.att_traj;
        end
    end
end