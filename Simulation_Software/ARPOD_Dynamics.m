% gravitational constant in km^2/s^2 from chance-constr MPC

classdef ARPOD_Dynamics
    properties 
        traj
    end
    methods (Static)
        function [A,B] = linearHCWDynamics(T)
            mu_GM = ARPOD_Mission.mu;
            R = ARPOD_Mission.a;

            n = sqrt(mu_GM / (R.^3) );
            A = zeros(6,6);
            B = zeros(6,3);
            S = sin(n * T);
            C = cos(n * T);

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
        function dtraj = nonlinearMotion(t,traj0, mu_GM, R, u)
            x = traj0(1);
            y = traj0(2);
            z = traj0(3);
            
            ut = u(t); %function
            ux = ut(1); %indexing thrusters
            uy = ut(2);
            uz = ut(3);

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
        function traj = nonlinearMotionSolver(traj0, R, mu_GM, u, tstep)
            motion = @(t,traj) ARPOD_Dynamics.nonlinearMotion(t, traj, mu_GM, R, u);
            [ts, trajs] = ode45(motion, 0:tstep, traj0);
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
        function obj = initDynamics(obj, traj0)
            obj.traj = traj0;
        end
        function obj = nextStep(obj, control_input, tstep, R, mu_GM, noise_model, thruster_type)
            %{
                Thruster Types:
                    Discrete-Time
                    Impulsive
            %}
            if thruster_type == 1
                ut = @(t) control_input;
            else
                ut = @(t) [0;0;0];
                obj.traj(4:6) = obj.traj(4:6) + control_input;
            end
            obj.traj = ARPOD_Dynamics.nonlinearMotionSolver(obj.traj, R, mu_GM, ut, tstep) + noise_model();
        end
        function traj = currentTraj(obj)
            traj = obj.traj;
        end
    end
end