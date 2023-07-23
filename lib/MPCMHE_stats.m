classdef MPCMHE_stats
    properties
        predictedTraj

        trackTraj
        trackEstTraj

        fuelconsumed
        trackFuelConsumption

        trackControlDisturbance
        trackTrueControlDisturbance
        timestamps

        trackPhase
        gifctr % gif ctr
    end
    methods
        function obj = init(obj, window_Realstates, window_MHEstates, windowMPCstates, window_control, tstep, mass)
            obj.gifctr = 0;

            obj.predictedTraj = windowMPCstates;
            obj.trackTraj = window_Realstates;
            obj.trackEstTraj = window_MHEstates;

            obj.timestamp = 0;

            obj.fuelconsumed = 0;
            obj.trackFuelConsumption = 0;
            [dim,n] = size(window_control);
            for i = 1:n
                obj.fuelconsumed = obj.fuelconsumed + mass*window_control(:,i)*tstep;
                obj.trackFuelConsumption = [obj.trackFuelConsumption, obj.fuelconsumed];

                obj.timestamp = [obj.timestamp, obj.timestamp(:,end)+tstep];
            end

            obj.trackPhase = zeros(1,length(window_MHEstates));
            for i = 1:n
                obj.trackPhase(:,i) = ARPOD_Mission.calculatePhase(window_MHEstates(:,i));
            end
        end
        function obj = insertInfo(obj, control, trueTraj, window_estTraj, predictedTraj, window_disturbance, tstep, phase, mass)
            obj.trackPhase = [obj.trackPhase, phase];

            [dim,n] = size(window_disturbance);
            obj.trackControlDisturbance(:,end-(n-1)-1:end) = window_disturbance(:,1:end-1);
            obj.trackControlDisturbance = [obj.trackControlDisturbance, window_disturbance(:,end)];

            [dim,n] = size(window_estTraj);
            obj.trackEstTraj(:,end-(n-1)-1:end) = window_estTraj(:,1:end-1);
            obj.trackEstTraj = [obj.trackEstTraj, window_estTraj(:,end)];

            obj.fuelconsumed = obj.fuelconsumed + mass*control*tstep;
            obj.trackFuelConsumption = [obj.trackFuelConsumption, obj.fuelconsumed];

            obj.trackTraj = [obj.trackTraj, trueTraj];

            obj.predictedTraj = predictedTraj;
            obj.timestamp = [obj.timestamp, obj.timestamp(:,end)+tstep];
        end
        function obj = graph(obj)
            theta1 = ARPOD_Mission.theta*pi/180;
            theta2 = theta1;

            figure(1), clf

            filename = 'vipercylon2.gif'; %continuing tradition
            delaytime = 0.25;

            c = ARPOD_Benchmark.rho_d;
            %drawing pillars of pyramid
            plot3([0,-sin(theta2/2)*c], [0,-sin(theta1/2)*c], [0,-cos(theta1/2)*c], 'g');
            axis equal
            hold on

            %continuing to draw pillars
            plot3([0,-sin(theta2/2)*c], [0,sin(theta1/2)*c], [0,-cos(theta1/2)*c], 'g');
            plot3([0,sin(theta2/2)*c], [0,-sin(theta1/2)*c], [0,-cos(theta1/2)*c], 'g');
            plot3([0,sin(theta2/2)*c], [0,sin(theta1/2)*c], [0,-cos(theta1/2)*c], 'g');
            %drawing base of pyramid
            plot3([-sin(theta2/2)*c, sin(theta2/2)*c], [-sin(theta1/2)*c,-sin(theta1/2)*c], [-cos(theta1/2)*c,-cos(theta1/2)*c], 'g');
            plot3([-sin(theta2/2)*c, sin(theta2/2)*c], [sin(theta1/2)*c, sin(theta1/2)*c], [-cos(theta1/2)*c,-cos(theta1/2)*c], 'g');
            plot3([-sin(theta2/2)*c, -sin(theta2/2)*c], [-sin(theta1/2)*c,sin(theta1/2)*c], [-cos(theta1/2)*c,-cos(theta1/2)*c], 'g');
            plot3([sin(theta2/2)*c, sin(theta2/2)*c], [-sin(theta1/2)*c, sin(theta1/2)*c], [-cos(theta1/2)*c,-cos(theta1/2)*c], 'g');

            %draw chaser
            plot3(obj.trackEstTraj(1,:), obj.trackEstTraj(2,:), obj.trackEstTraj(3,:), 'b');
            plot3(obj.predictedTraj(1,:), obj.predictedTraj(2,:), obj.predictedTraj(3,:), 'b--');
            plot3(obj.trackTraj(1,:), obj.trackTraj(2,:), obj.trackTraj(3,:), 'go');
            %draw target
            plot3(0,0,0,'rd'); %draw target as diamond

            %incomplete
            %stil needs to show phase, fuel, 2d position graphs,
            %ControlDisturbances

            

            drawnow
            [A,map] = rgb2ind(frame2im(getframe(1)),256);
            if obj.gifctr == 0
                imwrite(A,map,filename,'gif','LoopCount',1,'DelayTime',delaytime);
            else
                obj.gifctr = obj.gifctr + 1;
                imwrite(A,map,filename,'fig','WriteMode','append','Delaytime',delaytime);
            end
        end
    end
end