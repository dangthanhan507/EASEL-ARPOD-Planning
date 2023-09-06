classdef Graph_Util
    properties 
        benchmark
    end
    %helper methods are static
    methods (Static)
        function [pos,vel] = decomposeTraj(trajs)
            [state_size, n_trajs] = size(trajs);
            pos = trajs(1:state_size/2,1:end);
            vel = trajs(state_size/2:end, 1:end);
        end
        function void = scatterPos(pos,color)
            x = pos(1,:);
            y = pos(2,:);
            z = pos(3,:);
            scatter3(x,y,z,color,'filled');
        end
        function void = plotPos(pos, color)
            x = pos(1,:);
            y = pos(2,:);
            z = pos(3,:);
            plot3(x,y,z,color);
        end
        function void = plot3d_to_2d(t,x,y,z, colors)
            disp(colors(1))
            subplot(3,1,1)
            plot(t,x,colors(1));
            subplot(3,1,2)
            plot(t,y,colors(2));
            subplot(3,1,3)
            plot(t,z,colors(3));
        end
        function tstamps = steps2time(nsteps,tstep)
            tstamps = 0:tstep:(nsteps*tstep);
            tstamps = tstamps(:,1:end-1);
        end
        function void = plotSpacecraft(mhecell, mpccell, state)
            %{  
                mhecell: 6xN representing current mhe estimates
                mpccell: 6xN representing current mpc future expected states
                state    : 3x1 representing current spacecraft state
            %}
            [state_dim, nothing] = size(state);
            mhe_pos = mhecell(1:3,:);
            mpc_pos = mpccell(1:3,:);

            Graph_Util.plotPos(mhe_pos, 'b--o');
            Graph_Util.plotPos(mpc_pos, 'r--o');
            Graph_Util.plotPos(state(1:3,:), 'go');
        end
        function void = plotVecError(ts, gt, est, color)
            error = vecnorm(gt - est);
            plot(ts,error,color);
        end
        function void = plotFuel(controls,mass,tstep)
            [control_dim, horizon] = size(controls);
            fuel = mass*tstep*abs(controls);
            ts   = linspace(0,horizon*tstep,horizon);

            fuel = sum(fuel,1);
            plot(ts,fuel)
        end
    end
    methods
        function obj = init(obj, benchmark);
            obj.benchmark = benchmark;
        end
        function void = graphErrors(obj)
            [control_dim, nsteps] = size(obj.benchmark.control_vectors);
            ts = Graph_Util.steps2time(nsteps, obj.benchmark.tstep);
            figure;
            subplot(3,1,1);
            Graph_Util.plotVecError(ts, obj.benchmark.true_trajs, obj.benchmark.estimated_trajs, 'r-');
            subplot(3,1,2);
            Graph_Util.plotVecError(ts, obj.benchmark.true_control, obj.benchmark.control_vectors, 'b-');
            subplot(3,1,3);
            Graph_Util.plotVecError(ts, obj.benchmark.true_att_trajs, obj.benchmark.est_att_trajs, 'g-');
        end
        function void = graphSpacecraftPos(obj, tidx)
            mheXs = obj.benchmark.mheXs{tidx};
            mpcXs = obj.benchmark.mpcXs{tidx};

            true_trajs = obj.benchmark.true_trajs;
            est_trajs = obj.benchmark.estimated_trajs;

            [true_pos, true_vel] = Graph_Util.decomposeTraj(true_trajs);
            [est_pos, est_vel]   = Graph_Util.decomposeTraj(est_trajs);

            figure;
            %plot full traj
            Graph_Util.plotPos(true_pos,'r');
            hold on
            Graph_Util.plotPos(est_pos,'b');

            %plot spacecraft currently
            Graph_Util.plotSpacecraft(mheXs, mpcXs, est_trajs(:,tidx));
            hold off
        end
        function void = graphFuel(obj)
            controls = obj.benchmark.control_vectors;
            mass = ARPOD_Mission.m_c;
            tstep = obj.benchmark.tstep;
            figure
            Graph_Util.plotFuel(controls,mass,tstep);
        end
        function void = graphTrajs2d(obj)
            true_trajs = obj.benchmark.true_trajs;
            est_trajs  = obj.benchmark.estimated_trajs;
            [true_pos, true_vel] = Graph_Util.decomposeTraj(true_trajs);
            [est_pos, est_vel]   = Graph_Util.decomposeTraj(est_trajs);
            [control_dim, nsteps] = size(obj.benchmark.control_vectors);
            t = Graph_Util.steps2time(nsteps, obj.benchmark.tstep);

            figure;
            %plot full traj
            subplot(3,1,1)
            plot(t,true_pos(1,:),"r-");
            hold on
            plot(t,est_pos(1,:),"ro");
            hold off
            subplot(3,1,2)
            plot(t,true_pos(2,:),"g-");
            hold on
            plot(t,est_pos(2,:),"go");
            hold off
            subplot(3,1,3)
            plot(t,true_pos(3,:),"b-");
            hold on 
            plot(t,est_pos(3,:),"bo");
            hold off
        end
        function void = graphTrajs(obj)
            true_trajs = obj.benchmark.true_trajs;
            est_trajs = obj.benchmark.estimated_trajs;

            [true_pos, true_vel] = Graph_Util.decomposeTraj(true_trajs);
            [est_pos, est_vel]   = Graph_Util.decomposeTraj(est_trajs);

            figure;
            %plot full traj
            Graph_Util.plotPos(true_pos,'r');
            hold on
            Graph_Util.plotPos(est_pos,'b');
            hold off
        end
    end
end