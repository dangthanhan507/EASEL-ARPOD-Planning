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
        function void = plot3d_to_2d(t,x,y,z)
            subplot(3,1,1)
            plot(t,x,'r');
            subplot(3,1,2)
            plot(t,y,'g');
            subplot(3,1,3)
            plot(t,z,'b');
        end
        function tstamps = steps2time(nsteps,tstep)
            tstamps = 0:tstep:(nsteps*tstep);
        end
        function void = plotSpacecraft(mhecell, mpccell, state)
            %{  
                mhecell: 6xN representing current mhe estimates
                mpccell: 6xN representing current mpc future expected states
                state    : 3x1 representing current spacecraft state
            %}
            [state_dim, nothing] = state;
            mhe_pos = mhecell(1:state_dim/2,:);
            mpc_pos = mpccell(1:state_dim/2,:);

            plotPos(mhe_pos, 'b--o');
            plotPos(mpc_pos, 'r--o');
            plotPos(state(1:state_dim/2,:), 'go');
        end
        function void = plotVecError(ts, gt, est, color)
            error = vecnorm(gt - est);
            plot(ts,error,color);
        end
        function fuel = getFuel(control,mass,tstep)
            fuel = mass*control*tstep;
        end
    end
    methods (Static)
        function obj = init(obj, benchmark);
            obj.benchmark = benchmark;
        end
        function void = graphErrors(obj)
            [control_dim, nsteps] = obj.benchmark.control_vectors;
            ts = Graph_Util.steps2time(nsteps, obj.benchmark.tstep);
            figure;
            hold on
            subplot(2,1,1);
            Graph_Util.plotVecError(ts, obj.benchmark.true_trajs, obj.benchmark.estimated_trajs);
            subplot(2,1,2);
            Graph_Util.plotVecError(ts, obj.benchmark.true_control, obj.benchmark.control_vectors);
            % Graph_Util.plotVecError(ts, obj.benchmark.true_att_trajs, obj.benchmark.est_att_trajs);
            hold off
        end
        function void = graphSpacecraftPos(obj, tidx)
            mheXs = obj.benchmark.mheXs{tidx};
            mpcXs = obj.benchmark.mpcXs{tidx};

            true_trajs = obj.true_trajs;
            est_trajs = obj.estimated_trajs;

            [true_pos, true_vel] = Graph_Util.decomposeTraj(true_trajs);
            [est_pos, est_vel]   = Graph_Util.decomposeTraj(est_trajs);

            figure;
            hold on
            %plot full traj
            Graph_Util.scatterPos(true_pos,'r');
            Graph_Util.scatterPos(est_pos,'b');

            %plot spacecraft currently
            Graph_Util.plotSpacecraft(mheXs, mpcXs, est_trajs(tidx,:));
            hold off
        end
        function void = graphTrajs(obj)
            true_trajs = obj.true_trajs;
            est_trajs = obj.estimated_trajs;

            [true_pos, true_vel] = Graph_Util.decomposeTraj(true_trajs);
            [est_pos, est_vel]   = Graph_Util.decomposeTraj(est_trajs);

            figure;
            hold on
            %plot full traj
            Graph_Util.scatterPos(true_pos,'r');
            Graph_Util.scatterPos(est_pos,'b');
            hold off
        end
    end
end