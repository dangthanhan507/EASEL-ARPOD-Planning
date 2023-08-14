classdef GraphUtil
    methods (Static)

        %%%%% HELPER CODE %%%%%
        function void = setLatex()
            set(groot,'defaultAxesTickLabelInterpreter','latex');  
            set(groot,'defaulttextinterpreter','latex');
            set(groot,'defaultLegendInterpreter','latex');
        end
        function [r,g,b] = randomize_color()
            [r,g,b] = [rand,rand,rand];
        end
        function void = graphPtDistribution(x,y,z,color, first_graph)
            %{
                x: list of x trajectories (1,length)
                y: list of y trajectories (1,length)
                z: list of z trajectories (1,length)
                color: char -> 'r' is red 'b' is blue 'bl' is black
            %}
            scatter3(x,y,z, color, 'filled');
            if first_graph
                %axis equal
                hold on
            end
        end
        function void = graph2dPt(x,y,color,first_graph)
            scatter(x,y,color);
            if first_graph
                hold on
            end
        end
        function void = graphTrajectory(x,y,z,color,first_graph)
            plot3(x,y,z,color);
            if first_graph
                %axis equal
                hold on
            end
        end
        function void = graphTranslation(x,y,z,first_graph)
            t = length(x);
            ts = linspace(0,t,t);
            plot(ts,x,'r');
            if first_graph
                %axis equal
                hold on
            end
            plot(ts,y,'g');
            plot(ts,z,'b');

        end
        %%%%% END OF HELPER CODE %%%%%

        %%%%% START OF UTIL CODE %%%%%
        function void = graphInitialDistribution(arpod_stats)
            %{
                arpod_stats: Arpod_Statistics Data Structure array
            %}
            xs = zeros(length(arpod_stats));
            ys = zeros(length(arpod_stats));
            zs = zeros(length(arpod_stats));
            for idx = 1:length(arpod_stats)
                tracked_trajs = arpod_stats{idx}.trackTraj; %retrieve trajectories
                xs(idx) = tracked_trajs(1,1);
                ys(idx) = tracked_trajs(2,1);
                zs(idx) = tracked_trajs(3,1);
            end
            axis equal
            GraphUtil.graphPtDistribution(xs,ys,zs,'r', true);
            xlabel("X [km]");
            ylabel("Y [km]");
            zlabel("Z [km]");
            hold off
        end

        function void = graphFinalDistribution(arpod_stats)
            %{
                arpod_stats: Arpod_Statistics Data Structure array
            %}
            xs = zeros(length(arpod_stats));
            ys = zeros(length(arpod_stats));
            zs = zeros(length(arpod_stats));
            for idx = 1:length(arpod_stats)
                tracked_trajs = arpod_stats{idx}.trackTraj; %retrieve trajectories
                xs(idx) = tracked_trajs(1,end);
                ys(idx) = tracked_trajs(2,end);
                zs(idx) = tracked_trajs(3,end);
            end
            axis equal
            GraphUtil.graphPtDistribution(xs,ys,zs,'b', true);
            xlabel("X [km]");
            ylabel("Y [km]");
            zlabel("Z [km]");
            hold off
        end 

        function void = graphStartFinishDistribution(arpod_stats, first_graph)
            xs = zeros(length(arpod_stats));
            ys = zeros(length(arpod_stats));
            zs = zeros(length(arpod_stats));
            for idx = 1:length(arpod_stats)
                tracked_trajs = arpod_stats{idx}.trackTraj; %retrieve trajectories
                xs(idx) = tracked_trajs(1,end);
                ys(idx) = tracked_trajs(2,end);
                zs(idx) = tracked_trajs(3,end);
            end

            scatter3(xs,ys,zs, 'b', 'filled');
            if first_graph
                hold on
            end

            xs = zeros(length(arpod_stats));
            ys = zeros(length(arpod_stats));
            zs = zeros(length(arpod_stats));
            for idx = 1:length(arpod_stats)
                tracked_trajs = arpod_stats{idx}.trackTraj; %retrieve trajectories
                xs(idx) = tracked_trajs(1,1);
                ys(idx) = tracked_trajs(2,1);
                zs(idx) = tracked_trajs(3,1);
            end

            if first_graph
                hold on
            end
            scatter3(xs,ys,zs, 'r', 'filled');

            xlabel("X [km]");
            ylabel("Y [km]");
            zlabel("Z [km]");
            hold off
        end

        function void = graphTranslationDist(arpod_stats, first_graph)  
            for idx = 1:length(arpod_stats)
                track_trajs = arpod_stats{idx}.trackTraj;
                xs = track_trajs(1,:);
                xs = xs(:);
                ys = track_trajs(2,:);
                ys = ys(:);
                zs = track_trajs(3,:);
                zs = zs(:);
                GraphUtil.graphTranslation(xs,ys,zs,first_graph);
                if first_graph
                    first_graph = false;
                end
            end
            hold off;
        end

        function void = graphTrajectoryDist(arpod_stats)
            first_graph = true;
            for idx = 1:length(arpod_stats)
                track_trajs = arpod_stats{idx}.trackTraj;
                xs = track_trajs(1,:);
                xs = xs(:);
                ys = track_trajs(2,:);
                ys = ys(:);
                zs = track_trajs(3,:);
                zs = zs(:);
                plot3(xs,ys,zs,'Color',[rand,rand,rand]);
                if first_graph
                    %axis equal
                    hold on 
                    first_graph = false;
                end
            end
            GraphUtil.graphTarget();
            hold off
        end

        function void = graph2DViewMission(arpod_stat, title_)
            trackTraj = arpod_stat.trackTraj;
            trackU = [arpod_stat.trackU, arpod_stat.trackU(:,end)];
            xs = trackTraj(1,:);
            ys = trackTraj(2,:);
            zs = trackTraj(3,:);
            
            vxs = trackTraj(4,:);
            vys = trackTraj(5,:);
            vzs = trackTraj(6,:);

            uxs = trackU(1,:);
            uys = trackU(2,:);
            uzs = trackU(3,:);

            ts = arpod_stat.timestamps;
            subplot(3,1,1)
            plot(ts, xs, 'r','DisplayName','x');
            hold on
            plot(ts,ys,'g','DisplayName','y');
            plot(ts,zs,'b','DisplayName','z');
            legend('x','y','z');c
            hold off

            subplot(3,1,2)
            plot(ts, vxs, 'r','DisplayName','vx');
            hold on
            plot(ts,vys,'g','DisplayName','vy');
            plot(ts,vzs,'b','DisplayName','vz');
            legend('vx','vy','vz');
            hold off

            subplot(3,1,3)
            plot(ts, uxs, 'r','DisplayName','ux');
            hold on
            plot(ts,uys,'g','DisplayName','uy');
            plot(ts,uzs,'b','DisplayName','uz');
            legend('ux','uy','uz');
            hold off
        end

        function void = graph2DViewMissions(arpod_stats,title_,folder,save,xlim_)
            
            figure;
            first_graph = true;
            for i = 1:length(arpod_stats)
                arpod_stat = arpod_stats{i};
                ts = arpod_stat.timestamps;
                trackTraj = arpod_stat.trackTraj;
                xs = trackTraj(1,:);
                ys = trackTraj(2,:);
                zs = trackTraj(3,:);

                plot(ts, xs, 'r');
                if first_graph
                    first_graph = false;
                    hold on
                end
                plot(ts,ys,'g');
                plot(ts,zs,'b');
            end
            legend('x','y','z');
            xlabel('Time [s]');
            ylabel('Pos [km]');
            xlim([0,xlim_]);
            title(title_+" Position");
            hold off

            if save
                saveas(gcf,folder+title_+'_2d_pos.eps','epsc');
            end

            figure;
            first_graph = true;
            for i = 1:length(arpod_stats)
                arpod_stat = arpod_stats{i};
                ts = arpod_stat.timestamps;
                trackTraj = arpod_stat.trackTraj;
                vxs = trackTraj(4,:);
                vys = trackTraj(5,:);
                vzs = trackTraj(6,:);

                plot(ts, vxs, 'r');
                if first_graph
                    first_graph = false;
                    hold on
                end
                plot(ts,vys,'g');
                plot(ts,vzs,'b');
            end
            legend('vx','vy','vz');
            xlabel('Time [s]');
            ylabel('Velocity [km/s]');
            title(title_+" Velocity");
            xlim([0,xlim_]);
            hold off
            if save
                saveas(gcf,folder+title_+'_2d_vel.eps','epsc');
            end

            figure;
            first_graph = true;
            for i = 1:length(arpod_stats)
                arpod_stat = arpod_stats{i};
                ts = arpod_stat.timestamps;
                trackU = [arpod_stat.trackU, arpod_stat.trackU(:,end)];
                uxs = trackU(1,:);
                uys = trackU(2,:);
                uzs = trackU(3,:);

                plot(ts, uxs, 'r');
                if first_graph
                    first_graph = false;
                    hold on
                end
                plot(ts,uys,'g');
                plot(ts,uzs,'b');
            end
            legend('ux','uy','uz');
            xlabel('Time [s]');
            ylabel('Thrust [km/s^2]');
            title(title_+" Thrust");
            xlim([0,xlim_]);
            hold off
            if save
                saveas(gcf,folder+title_+'_2d_thr.eps','epsc');
            end
        end

        function percentages = normalizeTrajectories(track_trajs)
            num_steps = length(track_trajs);
            distance = 0;
            percentages = zeros(1,num_steps);
            for idx = 2:num_steps
                %distance = distance + sqrt( sum( ((track_trajs(idx) - track_trajs(idx-1)).^2).' ) );
                distance = distance + norm(track_trajs(idx) - track_trajs(idx-1));
                percentages(idx) = distance;
            end
            percentages = percentages / distance;
        end
        function void = graphErrorSpecial(arpod_stats)
            first_graph = true;
            sum_trajs = 0;
            percent_trajs = [];
            errors_ys = [];
            for idx = 1:length(arpod_stats)
                track_trajs = arpod_stats{idx}.trackTraj;
                %normalize trajectories (get x axis stuff)
                percent_traj = GraphUtil.normalizeTrajectories(track_trajs);

                %get y axis stuff
                errors_y = sum( (arpod_stats{idx}.trackEstTraj - arpod_stats{idx}.trackTraj).^2 );

                percent_trajs = [percent_trajs, percent_traj];
                errors_ys = [errors_ys, errors_y];
            end
            GraphUtil.graph2dPt(percent_traj, errors_y,'r',first_graph);
            hold off
        end
        
        function void = graphDeltaFuelSpecial(arpod_stats)
            first_graph = true;
            sum_trajs = 0;
            percent_trajs = [];
            fuel_ys = [];
            for idx = 1:length(arpod_stats)
                track_trajs = arpod_stats{idx}.trackTraj;
                track_fuel = arpod_stats{idx}.trackFuelConsumption;
                %normalize trajectories (get x axis stuff)
                percent_traj = GraphUtil.normalizeTrajectories(track_trajs);
                kernel = [1,-1];
                %get y axis stuff
                track_fuel = conv(track_fuel,kernel,"same");
                track_fuel = track_fuel(:,1:end-1);

                percent_trajs = [percent_trajs, percent_traj];
                fuel_ys = [fuel_ys, track_fuel];
            end
            GraphUtil.graph2dPt(percent_traj, errors_y,'r',first_graph);
            hold off
        end

        function void = graphDeltaErrorSpecial(arpod_stats)
            first_graph = true;
            sum_trajs = 0;
            percent_trajs = [];
            errors_ys = [];
            for idx = 1:length(arpod_stats)
                track_trajs = arpod_stats{idx}.trackTraj;
                %normalize trajectories (get x axis stuff)
                percent_traj = GraphUtil.normalizeTrajectories(track_trajs);
                kernel = [1,-1];
                
                %get y axis stuff
                errors_y = sum( (arpod_stats{idx}.trackEstTraj - arpod_stats{idx}.trackTraj).^2 );
                errors_y = conv(errors_y,kernel,"same"); %run convolution to get delta

                percent_trajs = [percent_trajs, percent_traj];
                errors_ys = [errors_ys, errors_y];
            end
            GraphUtil.graph2dPt(percent_traj, errors_y,'r',first_graph);
            hold off
        end

        function void = graphFuelBar(list_name_stats)
            %this will create a boxcat graph of the fuel consumption
            
            list_fuel_stats = [];
            for idx = 1:length(list_name_stats)
                load(list_name_stats(idx), "-mat", "savedstats");
                arpod_stats = savedstats;
                list_fuel_stats = [list_fuel_stats; zeros(1,length(arpod_stats))];
            end
            for yolo_idx = 1:length(list_name_stats)
                load(list_name_stats(yolo_idx), "-mat", "savedstats");
                arpod_stats = savedstats;
                for idx_idx = 1:length(arpod_stats)
                    list_fuel_stats(yolo_idx,idx_idx) = arpod_stats{idx_idx}.trackFuelConsumption(end) * ARPOD_Benchmark.m_c;
                end
            end
            boxchart(list_fuel_stats.');
        end

        function void = graphFuelBarPhase(list_name_stats, phase)
            %this will create a boxcat graph of the fuel consumption
            
            list_fuel_stats = [];
            for idx = 1:length(list_name_stats)
                load(list_name_stats(idx), "-mat", "savedstats");
                arpod_stats = savedstats;
                list_fuel_stats = [list_fuel_stats; zeros(1,length(arpod_stats))];
            end
            for yolo_idx = 1:length(list_name_stats)
                load(list_name_stats(yolo_idx), "-mat", "savedstats");
                arpod_stats = savedstats;
                for idx_idx = 1:length(arpod_stats)
                    idxPhase = find(arpod_stats{idx_idx}.trackPhase==phase);
                    fuel_phase = arpod_stats{idx_idx}.trackFuelConsumption(idxPhase);
                    if length(fuel_phase) ~= 0
                        list_fuel_stats(yolo_idx,idx_idx) = fuel_phase(end) * ARPOD_Benchmark.m_c;
                    end
                end
            end
            boxchart(list_fuel_stats.');
        end

        function void = graphRuntimeBar(list_name_stats)
            list_runtime_stats = [];
            for idx = 1:length(list_name_stats)
                load(list_name_stats(idx), "-mat", "savedstats");
                arpod_stats = savedstats;
                runtime_stats = [];
                for idx_idx = 1:length(arpod_stats)
                    runtime_stats = [runtime_stats, arpod_stats{idx_idx}.estimation_ts*1000];
                end
                disp(list_name_stats(idx) + "," + mean(runtime_stats*1000));
                isout = isoutlier(runtime_stats, 'quartiles');
                runtime_stats(isout) = NaN;
                if length(list_runtime_stats) > 0
                    if length(list_runtime_stats) > length(runtime_stats)
                        dist = length(list_runtime_stats) - length(runtime_stats);
                        runtime_stats = [runtime_stats, NaN(1,dist)];
                    elseif length(list_runtime_stats) < length(runtime_stats)
                        [row,col] = size(list_runtime_stats);
                        list_runtime_stats = [list_runtime_stats, NaN(row,length(runtime_stats)-col)];
                    end
                end
                list_runtime_stats = [list_runtime_stats; runtime_stats];
            end
            boxchart(list_runtime_stats.');
        end
        function void = graphRuntimeBarPhase(list_name_stats, phase)
            list_runtime_stats = [];
            for idx = 1:length(list_name_stats)
                load(list_name_stats(idx), "-mat", "savedstats");
                arpod_stats = savedstats;

                runtime_stats = [];
                for idx_idx = 1:length(arpod_stats)
                    idxPhase = find(arpod_stats{idx_idx}.trackPhase==phase);
                    time_phase = arpod_stats{idx_idx}.estimation_ts(idxPhase);
                    runtime_stats = [runtime_stats, time_phase*1000];
                end
                isout = isoutlier(runtime_stats, 'quartiles');
                runtime_stats(isout) = NaN;
                if length(list_runtime_stats) > 0
                    if length(list_runtime_stats) > length(runtime_stats)
                        dist = length(list_runtime_stats) - length(runtime_stats);
                        runtime_stats = [runtime_stats, NaN(1,dist)];
                    elseif length(list_runtime_stats) < length(runtime_stats)
                        [row,col] = size(list_runtime_stats);
                        list_runtime_stats = [list_runtime_stats, NaN(row,length(runtime_stats)-col)];
                    end
                end
                list_runtime_stats = [list_runtime_stats; runtime_stats];
            end
            boxchart(list_runtime_stats.');
        end

        function rmsef = calculateRMSE(est,gt)
            % mxn where m = features and n = # of data points
            %gt = ground truth
            %est = estimated 
            rmsef = rmse(est.',gt.');
        end
        function void = DockSuccessBar(list_name_stats, pos_bool)

            ys = [];
            for idx__ = 1:length(list_name_stats)
                load(list_name_stats(idx__), "-mat", "savedstats");
                arpod_stats = savedstats;

                error_ys = [];
                for idx_idx = 1:length(arpod_stats)
                    est = arpod_stats{idx_idx}.trackEstTraj;
                    gt = arpod_stats{idx_idx}.trackTraj;
                    if pos_bool
                        gt = gt(1:3,:);
                        est = est(1:3,:);
                    else
                        gt = gt(4:6,:);
                        est = est(4:6,:);
                    end
                    % gt = vecnorm(gt);
                    % est = vecnorm(est);
                    sos = sqrt(mean(vecnorm(gt - est)).^2);

                    % sos = GraphUtil.calculateRMSE(est, gt).';
                    error_ys = [error_ys, sos];
                end
                isout = isoutlier(error_ys, 'quartiles');
                error_ys(isout) = NaN;
                if length(ys) > 0
                    if length(ys) > length(error_ys)
                        dist = length(ys) - length(error_ys);
                        error_ys = [error_ys, NaN(1,dist)];
                    elseif length(ys) < length(error_ys)
                        [row,col] = size(ys);
                        ys = [ys, NaN(row,length(error_ys)-col)];
                    end
                end
                ys = [ys; error_ys];
            end
            boxchart(ys.');
        end
        function void = ErrorPhaseBar(list_name_stats, phase, pos_bool)

            ys = [];
            for idx__ = 1:length(list_name_stats)
                load(list_name_stats(idx__), "-mat", "savedstats");
                arpod_stats = savedstats;

                error_ys = [];
                for idx_idx = 1:length(arpod_stats)
                    idxPhase = find(arpod_stats{idx_idx}.trackPhase==phase);

                    est = arpod_stats{idx_idx}.trackEstTraj;
                    est = est(:,idxPhase);
                    gt = arpod_stats{idx_idx}.trackTraj;
                    gt = gt(:,idxPhase);
                    if pos_bool
                        gt = gt(1:3,:);
                        est = est(1:3,:);
                    else
                        gt = gt(4:6,:);
                        est = est(4:6,:);
                    end

                    % gt = vecnorm(gt);
                    % est = vecnorm(est);

                    % sos = GraphUtil.calculateRMSE(est, gt).';

                    % err = vecnorm(gt - est);
                    sos = sqrt(mean(vecnorm(gt - est)).^2);
                    error_ys = [error_ys, sos];
                end
                isout = isoutlier(error_ys, 'quartiles');
                error_ys(isout) = NaN;
                if length(ys) > 0
                    if length(ys) > length(error_ys)
                        dist = length(ys) - length(error_ys);
                        error_ys = [error_ys, NaN(1,dist)];
                    elseif length(ys) < length(error_ys)
                        [row,col] = size(ys);
                        ys = [ys, NaN(row,length(error_ys)-col)];
                    end
                end
                ys = [ys; error_ys];
            end
            boxchart(ys.');
        end
        %%%%% END OF UTIL CODE %%%%%
    end
end