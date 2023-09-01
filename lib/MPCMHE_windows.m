classdef MPCMHE_windows
    properties
        x0

        disturbType

        window_mhestates
        window_mhecontrols
        window_measurements
        window_controlDisturbances
        window_measError

        window_mpcstates
        window_mpccontrols
    end

    methods
        function obj = init(obj, x0, window_states, window_measurements, window_controls, forwardT, disturbType)
            obj.x0 = x0;

            [control_dim, backwardT] = size(window_controls);
            [state_dim, backwardT] = size(window_states);
            [meas_dim, backwardT] = size(window_measurements);


            %setup mhe windows 
            obj.window_mhestates = window_states;
            obj.window_mhecontrols = window_controls;
            obj.window_measurements = window_measurements;

            obj.disturbType = disturbType;
            %setup mhe disturbance window
            if disturbType == 0
                obj.window_controlDisturbances = zeros(control_dim,1);
            elseif disturbType == 1
                obj.window_controlDisturbances = zeros(control_dim,backwardT + forwardT);
            elseif disturbType == 2
                obj.window_controlDisturbances = zeros(control_dim,1);
            elseif disturbType == 3
                obj.window_controlDisturbances = zeros(control_dim,control_dim);
            else
                obj.window_controlDisturbances = [eye(control_dim,control_dim). zeros(control_dim,1)];
            end

            obj.window_measError = zeros(meas_dim, backwardT);

            %setup mpc windows
            obj.window_mpccontrols = zeros(control_dim, forwardT);
            obj.window_mpcstates = zeros(state_dim, forwardT);
        end

        function obj = setWindows(mheDs, vs, mheXs, mpcXs, mpcUs)
            obj.window_controlDisturbances = mheDs;
            obj.window_measError   = vs;
            obj.window_mhestates   = mheXs;
            obj.window_mpcstates   = mpcXs;
            obj.window_mpccontrols = mpcUs;
        end
        function obj = shiftWindows(obj, measurement, control)
            obj.x0 = obj.window_mhestates(:,1);

            %shift and use obj.window_mpcstates as a "predicted variable" for warm-start
            obj.window_mhestates = [obj.window_mhestates(:,2:end), obj.window_mpcstates(:,1)];
            %propagate zero control input thrust as warm-start for last element
            obj.window_mpcstates = [obj.window_mpcstates(:,2:end), obj.window_mpcstates(:,end)];

            %assume last element no control disturbance for warm-start

            %assume no measError on last element for warm-start
            [dim,n] = size(obj.window_measError);
            obj.window_measError = [obj.window_measError(:,2:end), zeros(dim,1)];

            %push mhecontrols with new control input
            %assume no new thrust at tail end of mpc
            [dim,n] = size(obj.window_mpccontrols);
            obj.window_mhecontrols = [obj.window_mhecontrols(:,2:end), control];
            obj.window_mpccontrols = [obj.window_mpccontrols(:,2:end), zeros(dim,1)];

            %push new measurement into mhe measurements
            obj.window_measurements = [obj.window_measurements(:,2:end), measurement];

            if obj.disturbType == 0
                [control_dim, backwardT] = size(obj.window_mhecontrols);
                obj.window_controlDisturbances = zeros(control_dim,1);
            elseif obj.disturbType == 1
                obj.window_controlDisturbances = [obj.window_controlDisturbances(:,2:end), obj.window_controlDisturbances(:,end)];
            elseif obj.disturbType == 2
                %nothin for TI
            elseif obj.disturbType == 3
                %nothin for TI
            else
                %nothin for TI
            end
        end

        function [mheXs, mheVs, mpcXs, mpcUs] = getOptWindows(obj)
            mheXs = obj.window_mhestates;
            mheVs = obj.window_measError;
            mpcXs = obj.window_mpcstates;
            mpcUs = obj.window_mpccontrols;
        end
    end
end