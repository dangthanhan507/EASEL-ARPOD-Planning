classdef MPCMHE_6dofutils
    methods (Static)
        function dynamicConstraints = createFixedBodyLTIDynamics(A, B, Tx0, Tx, uback, uforward, disturbance, disturbType)

            xk = [Tx0, Tx(:,1:end-1)];  
            uk = [uback, uforward];

            [control_dim, n] = size(uk);
            if disturbType == 0
                %leave empty
            elseif disturbType == 1
                %time-varying additive
                for i = 1:n
                    uk(:,i) = uk(:,i) + disturbance(:,i);
                end
            elseif disturbType == 2
                %time-invariant additive
                for i = 1:n
                    uk(:,i) = uk(:,i) + disturbance;
                end
            elseif disturbType == 3
                %time-invariant matrix transform
                for i = 1:n
                    uk(:,i) = disturbance*uk(:,i);
                end
            else
                %time-invariant matrix transform + additive
                dMatrix = disturbance(1:control_dim, 1:control_dim);
                dAdd    = disturbance(1:control_dim, end);
                for i = 1:n
                    uk(:,i) = dMatrix*uk(:,i) + dAdd;
                end
            end

            % 6dof vs normal spacecraft
            %{
                u is 6x1
                [ 1 -1  0  0  0  0] [ux+]
                [ 0  0  1 -1  0  0] [ux-]    [ux]
                [ 0  0  0  0  1 -1] [uy+]  = [uy]
                                    [uy-]    [uz]
                                    [uz+]
                                    [uz-]
            %}
            D = [1, -1,  0,  0,  0,  0;
                 0,  0,  1, -1,  0,  0;
                 0,  0,  0,  0,  1, -1];
            dynamicConstraints = (Tx == A*xk + B*D*uk);
        end
        function opt_properties = setupOptimizationCells(costFunction, dynamicConstraints, sensorConstraints, dMax, uMax, vMax, Tx0, Tx, Td, Tuback, Tyback, Tvback, Tuforward)
            %only difference is we set uMin to zero
            opt_properties = MPCMHE_properties;
            opt_properties = opt_properties.init(costFunction,  {Tuforward}, ...
                                                                {Td, Tvback,Tx0}, ...
                                                                {Tx},...
                                                                {Tuforward <= uMax*Tones(size(Tuforward)), Tuforward >= 0*Tones(size(Tuforward))},...
                                                                {Td<=dMax*Tones(size(Td)), Td>=-dMax*Tones(size(Td)),...
                                                                    Tvback<=vMax*Tones(size(Tvback)), Tvback>=-vMax*Tones(size(Tvback)),...
                                                                    sensorConstraints},...
                                                                {dynamicConstraints},...
                                                                {Tuforward,Td,Tx0,Tx,Tvback},...
                                                                {Tuback,Tyback});
        end
    end
end