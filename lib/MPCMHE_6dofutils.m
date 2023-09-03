classdef MPCMHE_6dofutils
    methods (Static)
        function dynamicConstraints = createFixedBodyDynamicConstraints(A, B, Tx0, Tx, uback, uforward, disturbance, disturbType)

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

            %setup 
            dynamicConstraints = (Tx == A*xk + B*uk)
        end
    end
end