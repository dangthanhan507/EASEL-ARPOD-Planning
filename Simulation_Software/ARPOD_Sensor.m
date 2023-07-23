classdef ARPOD_Sensor
    methods (Static)
        function z_t = measure(state) %add flags for measurement
            x = state(1,:);
            y = state(2,:);
            z = state(3,:);


            norm = sqrt(x*x+y*y+z*z);
            e1 = atan(y/x);
            e2 = asin(z/norm);
            e3 = norm;
            z_t = [e1;e2;e3];
        end
        function z_t = Tmeasure(state) %TensCalc Supported Measurement function
            x = state(1,:);
            y = state(2,:);
            z = state(3,:);


            norm = sqrt( (x.*x) + (y.*y) + (z.*z) );

            %e1 = atan(y ./ );
            e1 = Tatan2(y,x);
            %e1 = y .* x;
            %e2 = asin(z/norm);
            % asin not supported, so do it in terms of atan2
            z_norm = z ./ norm;
            e2 = y .* x;
            %e2 = Tatan2(z_norm, sqrt((1+z_norm)*(1-z_norm)));
            e3 = norm;
            z_t = [e1;e2;e3];
        end
        function dz_t = measureJacobian(state) %take jacobian of measurement function (used on EKF)
            x = state(1);
            y = state(2);
            z = state(3);

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

            dz_t = jacobian;
        end
        function z_t = sense(state, noise_model, phase)
            if phase == 1
                z_t = ARPOD_Sensor.measure(state);
                z_t = z_t + noise_model();
                z_t = z_t(1:2,:);
            else
                z_t = ARPOD_Sensor.measure(state);
                z_t = z_t + noise_model();
            end
        end
        function z_t = fullSense(state, noise_model)
            z_t = state + noise_model();
        end
    end
end