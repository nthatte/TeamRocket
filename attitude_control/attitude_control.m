% Student supplied function to compute the control input at each instant of time
% Input parameters:
%   t  -  Time (in seconds) since start of simulation
%   x - state of the rocket, x=[y,z,th,psi,dy,dz,dth,dpsi,m]^T
%   consts - structure that contains various system constants
%   ctrl  -  any student defined control parameters
% Output parameters:
%   u  -  [thrust; torque] - two inputs to the rocket
function u = attitude_control(t, x, consts, ctrl)
 
    u = [0;0];



            u(1) = x(end)*consts.g/consts.gamma/cos(x(3)+x(4));
            u(1) = x(end)*consts.g/consts.gamma*2;
            if abs(x(3))>pi/2 
                wn1 = ctrl.wn1h;
            else
                wn1 = ctrl.wn1l;
            end
            wn2 = ctrl.wn2;
            phi_d = (-x(3)*wn1^2 - x(7)*2*wn1)*consts.J/(consts.L*consts.gamma)/u(1);
            phi_d = asin(-phi_d);
            u(2) = ((-x(4)+phi_d)*wn2^2 - x(8)*2*wn2)*consts.JT;
 


            

        
end
