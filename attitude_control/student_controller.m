% Student supplied function to compute the control input at each instant of time
% Input parameters:
%   t  -  Time (in seconds) since start of simulation
%   x - state of the rocket, x=[y,z,th,psi,dy,dz,dth,dpsi,m]^T
%   consts - structure that contains various system constants
%   ctrl  -  any student defined control parameters
% Output parameters:
%   u  -  [thrust; torque] - two inputs to the rocket
function u = student_controller(t, x, consts, ctrl)
 
    u = [0;0];
    dt = ctrl.dt;
    index = floor((t-ctrl.conv_t)/dt)+1;

    if t<ctrl.conv_t
        
        u = attitude_control(t, x, consts, ctrl);
       
    elseif t >= (max(ctrl.time))          
         u = -ctrl.K*(x-ctrl.xeq) + [x(end)*consts.g/consts.gamma; 0];
    else
         u =  -ctrl.K*(x-ctrl.xtraj(index,:)') + [x(end)*consts.g/consts.gamma; 0];   
    end

end
