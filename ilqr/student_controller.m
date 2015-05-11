% Student supplied function to compute the control input at each instant of time
% Input parameters:
%   t  -  Time (in seconds) since start of simulation
%   x - state of the rocket, x=[y,z,th,psi,dy,dz,dth,dpsi,m]^T
%   consts - structure that contains various system constants
%   ctrl  -  any student defined control parameters
% Output parameters:
%   u  -  [thrust; torque] - two inputs to the rocket
function u = student_controller(t, x, consts, ctrl)
    % Replace line below with your controller
    % Ex: u = -ctrl.K*x ;
    
    %u = ctrl.func(x);
    index = floor(t/0.01)+1;
    if index <= ctrl.num_pts - 1000
        u = -ctrl.K_t{index}*[x; 1] + ctrl.utraj_des(index,:)';
        %u = -ctrl.K_t{index}*[x; 1] + [x(end)*consts.g/consts.gamma; 0];
    else
        %u = ctrl.K_t{end}*[x; 1] + ctrl.utraj_des(end,:)';
        %u = [0; 0]; 
        %u = -ctrl.K_t{end}*[x; 1] + ctrl.utraj_des(end,:)';
        %u = -ctrl.K_end*(x-ctrl.xeq) + ctrl.ueq;
        %u = -ctrl.K_t{end}*[x; 1] + [x(end)*consts.g/consts.gamma; 0];
        u = -ctrl.K_end*(x) + [x(end)*consts.g/consts.gamma; 0];
    end
end
