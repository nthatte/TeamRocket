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
    index = floor(t/0.1)+1;
    if index < ctrl.num_pts
        u = -ctrl.K_t{index}*[x; 1] + ctrl.utraj_des(index,:)';
    else
        u = -ctrl.K_t{end}*[x; 1] + ctrl.utraj_des(end,:)';
    end
end

