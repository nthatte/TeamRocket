% Function to setup and perform any one-time computations
% Input parameters
%   x - state of the rocket, x=[y,z,th,psi,dy,dz,dth,dpsi,m]^T
%   consts - structure that contains various system constants
% Output parameters
%   ctrl  -  any student defined control parameters
function ctrl = student_setup(x0, consts)
    % Replace line below with one time computation if needed.
    % Ex: ctrl.K = lqr(A, B, Q, R) ;
    syms x y z theta phi dy dz dtheta dphi m ft tau
    
    M = diag([m, m, consts.J, consts.JT, 1]);
    x = [y; z; theta; phi; dy; dz; dtheta; dphi; m];
    fx = [dy; dz; dtheta; dphi; inv(M)*[0; -m*consts.g; 0; 0; 0]];
    gx = [0, 0; 
          0, 0; 
          0, 0; 
          0, 0; 
          inv(M)*[-consts.gamma*sin(phi + theta),  0;
                   consts.gamma*cos(phi + theta),  0; 
                  -consts.L*consts.gamma*sin(phi), 0;
                   0, 1;
                  -1, 0]];
    u = [ft; tau];

    dx = fx + gx*u;
    
    ctrl.xeq = zeros(length(x0),1);
    ctrl.xeq(2) = 10;
    ctrl.xeq(end) = x0(end);
    ueq = [ctrl.xeq(end)*consts.g/consts.gamma; 0];

    dt = 0.1;
    [time, xtraj, utraj] = generate_trajectory(x0, ctrl.xeq, 5, dt, consts);
    ctrl.num_pts = length(time);
    ctrl.xtraj = xtraj;
    ctrl.utraj = utraj;
    new_xtraj = inf(size(xtraj));
    new_utraj = inf(size(utraj));
    new_xtraj(1,:) = xtraj(1,:);
    new_utraj(1,:) = utraj(1,:);

    ctrl.num_pts = length(time);
    ctrl.func = cell(ctrl.num_pts, 1);
    ctrl.K_t  = cell(ctrl.num_pts, 1);

    %lqr linearized about trajectory
    Q = diag([1, 1, 1, 0, 1, 10, 1, 0, 0]);
    R = eye(length(u));
    A_x = matlabFunction(jacobian(dx,x), 'vars', [x; u]);
    A_x = @(s) A_x(s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9), s(10), s(11));
    B_x = matlabFunction(jacobian(dx,u), 'vars', [x; u]);
    B_x = @(s) B_x(s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9), s(10), s(11));


    traj_error_thresh = 5;
    num_iter = 25;
    magic_factor = 0.9;
    odeopts = odeset;
    for j = 1:num_iter
        for i = 1:ctrl.num_pts
            A = A_x([xtraj(i,:)'; utraj(i,:)']);
            B = B_x([xtraj(i,:)'; utraj(i,:)']);
            ctrl.K_t{i} = lqr(A, B, Q, R);
            if mod(i,100) == 0
                i
            end
            alpha_x = @(s) -ctrl.K_t{i}*(s - ctrl.xtraj(i,:)') + ctrl.utraj(i,:)';
            ustep = alpha_x(new_xtraj(i,:)');
            new_utraj(i,:) = ctrl.utraj(i,:) + magic_factor*(ustep' - ctrl.utraj(i,:));
            if i ~= ctrl.num_pts
                [~, xstep] = ode45(@odefun_rocket, [0 dt], new_xtraj(i,:), odeopts, consts, alpha_x) ;
                %{
                disp('xstep')
                disp(xstep(end,:))
                disp('xtraj')
                disp(ctrl.xtraj(i+1,:))
                %}
                new_xtraj(i+1,:) = ctrl.xtraj(i+1,:) + magic_factor*(xstep(end,:) - ctrl.xtraj(i+1,:));
            end
        end

        figure(100)
        subplot(211)
        plot(time, ctrl.xtraj, '-', time, new_xtraj, '--')
        subplot(212)
        plot(time, ctrl.utraj, '-', time, new_utraj, '--')
        norm(ctrl.xtraj - new_xtraj)
        if norm(ctrl.xtraj - new_xtraj) < traj_error_thresh
            break;
        end
            
        ctrl.xtraj = new_xtraj;
        ctrl.utraj = new_utraj;

    end
end

function [dx u] = odefun_rocket(t, x, consts, alpha_x)
    y = x(1) ;
    z = x(2) ;
    th = x(3) ; 
    psi = x(4) ;
    
    dy = x(5) ;
    dz = x(6) ;
    dth = x(7) ;
    dpsi = x(8) ;
    
    m = x(9) ;
    
    f_vec = [   dy
               dz
              dth
             dpsi
                0
               -consts.g
                0
                0
                0] ;
            
    g_vec = [0,    0 ;
             0,    0 ;
             0,    0 ;
             0,    0 ;
            -consts.gamma*sin(psi+th)/m,    0 ;
             consts.gamma*cos(psi+th)/m,    0 ;
             -consts.L*consts.gamma*sin(psi)/consts.J,    0 ;
             0, 1/consts.JT ;
            -1,    0] ;
    
    % call student controller
    u = alpha_x(x);

    % Check if fuel is over
    if(m <= consts.m_nofuel)
        u(1) = 0 ;
    end
   
    % Thrust / Torque saturations
    u(1) = min(max(u(1), 0), consts.max.fT) ;
    u(2) = min(max(u(2), -consts.max.tau), consts.max.tau) ;
    
    dx = f_vec + g_vec*u ;
end
