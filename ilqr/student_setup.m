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
    
    %derive dynamics
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
    num_states = length(x);
    num_inputs = length(u);

    dx = fx + gx*u;

    %discrete time dynamics model
    dt = 0.01;
    fxu_discrete = x + dt*dx;
    fxu_discrete = matlabFunction(fxu_discrete, 'vars', [x; u]);
    fxu_discrete = @(s) fxu_discrete(s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9), s(10), s(11));

    %define equilibrium point and trajectroy from initial condition to eq
    ctrl.xeq = zeros(num_states,1);
    ctrl.xeq(2) = 9;
    ctrl.xeq(end) = x0(end);
    ctrl.ueq = [x0(end)*consts.g/consts.gamma, 0]';

    [time, xtraj, utraj] = generate_trajectory(x0, ctrl.xeq, 5, dt, consts);
    ctrl.num_pts = length(time);
    ctrl.xtraj_des = xtraj;
    ctrl.utraj_des = utraj;
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
    A_cont = matlabFunction(jacobian(dx,x), 'vars', [x; u]);
    A_cont = @(s) A_cont(s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9), s(10), s(11));
    B_cont = matlabFunction(jacobian(dx,u), 'vars', [x; u]);
    B_cont = @(s) B_cont(s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9), s(10), s(11));

    %compute infinite time horizon control to use after last timestep
    A = A_cont([ctrl.xeq; ctrl.ueq]);
    B = B_cont([ctrl.xeq; ctrl.ueq]);
    Q = diag([1, 1, 10, 1, 1, 10, 1, 1, 0]);
    R = eye(num_inputs);
    ctrl.K_end = lqr(A, B, Q, R);

    traj_error_thresh = .001;
    num_iter = 10;
    magic_factor_mult = 0.1/num_iter;
    %magic_factor_mult = 0.1;
    odeopts = odeset;
    for j = 1:num_iter
        j
        magic_factor = magic_factor_mult*j
        %backwards pass
        for i = ctrl.num_pts:-1:1
            Q = diag([1, 1, 10, 1, 1, 10, 10, 1, 0])/dt;
            q = ctrl.xtraj_des(i,:)*Q;
            Q = [Q, -q'; -q, 1];
            R = eye(num_inputs)/dt;

            A_discrete = eye(num_states) + A_cont([ctrl.xtraj(i,:)'; ctrl.utraj(i,:)'])*dt;
            B_discrete = B_cont([ctrl.xtraj(i,:)'; ctrl.utraj(i,:)'])*dt;
            A = [A_discrete, (eye(num_states) - A_discrete)*ctrl.xtraj(i,:)' ...
                    + B_discrete*(ctrl.utraj_des(i,:)' - ctrl.utraj(i,:)');
                 zeros(1, num_states), 1];
            B = [B_discrete; zeros(1, num_inputs)];

            if i == ctrl.num_pts
                Ps = Q;
                K = zeros(num_inputs, num_states+1);
            else
                K = (R + B'*Ps*B)\(B'*Ps*A);
            end
            tmp = A-B*K;
            Ps = Q + K'*R*K + tmp'*Ps*tmp;
            ctrl.K_t{i} = K;

            %{
            if mod(i,100) == 0
                i
            end
            %}
        end

        %forward rollout
        for i = 1:ctrl.num_pts-1
            new_utraj(i,:) = (ctrl.K_t{i}*[new_xtraj(i,:)'; 1] ...
                + ctrl.utraj_des(i,:)')';
            new_utraj(i,:) = ctrl.utraj_des(i,:) ...
                + magic_factor*(new_utraj(i,:) - ctrl.utraj_des(i,:));

            new_xtraj(i+1,:) = fxu_discrete([new_xtraj(i,:)';
                new_utraj(i,:)']);
            new_xtraj(i+1,:) = ctrl.xtraj_des(i+1,:) ...
                + magic_factor*(new_xtraj(i+1,:) - ctrl.xtraj_des(i+1,:));
        end
        new_utraj(end,:) = (ctrl.K_t{end}*[new_xtraj(end,:)'; 1] ...
            + ctrl.utraj_des(end,:)')';
        new_utraj(end,:) = ctrl.utraj_des(end,:) ...
            + magic_factor*(new_utraj(end,:) - ctrl.utraj_des(end,:));

        %{
        [~, new_xtraj] = ode45(@odefun_rocket, time, new_xtraj(1,:)', odeopts, consts, ctrl);
        for i = 1:ctrl.num_pts
            new_utraj(i,:) = (ctrl.K_t{i}*[new_xtraj(i,:)'; 1] + ctrl.utraj_des(i,:)')';
        end
        
        % apply magic factor
        new_utraj = ctrl.utraj_des + magic_factor*(new_utraj - ctrl.utraj_des);
        new_xtraj = ctrl.xtraj_des + magic_factor*(new_xtraj - ctrl.xtraj_des);
        %}

        %plot old, new, and desired trajectories
        figure(100)
        subplot(211)
        plot(time, ctrl.xtraj)
        hold on
        plot(time, new_xtraj,'--')
        plot(time, ctrl.xtraj_des,'-.')
        hold off
        subplot(212)
        plot(time, ctrl.utraj)
        hold on
        plot(time, new_utraj,'--')
        plot(time, ctrl.utraj_des,'-.')
        hold off

        pause(0.1)

        %check for convergence else perform another iteration
        norm(ctrl.xtraj - new_xtraj)
        if norm(ctrl.xtraj - new_xtraj) < traj_error_thresh
            break;
        end
            
        ctrl.xtraj = new_xtraj;
        ctrl.utraj = new_utraj;
    end
end

function [dx u] = odefun_rocket(t, x, consts, ctrl)
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
    index = floor(t/0.01)+1;
    u = -ctrl.K_t{index}*[x; 1] + ctrl.utraj_des(index,:)';
    %u = student_controller(t, x, consts, ctrl) ;

    % Check if fuel is over
    if(m <= consts.m_nofuel)
        u(1) = 0 ;
    end
   
    % Thrust / Torque saturations
    u(1) = min(max(u(1), 0), consts.max.fT) ;
    u(2) = min(max(u(2), -consts.max.tau), consts.max.tau) ;
    
    dx = f_vec + g_vec*u ;
end
