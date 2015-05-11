% Function to setup and perform any one-time computations
% Input parameters
%   x - state of the rocket, x=[y,z,th,psi,dy,dz,dth,dpsi,m]^T
%   consts - structure that contains various system constants
% Output parameters
%   ctrl  -  any student defined control parameters
function ctrl = student_setup(x0, consts)

    % attitude control gain
    %it takes 2 second for pi/2 to converge, 10 second for pi to converge(upside down)
       
    if x0(3)<0.3
        ctrl.conv_t = 0.1;
    elseif x0(3)<pi/1.9
        ctrl.conv_t = 2;
    else
        ctrl.conv_t = 5;   
    end
    
    ctrl.wn1h = 1;
    ctrl.wn1l = 0.5;
    ctrl.wn2 = 7;
    dt = 0.01;
    ctrl.dt = dt;
    
    %forward roll out attitude control
    odeopts = odeset;
    [t x] = ode45(@odefun_rocket, [0 ctrl.conv_t], x0, odeopts, consts, ctrl);
    while x(end,3)>0.3
        ctrl.conv_t = ctrl.conv_t*1.5;
        [t x] = ode45(@odefun_rocket, [0 ctrl.conv_t], x0, odeopts, consts, ctrl);
    end
    start_point = x(end,:)';
    ctrl.start_t = t(end);

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
%     fxu_discrete = x + dt*dx;
%     fxu_discrete = matlabFunction(fxu_discrete, 'vars', [x; u]);
%     fxu_discrete = @(s) fxu_discrete(s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9), s(10), s(11));

    %define equilibrium point and trajectroy from initial condition to eq
    ctrl.xeq = zeros(num_states,1);
    ctrl.xeq(2) = 10;
    ctrl.xeq(end) = x0(end);
    ctrl.ueq = [x0(end)*consts.g/consts.gamma, 0]';

    T = 5;
    [time, xtraj, utraj] = generate_trajectory(start_point, ctrl.xeq, T, dt, consts);
    while max(time)>=(60-ctrl.conv_t)
        T = T/2;
        [time, xtraj, utraj] = generate_trajectory(start_point, ctrl.xeq, T, dt, consts);
    end
    while max(time)<40
        T  = T*1.5;
        [time, xtraj, utraj] = generate_trajectory(start_point, ctrl.xeq, T, dt, consts);
    end

    if any(xtraj(:,2) < 0)
        midpoint = (start_point + ctrl.xeq)/2;
        midpoint(5:8) = 0;
        [time1, xtraj1, utraj1] = generate_trajectory(start_point, midpoint, T/2, dt, consts);
        [time2, xtraj2, utraj2] = generate_trajectory(midpoint, ctrl.xeq, T/2, dt, consts);
        time = [time1'; time2(2:end)'+time1(end)'];
        xtraj = [xtraj1; xtraj2(2:end,:)];
        utraj = [utraj1; utraj2(2:end,:)];
    end
   
        
    ctrl.time = time+ctrl.conv_t;
    ctrl.xtraj = xtraj;
    ctrl.utraj = utraj;

    %lqr linearized about trajectory
    A_cont = matlabFunction(jacobian(dx,x), 'vars', [x; u]);
    A_cont = @(s) A_cont(s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9), s(10), s(11));
    B_cont = matlabFunction(jacobian(dx,u), 'vars', [x; u]);
    B_cont = @(s) B_cont(s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9), s(10), s(11));

    %compute infinite time horizon control to use after last timestep
    A = A_cont([ctrl.xeq; ctrl.ueq]);
    B = B_cont([ctrl.xeq; ctrl.ueq]);
    Q = diag([1, 10, 10, 1, 10, 10, 1, 1, 0]);
    R = eye(num_inputs);

    ctrl.K = lqr(A, B, Q, R);

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
    u = attitude_control(t, x, consts, ctrl) ;

    % Check if fuel is over
    if(m <= consts.m_nofuel)
        u(1) = 0 ;
    end
   
    % Thrust / Torque saturations
    u(1) = min(max(u(1), 0), consts.max.fT) ;
    u(2) = min(max(u(2), -consts.max.tau), consts.max.tau) ;
    
    dx = f_vec + g_vec*u ;
end
