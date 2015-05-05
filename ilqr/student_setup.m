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
    
    % normal lqr
    xeq = zeros(length(x0),1);
    xeq(2) = 10;
    xeq(end) = x0(end);
    ueq = [xeq(end)*consts.g/consts.gamma; 0];

    [time, xtraj, utraj] = generate_trajectory(x0, xeq, 10, 0.1, consts);
    ctrl.num_pts = length(time);
    ctrl.xtraj = xtraj;
    ctrl.utraj = utraj;

    Q = diag([1, 1, 1, 0, 1, 10, 1, 0, 0]);
    R = eye(length(u));

    ctrl.num_pts = length(time);
    ctrl.func = cell(ctrl.num_pts, 1);
    ctrl.K_t = cell(ctrl.num_pts, 1);

    %lqr linearized about eq following trajectory
    %{
    A = jacobian(dx,x);
    B = jacobian(dx,u);
    A = subs(A,x,xeq);
    B = subs(B,x,xeq);
    A = double(subs(A,u,ueq));
    B = double(subs(B,u,ueq));


    K = lqr(A, B, Q, R);
    for i = 1:ctrl.num_pts
        ctrl.K_t{i} = K;
    end
    %}

    %lqr linearized about trajectory
    A_x = matlabFunction(jacobian(dx,x), 'vars', [x; u]);
    A_x = @(s) A_x(s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9), s(10), s(11));
    B_x = matlabFunction(jacobian(dx,u), 'vars', [x; u]);
    B_x = @(s) B_x(s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9), s(10), s(11));

    ctrl.num_pts
    for i = 1:ctrl.num_pts
        A = A_x([xtraj(i,:)'; utraj(i,:)']);
        B = B_x([xtraj(i,:)'; utraj(i,:)']);
        ctrl.K_t{i} = lqr(A, B, Q, R);
        if mod(i,100) == 0
            i
        end
    end
end
