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
    
    xeq = zeros(length(x0),1);
    xeq(2) = 10;
    xeq(end) = x0(end);
    ueq = [xeq(end)*consts.g/consts.gamma; 0];

    A = jacobian(dx,x);
    B = jacobian(dx,u);
    A = subs(A,x,xeq);
    B = subs(B,x,xeq);
    A = double(subs(A,u,ueq));
    B = double(subs(B,u,ueq));

    Q = diag([1, 1, 1, 1, 1, 10, 1, 1, 0]);
    R = eye(length(u));
    [K, P] = lqr(A, B, Q, R);
    ctrl.K = K;
    %{
    lqr_ctrl = matlabFunction(-K*(x - xeq) + ueq, 'vars', x);
    ctrl.func = @(s) lqr_ctrl(s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9));
    %}

    [time, xtraj, utraj] = generate_trajectory(x0, xeq, 10, 0.01, consts);
    ctrl.xtraj = xtraj;
    ctrl.utraj = utraj;
    lqr_ctrl = matlabFunction(-K*x, 'vars', x);
    ctrl.func = @(s) lqr_ctrl(s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9));

    %{
    V = (x - xeq).'*P*(x - xeq);
    LfV = jacobian(V,x)*fx;
    LgV = jacobian(V,x)*gx;
    sontag = -(LfV + sqrt(LfV^2 + (LgV*LgV.')^2))*LgV/(LgV*LgV.').*(LgV ~= 0);
    sontag_func = matlabFunction(sontag,'vars',x);
    ctrl.func = @(s) sontag_func(s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8),s(9))';
    ctrl.xeq = xeq;
    %}
end
