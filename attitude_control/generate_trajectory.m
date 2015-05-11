function [time, x_t, u_t] = generate_trajectory(x0, xf, traj_time_scale, dt, consts);

    ueq = [x0(end)*consts.g/consts.gamma, 0];

    g = 9.81;
    height = x0(2);
    T = sqrt(2*height/g)*traj_time_scale;

    num_coeff = 6;
    A = [1, 0,   0,     0,      0,      0;
         1, T, T^2,   T^3,    T^4,    T^5;
         0, 1,   0,     0,      0,      0;
         0, 1, 2*T, 3*T^2,  4*T^3,  5*T^4;
         0, 0,   2,     0,      0,      0; 
         0, 0,   2,   6*T, 12*T^2, 20*T^3]; 

    n = 4;
    initial_pt  = x0(1:n)';
    initial_vel = x0(n+1:2*n)';
    final_pt  = xf(1:n)';
    final_vel = xf(n+1:2*n)';

    end_points = [initial_pt;
                  final_pt;
                  initial_vel;
                  final_vel;
                  zeros(1,n);
                  zeros(1,n)];

    coeff_x = A\end_points;
    coeff_dx = [coeff_x(2,:); 2*coeff_x(3,:); 3*coeff_x(4,:); ...
        4*coeff_x(5,:); 5*coeff_x(6,:)];

    time = 0:dt:T;
    x_t = zeros(length(time), 2*n+1);

    %position traj
    for i = 1:n;
        x_t(:,i) = polyval(flipud(coeff_x(:,i)), time);
    end
    %vel traj
    for i = 1:n;
        x_t(:,n+i) = polyval(flipud(coeff_dx(:,i)), time);
    end
    %mass traj
    x_t(:,end) = x0(end)*ones(length(time), 1);

    u_t = repmat(ueq,length(time),1);
    
    %if z trajectory goes below zero then refind z subject to constraints
    if any(x_t(:,2) < 0)
        %end point equality constraints
        Aeq = A(1:4,:);
        beq = end_points(1:4,2);

        %height inequality constraints: above quadratic from start/2 to end height
        Aineq = -[ones(size(time')), time', time'.^2, time'.^3, time'.^4, time'.^5];
        Ah = [0, 0, 1; T^2, T, 1; 2*T, 1, 0];
        bh = [x0(2)/2; xf(2); xf(6)];
        bound_coeffs = Ah\bh;
        bineq = -polyval(bound_coeffs, time);
        %bineq = -xf(2)*ones(size(Aineq,1),1);

        %minimum jerk cost
        H = zeros(num_coeff);
        for t = time
            w = [0, 0, 0, 6, 24*t, 60*t^2]';
            H = H + w*w';
        end
        f = zeros(num_coeff,1);

        coeff_z = quadprog(H, f, Aineq, bineq, Aeq, beq);
        coeff_dz = [coeff_z(2); 2*coeff_z(3); 3*coeff_z(4); ...
            4*coeff_z(5); 5*coeff_z(6)];
        x_t(:,2) = polyval(flipud(coeff_z), time);
        x_t(:,n+2) = polyval(flipud(coeff_dz), time);
    end
