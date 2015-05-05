%% Runs sim_rocket
clearvars; close all;

% System Constants
consts = get_consts() ;

% Default  ICs
x0 = [10; 150; 0; 0;
          0; 0; 0; 0;
          consts.m_nofuel+0.7*consts.max.m_fuel] ;
      
% Desired state
xeq = zeros(length(x0),1);
xeq(2) = 10;
xeq(end) = x0(end);


% New ICs
x_offset = [1;              % y
            10;             % z
            2*pi/180;       % theta
            10*pi/180;      % phi
            0.1;            % dy
            -0.1;           % dz
            1*pi/180;       % dtheta
            2*pi/180;       % dphi
            0];             % m

[time, x_t, u_t] = generate_trajectory(x0, xeq, 5, 0.1, consts);

figure(1)
plot(time, x_t)
figure(2)
plot(time, u_t)
