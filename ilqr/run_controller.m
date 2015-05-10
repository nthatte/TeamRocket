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
x_offset = [10;              % y
            150;             % z
            60*pi/180;       % theta
            00*pi/180;      % phi
            1;            % dy
            10;           % dz
            00*pi/180;       % dtheta
            0*pi/180;       % dphi
            0];             % m
        

% Run model (ICs close to xeq) (easy)
sim_rocket(xeq + x_offset)

% Run model (ICs = x0) (medium)
%sim_rocket(x0)

% Run model (ICs = x0_difficult)
% difficulty_factor = 10;
% x0 + x_offset*difficulty_factor;
% sim_rocket()


% Post-processing
tilefigs;
