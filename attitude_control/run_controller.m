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
xeq(2) = 1;
xeq(end) = x0(end);


% New ICs
x_offset = [50;              % y
            1500;             % z
            0.2;       % theta
            0;      % phi
            0;            % dy
            0;           % dz
            0;       % dtheta
            0;       % dphi
            consts.m_nofuel+1*consts.max.m_fuel];             % m
        

x_offset = [ 10 150 0.1 0   0 0 0 0  consts.m_nofuel+1*consts.max.m_fuel]';
x_offset = [ 100 25 0 0   0 0 0 0  consts.m_nofuel+1*consts.max.m_fuel]';
x_offset = [ 100 500 pi/2 0   0 0 0 0  consts.m_nofuel+1*consts.max.m_fuel]';
%x_offset = [ 10 1500 179*pi/180 0   0 0 0 0  consts.m_nofuel+1*consts.max.m_fuel]';

% x_offset = [ 60 480 0.2 0   -40 -10 0 0  consts.m_nofuel+1*consts.max.m_fuel]';
% Run model (ICs close to xeq) (easy)
sim_rocket(x_offset)

% Run model (ICs = x0) (medium)
%sim_rocket(x0)

% Run model (ICs = x0_difficult)
% difficulty_factor = 10;
% x0 + x_offset*difficulty_factor;
% sim_rocket()


% Post-processing
tilefigs;
