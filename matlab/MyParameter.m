%% Parameter setting %%

% % % % % % Confirm your parameters % % % % % % %

param.nx = 750;     % number of voxels
param.ny = 750;
param.nz = 450;

% MODIFIED
param.sx = 150;     % [mm]
param.sy = 150;     % [mm]
param.sz = 90;      % [mm]

% The real detector panel pixel density (number of pixels)
param.nu = 1500;    % number of pixels
param.nv = 1628;

% MODIFIED
param.su = 147;	    % u-directional detector size [mm]
param.sv = 159.544; % v-directional detector size [mm] 

% X-ray source and detector setting (MODIFIED)
param.DSD = 658.45; % source-to-detector distance [mm]
param.DSO = 409.70;	% source-to-object distance [mm]

% angle setting
param.dir = 1;              % gantry rotating direction (clock wise/ counter clockwise)
param.dang = 360/706;       % angular step size (deg)

param.deg = 0:param.dang:360-param.dang;    % MODIFIED
param.deg = param.deg*param.dir;
param.nProj = length(param.deg);
% disp(param.nProj);

param.parker = 0;   % data with 360 deg -> param.parker = 0 , data less than 360 deg -> param.parker=1 

% % % % % % Confirm your parameters % % % % % % %
 
% filter='ram-lak','cosine', 'hamming', 'hann' 
param.filter = 'ram-lak';       % high pass filter

param.dx = param.sx/param.nx;   % single voxel size
param.dy = param.sy/param.ny;
param.dz = param.sz/param.nz;
param.du = param.su/param.nu;
param.dv = param.sv/param.nv;

param.off_z = -40;  % [mm] 

param.off_u = -41.633;  % [mm] detector rotation shift (in real size)
param.off_v = -74.662;  % [mm]

% % % Geometry calculation % % %
param.xs = (-(param.nx-1)/2:1:(param.nx-1)/2)*param.dx;
param.ys = (-(param.ny-1)/2:1:(param.ny-1)/2)*param.dy;
param.zs = (-(param.nz-1)/2:1:(param.nz-1)/2)*param.dz + param.off_z;

param.us = (-(param.nu-1)/2:1:(param.nu-1)/2)*param.du + param.off_u;
param.vs = (-(param.nv-1)/2:1:(param.nv-1)/2)*param.dv + param.off_v;

param.interptype = 'nearest'; % 'linear', 'nearest'

% % % % % % Confirm your parameters % % % % % % %
% Only for Matlab version above 2013b with parallel computing toolbox: Speed dependent on your GPU
% You don't need to install CUDA, but install latest graphics driver.
% only Nvidia GPU cards can use this. otherwise please "param.gpu=0"
% This option is semi-GPU code using built-in Matlab GPU functions: several times faster than CPU
param.gpu = 1;

















