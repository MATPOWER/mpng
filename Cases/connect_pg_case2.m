function connect = connect_pg_case2
%% ========================================================================
% Input necessary to connect power and natural gas systems
% =========================================================================
% Authors:  (1) Wilson González Vanegas, M.Sc
%           (2) Sergio García Marín, M.Sc(c)
% (1) Automatic Research Group, Universidad Tecnológica de Pereira.
% (2) GIPEM, Universidad Nacional de Colombia - Sede Manizales.
% =========================================================================
%% ---------------------- Matpower gas case version ----------------------
connect.version = '1';

%% ---------------------- Matpower gas case version ----------------------

connect.power.time =  [24];  % number of hours for each period of time sum()=24
connect.power.ramp_time =  [];
connect.power.UC =  [];
connect.power.demands.pd = [
    0;
    5;
    ];
connect.power.demands.qd = [
    0;
    1;
    ];

%%
connect.power.cost = 1000; % power for non-suplied power demand

% matrix, nt -> colums, nz -> rows
connect.power.sr = [0];
%   genid       maxenergy
connect.power.energy = ...
    [
    1   24*100;
%     2   24*100;
    ];%


%% ---------------------- Matpower gas case version ----------------------
%                   idcomp      bus
connect.interc.comp =   [];     % compressor connected to bus
% connect.interc.comp =   [1      2];     % compressor connected to bus
% connect.interc.comp =  [1   5;      % compressor connected to bus
%                         2   7];

%                   idgen	node    eff
% connect.interc.term =  [1   48      10e-3]; % gas-fired unit connected to node
connect.interc.term =  []; % gas-fired unit connected to node


