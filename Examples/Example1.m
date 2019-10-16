% MPNG - Example 1. Optimal power and natural gas flow. 8 node natural gas 
%   system and a 9 bus power system, for 4 equal periods of time during a day

%% Define useful constants
define_constants                % power system constants
define_constants_gas            % natural gas system and connect struct constants
%% defining option vector
mpopt = mpoption;                   % initialize option struct
mpopt.opf.ac.solver = 'IPOPT';      % current stable sover
mpopt.ipopt.opts.max_iter = 1e5;    % max iterations
%% Load cases
mpc = loadcase('case9_new');    % 9 bus power system with a little modifications
mgc = ng_case8;                 % 8 node natural gas system
connect = connect_pg_case9;     % interconnection case, it mainly has to fit with the power case
%% make some modifications
% time
connect.power.time = [6 6 6 6];     % four periods of time, all the same length
% power system demands
factors = [1 1.1 1.2 0.9];                          % factors (external aid of MPNG)
connect.power.demands.pd = factors.*mpc.bus(:,PD);  % PD matrix
connect.power.demands.qd = factors.*mpc.bus(:,QD);  % QD matrix
% energy available in gen 2
connect.power.energy = zeros(1,2);          % initialize (one gen with max energy)
connect.power.energy(1,GEN_ID) = 2;         % id for max energy in gen 2 
connect.power.energy(1,MAX_ENER) = 5000;    % max energy is 1000 MWh/d 
% spinning reserve
connect.power.sr = [];                      % no spinning reserves
% interconnection through compressors
connect.interc.comp = zeros(1,2);           % initialize (one compressor working with power)
connect.interc.comp(1,COMP_ID) = 2;         % second compressor
connect.interc.comp(1,BUS_ID)  = 5;         % connected to bus # 5
mgc.comp(2,TYPE_C) = COMP_P;                % create concordance with gas case
% interconnection through gas-fired power plants
connect.interc.term = zeros(1,3);           % initialize (one gas-fired plant)
connect.interc.term(1,GEN_ID)   =  3;       % third power plant
connect.interc.term(1,NODE_ID)  =  7;       % connected to node # 7
connect.interc.term(1,EFF)      =  10e-3;   % power plant efficiency (MMSCFD/MWh)
%% putting all together
mpgc = mpc;                         % initialize MPNG case
mpgc.mgc = mgc;                     % adding gas case
mpgc.connect = connect;             % adding connection struct
%% running program 
res = mpng(mpgc,mpopt);             % running MPNG