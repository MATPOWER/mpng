% MPNG - Example 4. Optimal power and natural gas flow given unit commitment
%                   and generation ramp-rates. The system consists of modified
%                   versions of the 118-node power system (MATPOWER's 'case118') 
%                   and the 48-node meshed natural gas network included in the 
%                   'Cases' folder of MPNG ('ng_case_48').
%
% Simulation is customized following the case study presented in the paper: 
%
% "Unit Commitment With an Enhanced Natural Gas-Flow Model" - Sheng Chen,
% Antonio J. Conejo, Ramteen Sioshansi, and Zhinong Wei.
%
% https://doi.org/10.1109/TPWRS.2019.2908895
%
% Main author: Wilson González-Vanegas


%   MPNG Matpower - Natural Gas
%   Copyright (c) 2019-2022 - v0.99beta
%   Sergio García-Marín - Universidad Nacional de Colombia - Sede Manizales
%   Wilson González-Vanegas - Universidad Tecnológica de Pereira
%   Carlos E. Murillo-Sánchez - Universidad Nacional de Colombia - Sede Manizales
% 
%   This file is part of MPNG.
%   Covered by the 3-clause BSD License (see LICENSE file for details).


%% Define indexing constants
define_constants
define_constants_gas

%% Define solution options
mpopt = mpoption;                   % Initialize options struct
mpopt.opf.ac.solver = 'IPOPT';      % Current stable solver is IPOPT
mpopt.ipopt.opts.max_iter = 1e5;    % Set max iterations
mpopt.opf.start = 2;                % Pass the initial point in cases to IPOPT

%% Load gas & power cases and interconnection case
mpc = loadcase('case118');          % IEEE 118-bus case distributed with MATPOWER
mgc = ng_case48;                    % 48-node natural gas system distributed in MPNG
connect = connect_pg_case118;       % Interconnection case distributed in MPNG

%% Modify cases and define the interconnection set-up
% Define flow and power bases
mgc.fbase = 1e4;                         
mgc.wbase = 1e2;

% Remove bus_name field from power system case for island-creation purposes
mpc = rmfield(mpc,'bus_name');

% Branch power capacities
mpc.branch(:,RATE_A) = [175
                        175
                        500
                        175
                        175
                        175
                        500
                        500
                        500
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        500
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        500
                        500
                        500
                        175
                        175
                        500
                        175
                        500
                        175
                        175
                        140
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        500
                        500
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        500
                        175
                        175
                        500
                        500
                        500
                        500
                        500
                        500
                        500
                        175
                        175
                        500
                        175
                        500
                        175
                        175
                        500
                        500
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        500
                        175
                        175
                        175
                        175
                        175
                        175
                        500
                        500
                        175
                        500
                        500
                        200
                        200
                        175
                        175
                        175
                        500
                        500
                        175
                        175
                        500
                        500
                        500
                        175
                        500
                        500
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        200
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        500
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        175
                        500
                        175
                        175
                        175
                        500
                        175
                        175
                        175];
                    
% Minimum gen production
mpc.gen(:,PMIN) =  [10       % Change minimum gen production. Maximum gen 
                    10       % production as included in the case (no changes)
                    10
                    10
                    55
                    18.5
                    10
                    10
                    10
                    10
                    32
                    41.4
                    10
                    10.7
                    10
                    10
                    10
                    10
                    10
                    11.9
                    30.4
                    14.8
                    10
                    10
                    25.5
                    26
                    10
                    49.1
                    49.2
                    80.52
                    10
                    10
                    10
                    10
                    10
                    10
                    57.7
                    10
                    10.4
                    70.7
                    10
                    10
                    10
                    10
                    35.2
                    14
                    10
                    10
                    10
                    10
                    13.6
                    10
                    10
                    10];                    
                
% Generation ramp-rates
mpc.gen(:,RAMP_30)  =  [20                        
                        20
                        20
                        20
                        330
                        111
                        60
                        20
                        20
                        20
                        192
                        248.4
                        20
                        21.4
                        20
                        20
                        20
                        20
                        20
                        23.8
                        60.8
                        29.6
                        20
                        20
                        51
                        156
                        20
                        294.6
                        295.2
                        483.12
                        20
                        20
                        20
                        20
                        20
                        20
                        115.4
                        20
                        20.8
                        141.4
                        20
                        20
                        20
                        20
                        70.4
                        28
                        20
                        20
                        20
                        20
                        27.2
                        20
                        20
                        20];
mpc.gen(:,RAMP_30) = (1/2)*mpc.gen(:,RAMP_30); % MPNG multiplies internally by 2.0                   

% Generation costs (linear cost)
mpc.gencost(:,COST) = zeros(size(mpc.gen,1),1);
mpc.gencost(:,COST+1)  =   [41                            
                            41
                            41
                            41
                            3.5        % Ope.&Maint. of gas-fired plants
                            3.5        % Ope.&Maint. of gas-fired plants
                            3.5        % Ope.&Maint. of gas-fired plants
                            41
                            41
                            41
                            3.5       % Ope.&Maint. of gas-fired plants
                            3.5       % Ope.&Maint. of gas-fired plants
                            41
                            172.85699
                            41
                            41
                            41
                            41
                            41
                            82.631604
                            34.9019584
                            50.833284
                            41
                            41
                            36.4516055
                            3.5        % Ope.&Maint. of gas-fired plants
                            41
                            3.5        % Ope.&Maint. of gas-fired plants
                            3.5        % Ope.&Maint. of gas-fired plants
                            3.5        % Ope.&Maint. of gas-fired plants
                            41
                            41
                            41
                            41
                            41
                            41
                            32.0964588
                            41
                            280
                            31.6474715
                            41
                            41
                            41
                            41
                            33.96824
                            55
                            41
                            41
                            41
                            41
                            57.777808
                            41
                            41
                            41];                      

% Start-up/Shut-down generation costs
mpc.gencost(:,STARTUP) = (500)*ones(size(mpc.gen,1),1);
mpc.gencost(:,SHUTDOWN) = (20)*ones(size(mpc.gen,1),1);

% Periods
connect.power.time = ones(1,24);    % 24-hour horizon

% Power system demands
power_factors   =  [1.2142440167    % P/Q demand factors
                    1.0624635138
                    0.986573263
                    0.9638061875
                    0.9486281376
                    1.0017513129
                    1.1231757147
                    1.1838879156
                    1.2977232931
                    1.4115586693
                    1.4646818446
                    1.4798598945
                    1.4343257448
                    1.3280793929
                    1.267367192
                    1.2825452419
                    1.312901343
                    1.4874489201
                    1.7302977237
                    1.7227086981
                    1.7075306482
                    1.4358435494
                    1.2628137769
                    1.0897840044
                    ];                         
connect.power.demands.pd = power_factors'.*mpc.bus(:,PD);  % PD matrix
connect.power.demands.qd = power_factors'.*mpc.bus(:,QD);  % QD matrix

% Non-supplied power demand cost
connect.power.cost = 300;     % Power for non-suplied power demand   

% Daily energy usage for reservoir-based hydroelectric plants
connect.power.energy = [];    % Do not consider maximum hydro-energy constraints


% Pressure limits and initial point for pressures
mgc.node.info(:,PRMAX) = (1150)*ones(size(mgc.node.info,1),1);  % Maximum nodal pressure is 1150 psia
mgc.node.info(:,PRMIN) = (300)*ones(size(mgc.node.info,1),1);   % Minimum nodal pressure is 300 psia except
mgc.node.info([1 2],PRMIN) = [900 990];                         % in nodes 1 and 2 (900 psia and 990 psia, respectively) 
mgc.node.info(:,PR)  = [1150.00                                 % Change initial point for pressures                        
                        758.11 
                        853.03 
                        912.22 
                        914.23 
                        851.99 
                        823.28 
                        950.00 
                        1013.30 
                        963.78 
                        911.99 
                        887.21 
                        1034.10 
                        957.15 
                        876.82 
                        944.92 
                        945.70 
                        856.56 
                        873.45 
                        737.35 
                        796.32 
                        941.54 
                        852.36 
                        787.63 
                        904.75 
                        844.41 
                        843.79 
                        794.32 
                        730.84 
                        717.02 
                        716.78 
                        722.58 
                        728.34 
                        712.76 
                        704.80 
                        702.28 
                        734.97 
                        670.39 
                        652.58 
                        652.55 
                        659.10 
                        674.21 
                        702.22 
                        749.12 
                        781.83 
                        819.21 
                        780.46 
                        833.19];

% Compute initial pipeline flows (non-linear approx. of Weymouth's Eq.)
mgc.pipe(:,K_O) = 24*mgc.pipe(:,K_O);                                         % The factor 24 is to adjust the single-day window for the gas system
K = mgc.pipe(:,K_O)/(mgc.fbase/(mgc.pbase)); K = K.^2;                        % Extract Weymouth constants
pi_star = ((0.2*mgc.pipe(:,FMAX_O)/mgc.fbase).^2)./K;                         % Compute Pi* using the default 20% of maximum pipelines flow
mgc.pipe(:,FG_O) = wey_approx(K,mgc.node.info(mgc.pipe(:,F_NODE),PR),...      % Calculate pipelines flow using non-linear approx.
                                mgc.node.info(mgc.pipe(:,T_NODE),PR),...
                                pi_star);
mgc.pipe(:,FG_O) = mgc.pipe(:,FG_O)*mgc.fbase;                                % Back to real values from P.U values

% Pipeline transportation costs
mgc.pipe(:,COST_O) = zeros(size(mgc.pipe,1),1);  % No pipeline costs

% Maximum pipeline flows
mgc.pipe(:,FMAX_O) = 24*mgc.pipe(:,FMAX_O);     % The factor 24 is to adjust the single-day window for the gas system
mgc.pipe(:,FMIN_O) = 24*mgc.pipe(:,FMIN_O);     
    
% Natural gas demand
mgc.node.dem = [0                       % Mean gas demand
                0
                0
                0
                0
                0
                0
                0
                400
                0
                100
                0
                0
                0
                0
                50
                0
                0
                0
                0
                0
                0
                100
                0
                400
                0
                50
                0
                0
                30
                30
                0
                30
                30
                30
                30
                0
                30
                30
                30
                30
                40
                40
                40
                100
                200
                180
                0];
mgc.node.dem = 24*mgc.node.dem;                     % The factor 24 is to adjust the single-day window for the gas system
mgc.node.info(:,GD) = mgc.node.dem;                 % Use a single user (no-diversified demand)
mgc.node.info(:,NGD) = ones(numel(mgc.node.dem),1); 

gas_factor =  [0.6438232416                         % gas factors             
                0.6050579688
                0.5762837456
                0.6122515248
                0.621842932
                0.645821452
                0.7653144056
                0.9079865952
                1.1150011448
                1.1409778736
                1.085027996
                1.0027017464
                0.9103844472
                0.8844077184
                0.8688216808
                0.8600295568
                0.8692213224
                0.81247216
                0.7896925672
                0.7417355288
                0.7173573672
                0.7173573672
                0.7009720456
                0.6054576112];
gas_factor = mean(gas_factor);                 % MPNG works with 1-day horizon on the gas side. Thus, use the mean of the factors
mgc.node.dem = gas_factor*mgc.node.dem;        % Change demand accordingly for a single-user on the 48-node gas network
mgc.node.info(:,GD) = sum(mgc.node.dem,2);     % Change total nodal gas demand

% Non-supplied demand costs   
mgc.node.demcost = (6000)*ones(size(mgc.node.info,1),1);

% Well's limits and gas production costs
mgc.well(:,GMAX) = [850                    
                    600
                    280
                    280
                    280
                    280
                    150
                    150
                    550];
mgc.well(:,GMIN) = [600
                    400
                    0
                    0
                    0
                    0
                    0
                    0
                    0];
mgc.well(:,COST_G)  =  [1000               
                        1000
                        850
                        850
                        850
                        850
                        800
                        800
                        900];
mgc.well(:,GMAX) = 24*mgc.well(:,GMAX);                    % The factor 24 is to adjust the single-day window for the gas system
mgc.well(:,GMIN) = 24*mgc.well(:,GMIN);     

% Parameters for compressors
mgc.comp(:,B_C) = (1/24)*mgc.comp(:,B_C);                  % The factor 24 is to adjust the single-day window for the gas system
comp_eff = 0.03;                                           % The model from Shen Cheng et al. proposes: GasCons = comp_eff*GasTrans.
mgc.comp(:,BETA) = comp_eff./(mgc.comp(:,B_C)*0.5);        % Try to capture the linear relationship for a expected ratio of 1.5
mgc.comp(:,RATIO_C) = 1.5*ones(size(mgc.comp,1),1);        % Maximum compression ratio: 1.5 

% Compressors flow capacity and compression costs
mgc.comp(:,FMAX_C) = (2000)*ones(size(mgc.comp,1),1);      % 2000 MMSCFD maximum flow for all compressors
mgc.comp(:,FMAX_C) = 24*mgc.comp(:,FMAX_C);                % The factor 24 is to adjust the single-day window for the gas system
mgc.comp(:,COST_C) = zeros(size(mgc.comp,1),1);            % No compression costs

% spinning reserve
connect.power.sr = [];                      % no spinning reserves

% interconnection through compressors
connect.interc.comp = [];                   % no power-driven compressors

% interconnection through gas-fired power plants
connect.interc.term = zeros(9,3);           % initialize memory for 9 gas-fired plants
connect.interc.term(:,GEN_ID)   =  [5       % gen ids of gas-fired plants
                                    6
                                    7
                                    11
                                    12
                                    26
                                    28
                                    29
                                    30];
connect.interc.term(:,NODE_ID)  =  [45       % node ids of gas-fired plants
                                    33
                                    36
                                    31
                                    30
                                    40
                                    39
                                    38
                                    25];
connect.interc.term(:,EFF)      =  24*0.16*ones(9,1);   % power plant efficiencies
                                                        % The factor 24 is to adjust the single-day window for the gas system

%% Initial point for generation in gas-fired plants
mpc.gen(connect.interc.term(:,GEN_ID),PG) = mpc.gen(connect.interc.term(:,GEN_ID),PMIN);

% All generators ON for base case
connect.power.UC = ones(size(mpc.gen,1),24);
connect.power.ramp_time = ones(size(connect.power.UC));

%% Create MPNG-case with all structs -> Base case
mpgc = mpc;                         % initialize struct with the power system case
mpgc.mgc = mgc;                     % add gas case
mpgc.connect = connect;             % add connection struct
% run OPNGF for Base Case
fprintf('\n ======================================================')
fprintf('\n Running Base Case simulation with all generators ON... \n ')
fprintf('======================================================\n\n\n\n')

res_base_case = mpng(mpgc,mpopt);             % call MPNG

%% Now shut down a subset of gens to get a more economic operation
connect.power.UC = [
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1	1	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
0	0	0	0	0	0	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1	1	0	0	0
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	0	0	0	0	0	0	0	0	0	0	0	1	1	1	1	0	0	0
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1
];

%% Update MPNG-case with connect struct -> Economic case
mpgc.connect = connect;             % add connection struct
%% run OPNGF for Base Case
fprintf('\n\n\n ===================================')
fprintf('\n Running Economic Case simulation... \n ')
fprintf('===================================\n\n\n\n')

res_economic = mpng(mpgc,mpopt);    % call MPNG

%% Extract results and make some plots for comparing the results
UC_base_case = uc2image(res_base_case);
UC_economic = uc2image(res_economic);

f_base_case = res_base_case.f + susd_cost(res_base_case);
f_economic = res_economic.f + susd_cost(res_economic);

[pg_gas_base_case,pg_nongas_base_case,p_loss_base_case,p_supplied_base_case] = extractPs(res_base_case);
[pg_gas_eco_case,pg_nongas_eco_case,p_loss_eco_case,p_supplied_eco_case] = extractPs(res_economic);

p_nonserved_base_case = res_base_case.nsd.original.PD;  
p_nonserved_base_case = p_nonserved_base_case(p_nonserved_base_case(:,1) ~= 0);
p_nonserved_base_case = sum(p_nonserved_base_case - p_supplied_base_case);

p_nonserved_eco_case = res_economic.nsd.original.PD;  
p_nonserved_eco_case = p_nonserved_eco_case(p_nonserved_eco_case(:,1) ~= 0);
p_nonserved_eco_case = sum(p_nonserved_eco_case - p_supplied_eco_case);

p_served_base_case = sum(res_base_case.nsd.original.PD)-p_nonserved_base_case;
p_served_eco_case = sum(res_economic.nsd.original.PD)-p_nonserved_eco_case;

powers = [pg_nongas_base_case       pg_nongas_eco_case
          pg_gas_base_case          pg_gas_eco_case
          -p_loss_base_case         -p_loss_eco_case              
          -p_served_base_case       -p_served_eco_case];
                
figure('units','normalized','outerposition',[0 0 1 1]);
ax1 = subplot(2,6,[1,2]);
imagesc(UC_base_case,'AlphaData',0.8)
colormap(ax1,[0 1 1; 1 1 0]);
clabels = {'on non-gas','on gas'};
colorbar(ax1,'location', 'EastOutside','Ticks',[2 4],'TickLabels',clabels);
ax1.XLabel.String = 'Period (t)';
ax1.XAxisLocation = 'top';
ax1.FontSize = 11;
ax1.XTick = [1:24];
ax1.YTick = [1:3:54];
ax1.YLabel.String = 'Gen. Num. Base Case';
ax1.YAxisLocation = 'left';
ax1.YDir = 'reverse';

ax2 = subplot(2,6,[4,6]);
imagesc(UC_economic,'AlphaData',0.8)
colormap(ax2,[0 0 1; 0 1 1; 1 0 0; 1 1 0]);
clabels = {'off non-gas','on non-gas','off gas','on gas'};
colorbar(ax2,'location', 'EastOutside','Ticks',[1:4],'TickLabels',clabels);
ax2.XLabel.String = 'Period (t)';
ax2.XAxisLocation = 'top';
ax2.FontSize = 11;
ax2.XTick = [1:24];
ax2.YTick = [1:3:54];
ax2.YLabel.String = 'Gen. Num. Economic Case';
ax2.YAxisLocation = 'left';
ax2.YDir = 'reverse';

ax3 = subplot(2,6,[7,7]);
B1 = bar(categorical({'Base Case','Economic Case'}),[f_base_case f_economic]);
B1.FaceColor = 'flat';
B1.CData(1,:) = [1 0.1 0];
B1.CData(2,:) = [0 1 0.1];
ylabel(['Obj. Function (','$)'])
ax3.FontSize = 11;

ax4 = subplot(2,6,[9,9]);
B2 = bar(categorical({'Base Case','Economic Case'}),[p_nonserved_base_case p_nonserved_eco_case]);
B2.FaceColor = 'flat';
B2.CData(1,:) = [1 0.1 0];
B2.CData(2,:) = [0 1 0.1];
ylabel('Non-served power (MW)')
ax4.Title.String = {['\alpha_{\epsilon} = ',num2str(300),'$/MW'];... % Non-served active power cost
                    ['Max Gen.Cost = ',num2str(172.86),'$/MW'];...   % Maximum provided cost of generation
                    [' ']};        
ax4.FontSize = 11;

% ax5 = subplot(2,6,[10,12]);
% bar(categorical({'Base Case','Economic Case'}),powers,'stacked')
% legend({'Non Gas-fired Gen.','Gas-fired Gen.','Losses','Served Power'},...
%         'location','westoutside')
% ylabel('Active power (MW)')
% ax5.FontSize = 11;



%% Appendix A) ------ An auxiliary function to compute start-up/shut-down costs ------
function cost = susd_cost(res)
define_constants
nt = size(res.connect.power.time,2);
ng = size(res.genid.original,1)/nt;
if isempty(res.connect.power.UC)
    res.connect.power.UC = ones(ng,nt);
end
start_up_cost = res.gencost(1:ng,STARTUP);
shut_down_cost = res.gencost(1:ng,SHUTDOWN);
cost = 0;
for t = 1:nt-1
    T  = res.connect.power.UC(:,t);
    T1 = res.connect.power.UC(:,t+1);
    difference = T - T1;
    on_cost = sum(start_up_cost(difference<0));
    off_cost = sum(shut_down_cost(difference>0));
    cost = cost + on_cost + off_cost;
end
end

%% Appendix B) ------ An auxiliary function for mapping UC to colors by gen-type -------
function UCcolor = uc2image(res)
define_constants
define_constants_gas
id_gas_gen = res.connect.interc.term(:,GEN_ID);
nt = size(res.connect.power.time,2);
ng = size(res.genid.original,1)/nt;
if isempty(res.connect.power.UC)
    res.connect.power.UC = ones(ng,nt);
end
mask_gas_gen = zeros(ng,1);
mask_gas_gen(id_gas_gen) = 1;
mask_gas_gen = repmat(mask_gas_gen,1,nt);
UCcolor = ones(size(mask_gas_gen));
UCcolor(res.connect.power.UC == 1 & mask_gas_gen == 1) = 4;
UCcolor(res.connect.power.UC == 0 & mask_gas_gen == 1) = 3;
UCcolor(res.connect.power.UC == 1 & mask_gas_gen == 0) = 2;
UCcolor(res.connect.power.UC == 0 & mask_gas_gen == 0) = 1;
end

%% Appendix C) ------ An auxiliary function for extracting different active-power quantities ------
function [pg_gas,pg_nongas,p_loss,p_supplied] = extractPs(res)
define_constants
define_constants_gas
nt = size(res.connect.power.time,2);
ng = size(res.genid.original,1)/nt;
if isempty(res.connect.power.UC)
    res.connect.power.UC = ones(ng,nt);
end
id_gas_gen = res.connect.interc.term(:,GEN_ID);
id_gas_gen = repmat(id_gas_gen,nt,1) + reshape(ng*repmat([0:nt-1],numel(id_gas_gen),1),[],1);
id_nongas_gen = [1:ng]';
id_nongas_gen = id_nongas_gen(~ismember(id_nongas_gen,res.connect.interc.term(:,GEN_ID)));
id_nongas_gen = repmat(id_nongas_gen,nt,1) + reshape(ng*repmat([0:nt-1],numel(id_nongas_gen),1),[],1);
pg_gas = sum(res.gen(id_gas_gen,PG));
pg_nongas = sum(res.gen(id_nongas_gen,PG));
p_loss = sum(real(get_losses(res)));
p_supplied = -res.gen(end-res.nsd.N+1:end,PG);
end