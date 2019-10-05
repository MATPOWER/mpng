function mgc = NGcase48
%% ========================================================================
% 48-node natural gas system for optimal power & natural gas flow.
% =========================================================================
% Authors:  (1) Wilson González Vanegas, M.Sc(c) 
%           (2) Sergio García Marín, M.Sc(c)  
% (1) Automatic Research Group, Universidad Tecnológica de Pereira.
% (2) GIPEM, Universidad Nacional de Colombia - Sede Manizales.
% =========================================================================
%% ---------------------- Matpower gas case version ----------------------
mgc.version = '1'; 

%% ------------------------------ Gas bases ------------------------------

mgc.pbase = 500; % base pressure [psia]
mgc.fbase = 10;  % base flow [MMSCFD]: "Million Standard Cubic Feet per Day" 
mgc.wbase = 100; % base power [KW]

%% ------------------------------ Node data ------------------------------

% Notes:
%   NODE_I      = 1;    %% node number (1 to 29997)
%   NODE_TYPE	= 2;    %% node type (1: demand node, 2: extraction node)
%   PR          = 3;    %% pressure [psia]
%   PRMAX       = 4;    %% Maximun pressure [psia]
%   PRMIN       = 5;    %% Minimun pressure [psia]
%   OVP         = 6;    %% overpressure [psia]
%   UNP         = 7;    %% underpressure [psia]
%   COST_OVP	= 8;    %% overpressure cost [$/psia]
%   COST_UNP	= 9;    %% underpressure cost [$/psia]
%   GD          = 10;	%% Total gas demand [MMSCFD]
%   NGD         = 11;	%% Number of different type of demands (positive integer)
% 
%                                         
% node_id  type  p 	Pmax    Pmin	pi_+	pi_-	al_pi_+     al_pi_-     gd  #gd              

mgc.node.info = [
    1	2	1000        1450	900	0	0	100	100	0	2	;
    2	2	950         1400	900	0	0	100	100	0	2	;
    3	1	763.527 	1450	300	0	0	100	100	0	2	;
    4	2	935.161     1400	300	0	0	100	100	0	2	;
    5	2	933.285     1400	300	0	0	100	100	0	2	;
    6	2	897.055     1400	300	0	0	100	100	0	2	;
    7	2	878.369     1400	300	0	0	100	100	0	2	;
    8	1	833.448     1400	300	0	0	100	100	0	2	;
    9	1	931.494     1400	300	0	0	100	100	300	2	;
    10	1	1007.671	1400	300	0	0	100	100	0	2	;
    11	1	909.309     1400	300	0	0	100	100	100	2	;
    12	1	839.048     1400	300	0	0	100	100	0	2	;
    13	1	929.032     1400	300	0	0	100	100	0	2	;
    14	1	865.125     1400	300	0	0	100	100	0	2	;
    15	2	800.913     1400	300	0	0	100	100	0	2	;
    16	1	851.332     1400	300	0	0	100	100	50	2	;
    17	1	852.458     1400	300	0	0	100	100	0	2	;
    18	2	781.837     1400	300	0	0	100	100	0	2	;
    19	1	796.105     1400	300	0	0	100	100	0	2	;
    20	2	672.384     1400	300	0	0	100	100	0	2	;
    21	1	726.750     1400	300	0	0	100	100	0	2	;
    22	1	785.512     1400	300	0	0	100	100	0	2	;
    23	1	645.666     1400	300	0	0	100	100	100	2	;
    24	1	517.395     1400	300	0	0	100	100	0	2	;
    25	1	753.619     1400	300	0	0	100	100	400	2	;
    26	1	631.640     1400	300	0	0	100	100	0	2	;
    27	1	628.360     1400	300	0	0	100	100	80	2	;
    28	1	531.710     1400	300	0	0	100	100	0	2	;
    29	1	439.519     1400	300	0	0	100	100	0	2	;
    30	1	416.563     1400	300	0	0	100	100	70	2	;
    31	1	412.988     1400	300	0	0	100	100	70	2	;
    32	1	414.406     1400	300	0	0	100	100	0	2	;
    33	1	415.820     1400	300	0	0	100	100	60	2	;
    34	1	411.250     1400	300	0	0	100	100	50	2	;
    35	1	401.610     1400	300	0	0	100	100	60	2	;
    36	1	401.374     1400	300	0	0	100	100	60	2	;
    37	1	492.950     1400	300	0	0	100	100	0	2	;
    38	1	450.870     1400	300	0	0	100	100	40	2	;
    39	1	425.873     1400	300	0	0	100	100	30	2	;
    40	1	411.528     1400	300	0	0	100	100	30	2	;
    41	1	405.328     1400	300	0	0	100	100	40	2	;
    42	1	404.838     1400	300	0	0	100	100	40	2	;
    43	1	405.983     1400	300	0	0	100	100	40	2	;
    44	1	429.873     1400	300	0	0	100	100	40	2	;
    45	1	474.613     1400	300	0	0	100	100	100	2	;
    46	1	576.135     1400	300	0	0	100	100	100	2	;
    47	1	463.664     1400	300	0	0	100	100	200	2	;
    48	1	711.844     1400	300	0	0	100	100	0	2	;    
    ];

%   Res     Ind         Com     NGV         Ref     Pet      <---- [MMSCFD]
mgc.node.dem = [
    0	0	;
    0	0	;
    0	0	;
    0	0	;
    0	0	;
    0	0	;
    0	0	;
    0	0	;
    180	120	;
    0	0	;
    70	30	;
    0	0	;
    0	0	;
    0	0	;
    0	0	;
    35	15	;
    0	0	;
    0	0	;
    0	0	;
    0	0	;
    0	0	;
    0	0	;
    70	30	;
    0	0	;
    280	120	;
    0	0	;
    56	24	;
    0	0	;
    0	0	;
    49	21	;
    49	21	;
    0	0	;
    42	18	;
    35	15	;
    42	18	;
    42	18	;
    0	0	;
    28	12	;
    21	9	;
    21	9	;
    28	12	;
    28	12	;
    28	12	;
    28	12	;
    70	30	;
    70	30	;
    140	60	;
    0	0	;    
    ];

% al_Res    al_Ind	al_Com  al_NGV	al_Ref	al_Pet          <----[$/MMSCFD]
mgc.node.demcost = [
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    200	500	;
    ];
    
%% ------------------------------ Well data ------------------------------     

% Notes
%     WELL_NODE   = 1;    %% Well node
%     G           = 2;    %% Gas injection [MMSCFD]
%     PW          = 3;    %% Known pressure at well [psia]
%     GMAX        = 4;    %% Gmax, maximum gas output [MMSCFD]
%     GMIN        = 5;    %% Gmin, minimum gas output [MMSCFD]
%     WELL_STATUS = 6;    %% well status 
%     COST_G      = 7;    %% gas production cost [$/MMSCFD]

% node  I   pw  Imax   Imin     status      Cg
mgc.well = [
    1	200	950	950	0	1	10	;
    2	200	950	750	0	1	10	;
    4	200	950	350	0	1	10	;
    5	200	950	300	0	1	10	;
    6	200	950	300	0	1	10	;
    7	200	950	300	0	1	10	;
    15	100	950	150	0	1	10	;
    18	100	950	150	0	1	10	;
    20	450	950	600	0	1	10	;    
    ];

%% ---------------------------- Pipeline data ----------------------------
% Notes:
% 
%     F_NODE      = 1;    %% f, from node number
%     T_NODE      = 2;    %% t, to node number
%     FG_O        = 3;    %% gas flow in pipeline
%     K_O         = 4;    %% kewmouth constant
%     DIAM        = 5;    %% diameter [inch]
%     LNG         = 6;    %% length [km]
%     FMAX_O      = 7;    %% maximun gas flow
%     FMIN_O      = 8;    %% maximun gas flow
%     COST_O      = 9;    %% cost of gas transport
%
% fnode  tnode  FG_O	Kij        D   L   Fg_max Fg_min C_O
mgc.pipe = [
    1	8	300	1.3023	0	0	950     -950	10	;
    2	4	300	1.9474	0	0	1172	-1172	10	;
    4	7	300	1.8238	0	0	1098	-1098	10	;
    5	6	300	0.6032	0	0	363     -363	10	;
    6	7	300	1.8238	0	0	1098	-1098	10	;
    7	3	300	4.1012	0	0	2469	-2469	10	;
    9	11	300	1.3023	0	0	950     -950	10	;
    10	11	300	4.1012	0	0	2469	-2469	10	;
    11	12	300	8.3042	0	0	4999	-4999	10	;
    13	14	300	1.3023	0	0	950     -950	10	;
    14	19	300	1.3023	0	0	950     -950	10	;
    15	19	300	1.3023	0	0	950     -950	10	;
    19	20	300	1.3023	0	0	950     -950	10	;
    13	17	300	2.8505	0	0	1716	-1716	10	;
    17	16	300	1.3023	0	0	784     -784	10	;
    17	18	300	2.8505	0	0	1716	-1716	10	;
    18	20	300	2.8505	0	0	1716	-1716	10	;
    25	26	300	1.3023	0	0	950     -950	10	;
    26	27	300	1.5487	0	0	932     -932	10	;
    26	28	300	1.3023	0	0	784     -784	10	;
    28	29	300	0.6092	0	0	367     -367	10	;
    29	30	300	0.6441	0	0	388     -388	10	;
    30	31	300	0.6092	0	0	367     -367	10	;
    31	32	300	0.6441	0	0	388     -388	10	;
    32	33	300	0.6441	0	0	388     -388	10	;
    33	44	300	0.6441	0	0	388     -388	10	;
    29	34	300	0.6092	0	0	367     -367	10	;
    34	35	300	0.6441	0	0	388     -388	10	;
    35	36	300	0.6441	0	0	388     -388	10	;
    36	43	300	0.6441	0	0	388     -388	10	;
    28	37	300	0.6092	0	0	367     -367	10	;
    37	38	300	0.6092	0	0	367     -367	10	;
    38	39	300	0.6092	0	0	367     -367	10	;
    39	40	300	0.6092	0	0	367     -367	10	;
    40	41	300	0.6092	0	0	367     -367	10	;
    41	42	300	0.6092	0	0	367     -367	10	;
    43	42	300	0.6441	0	0	388     -388	10	;
    44	43	300	0.6441	0	0	388     -388	10	;
    45	44	300	1.4341	0	0	863     -863	10	;
    45	47	300	3.8938	0	0	2344	-2344	10	;
    46	45	300	2.7427	0	0	1651	-1651	10	;
    22	23	300	2.7427	0	0	1651	-1651	10	;
    23	24	300	2.7533	0	0	1657	-1657	10	;
    ];
%% ---------------------------- Compressor data ----------------------------
% 
% 
% Notes:
%   F_NODE      = 1;    %% f, from node number
%   T_NODE      = 2;    %% t, to node number
%   TYPE_C      = 3;    %% type of compressor
%   FG_C        = 4;    %% gas flow through compressor  [MMSCFD]\ 
%   PC_C        = 5;	%% power consumed by compressor [KW]     | -> VARIABLES    
%   GC_C        = 6;	%% gas consumed by compressor   [MMSCFD]/
%   RATIO_C     = 7;    %% compressor ratio
%   B_C         = 8;    %% B_C [KW/MMSCFD]
%   Z_C         = 9;    %% Z_C
%   ALPHA       = 10;	%% alpha [MMSCFD]      \    
%   BETA        = 11;	%% beta  [MMSCFD/KW]    | --> compresors gas compsumption parameters
%   GAMMA       = 12;	%% gamma [MMSCFD/KW^2] /
%   FMAX_C      = 13;   %% maximum flow through compressors [MMSCFD]
%   COST_C      = 14;   %% compressor cost [$/KW]
% 
% compressor type
% 
%   COMP_P      = 1;    %% compressor working with power
%   COMP_G      = 2;    %% compressor working with gas


% fnode  tnode  Type    fgc     pcc     gcc     ratio   bc  zc  al  be  ga      fmaxc   costc
mgc.comp = [
    8	9	2	0	0	0	1.2     4.8808	0.236	0	0.00025	0	5000	10	;
    3	10	2	0	0	0	1.3     4.8808	0.236	0	0.00025	0	5000	10	;
    12	13	2	0	0	0	1.2     4.8808	0.236	0	0.00025	0	5000	10	;
    20	21	2	0	0	0	1.1     4.8808	0.236	0	0.00025	0	5000	10	;
    20	48	2	0	0	0	1.15	4.8808	0.236	0	0.00025	0	5000	10	;
    21	22	2	0	0	0	1.2     4.8808	0.236	0	0.00025	0	5000	10	;
    48	25	2	0	0	0	1.1     4.8808	0.236	0	0.00025	0	5000	10	;
    24	46	2	0	0	0	1.05	4.8808	0.236	0	0.00025	0	5000	10	;
    ];
%

%% -------------------------------- Storage --------------------------------

% Notes:
%   STO_NODE	= 1;    %% storage node 
%   STO         = 2;    %% storage level -> (after running program) 
%   STO_0       = 3;    %% storage initial level 
%   STOMAX      = 4;    %% maximun storage
%   STOMIN      = 5;    %% minimun storage
%   FSTO        = 6;    %% storage outflow diference
%   FSTO_OUT	= 7;    %% storage outflow 
%   FSTO_IN     = 8;    %% storage inflow 
%   FSTOMAX     = 9;    %% maximun storage outflow diference     
%   FSTOMIN     = 10;	%% minimun storage outflow diference
%   S_STATUS	= 11;	%% storage status 
%   COST_STO	= 12;	%% storage cost 
%   COST_OUT	= 13;   %% storage outflow cost
%   COST_IN     = 14;   %% storage inflow cost

%node   V	V0 Vmax V_max   fs  fs+ fs- fsmax   fsmin   sta	C_V C_S+ C_S-
mgc.sto = [
    3	0   0   0   0	0   0   0   0   0   0   0   0   0;
    5 	0   0   0   0	0   0   0   0   0   0   0   0   0;
    7 	0   0   0   0	0   0   0   0   0   0   0   0   0;
    ];
%% 

%% Names

% mgc.nodenames = {
%     'Node 1','Node 2','Node 3','Node 4','Node 5','Node 6','Node 7'
%     };
% 
% mgc.wellnames = {
%     'Well_1'
%     };

             
                  