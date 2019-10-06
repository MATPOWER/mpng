function mgc = NGcase8
%% ========================================================================
% 8-node natural gas system for optimal power & natural gas flow.
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
mgc.wbase = 100; % base power [MW]

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
    1   2	435.1130    464  406  0    0   25  25     0       6;
    2	1	406.1055    435  377  0    0   25  25     0       6;
    3	1	580.1507    653  508  0    0   25  25     0       6;
    4	1	565.6470	624  522  0    0   25  25     0       6;
    5	1	507.6319	551  464  0    0   25  100     26.4149 6;
    6	1	551.1432	580  522  0    0   25  25     0       6;
    7	1	507.6319	551  464  0    0   25  25     12.2186 6;
    8	1	580.1507	653  508  0    0   25  25     0       6;
    ];

%   Res     Ind         Com     NGV         Ref     Pet      <---- [MMSCFD]
mgc.node.dem = [
    0         0         0         0         0         0;
    0         0         0         0         0         0;
    0         0         0         0         0         0;
    0         0         0         0         0         0;
    18.4904   7.9245	0         0         0         0;
    0         0         0         0         0         0;
    12.2186   0         0         0         0         0;
    0         0         0         0         0         0;
    ];

% al_Res    al_Ind	al_Com  al_NGV	al_Ref	al_Pet          <----[$/MMSCFD]
mgc.node.demcost = [
    0     0     0     0     0     0;
    0     0     0     0     0     0;
    0     0     0     0     0     0;
    0     0     0     0     0     0;
    5000   5000   0     0     0     0;
    0     0     0     0     0     0;
    5000   0     0     0     0     0;
    0     0     0     0     0     0;
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
   1	38.6335  445.9454	40     0     1    5;
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

% fnode  tnode  FG_O    Kij	D   L   Fg_max Fg_min C_O 
mgc.pipe = [
    1   2	1	0.1412      0	0   56.0000         0    0.0100;
    3   4	1	0.1214      0	0   35.0000  -35.0000    0.01000;
    3	6	1	0.1567      0	0   28.0000  -28.0000    0.01000;
    4	5	1	0.0604      0	0   42.0000         0    0.01000;
    4	6	1	0.0316      0	0   17.0000  -17.0000    0.01000;
    6	7	1	0.0736      0	0   38.0000         0    0.01000;
    ];
%% ---------------------------- Compressor data ----------------------------


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
    2      8     2      0       0       0       1.2     4.8808  0.236   0   0.011839   0   40     0;
    8      3     2      0       0       0       1.2     4.8808  0.236   0   0.011839   0   40     0;
    ];
% 4.8808

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

mgc.nodenames = {
    'Node 1','Node 2','Node 3','Node 4','Node 5','Node 6','Node 7'
    };

mgc.wellnames = {
    'Well_1'
    };

             
                  