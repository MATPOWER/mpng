function mgc = ng_case8
%% ========================================================================
% 8-node natural gas system for optimal power & natural gas flow.
% =========================================================================
% Authors:  (1) Wilson González Vanegas, M.Sc 
%           (2) Sergio García Marín, M.Sc(c)  
% (1) Automatic Research Group, Universidad Tecnológica de Pereira.
% (2) GIPEM, Universidad Nacional de Colombia - Sede Manizales.
% =========================================================================
%% ---------------------- Matpower gas case version ----------------------
mgc.version = '1'; 

%% ------------------------------ Gas bases ------------------------------

mgc.pbase = 500; % base pressure [psia]
mgc.fbase = 50;  % base flow [MMSCFD]
mgc.wbase = 100; % base power [MW]

%% ------------------------------ Node data ------------------------------                                    
% node_id  type  p 	Pmax    Pmin	pi_+	pi_-	al_pi_+     al_pi_-     gd  #gd              

mgc.node.info = [
    1   2	435.1130    650  406  0    0  10e3  10e3      0         2;
    2	1	406.1055    650  377  0    0  10e3  10e3      0         2;
    3	1	580.1507    653  508  0    0  10e3  10e3      0         2;
    4	1	565.6470	650  522  0    0  10e3  10e3      0         2;
    5	1	507.6319	650  464  0    0  10e3  10e3      0         2;
    6	1	551.1432	650  522  0    0  10e3  10e3      0         2;
    7	1	507.6319	650  464  0    0  10e3  10e3      12.2186   2;
    8	1	507.6319	650  464  0    0  10e3  10e3      26.4149   2;
    ];

%   Res     Ind         Com     NGV         Ref     Pet      <---- [MMSCFD]
mgc.node.dem = [
    0         0;
    0         0;
    0         0;
    0         0;
    0         0;
    0         0;
    12.2186   0;
    18.4904   7.9245;
    ];

% al_Res    al_Ind	al_Com  al_NGV	al_Ref	al_Pet          <----[$/MMSCFD]
mgc.node.demcost = 1e6*[
    0       0;
    0       0;
    0       0;
    0       0;
    0       0;
    0       0;
    50      0;
    50      25;
    ];
    
%% ------------------------------ Well data ------------------------------     
% node  I   pw  Imax   Imin     status      Cg
mgc.well = [
   1	38.6335  445.9454	80     0     1    5e3;
            ];
%% ---------------------------- Pipeline data ----------------------------
% fnode  tnode  FG_O    Kij	D   L   Fg_max Fg_min C_O 
mgc.pipe = [
    1   2	1	0.1412      0	0   80.00   -80.00      0.0050e3;
    4   5	1	0.1214      0	0   40.00   -40.00      0.0050e3;
    4	6	1	0.1567      0	0   40.00   -40.00      0.0050e3;
    5	6	1	0.0604      0	0   40.00	-40.00      0.0050e3;
    6	7	1	0.0736      0	0   40.00	-40.00      0.0050e3;
    5	8	1	0.0736      0	0   40.00   -40.00      0.0050e3;
    ];
%% ---------------------------- Compressor data ----------------------------
% fnode  tnode  Type    fgc     pcc     gcc     ratio   bc  zc  al  be  ga      fmaxc   costc
mgc.comp = [
    2      3     2      0       0       0       1.05    4.8808  0.236   0   0.00025    0   80     0.005e3;
    3      4     1      0       0       0       1.15    4.8808  0.236   0   0.00025    0   80     0.005e3;
    ];
%% -------------------------------- Storage --------------------------------
%   node   S  S0  Smax Smin  fs  fs+ fs- fsmax   fsmin   sta C_S C_S+ C_S-
mgc.sto = [
    %3	3   3   5   1	0   0   0   3   -3   1   10e3   30e3  20e3;
    %5 	0   0   0   0	0   0   0   0   0   0   0   0   0;
    %7 	0   0   0   0	0   0   0   0   0   0   0   0   0;
    ];

%% Names
mgc.nodenames = {
    'Node 1','Node 2','Node 3','Node 4','Node 5','Node 6','Node 7'
    };
mgc.wellnames = {
    'Well_1'
    };

             
                  