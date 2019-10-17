function connect = connect_pg_case9
%% ========================================================================
% Input necessary to connect power and natural gas systems
% =========================================================================
% Authors:  (1) Wilson Gonz�lez Vanegas, M.Sc 
%           (2) Sergio Garc�a Mar�n, M.Sc(c)  
% (1) Automatic Research Group, Universidad Tecnol�gica de Pereira.
% (2) GIPEM, Universidad Nacional de Colombia - Sede Manizales.
% =========================================================================
%% ---------------------- Matpower gas case version ----------------------
connect.version = '1'; 

%% ---------------------- Matpower gas case version ----------------------

connect.power.time =  [6 8 4 6];  % number of hours for each period of time sum()=24
connect.power.demands.pd = [     
         0         0         0         0;
         0         0         0         0;
         0         0         0         0;
         0         0         0         0;
   45.0000   67.5000   90.0000  112.5000;
         0         0         0         0;
   50.0000   75.0000  100.0000  125.0000;
         0         0         0         0;
   62.5000   93.7500  125.0000  156.2500;
   ];       
connect.power.demands.qd = [
         0         0         0         0;
         0         0         0         0;
         0         0         0         0;
         0         0         0         0;
   -15.0000   -22.5000   -30.0000   -37.5000;
         0         0         0         0;
   17.5000   26.2500   35.0000   43.7500;
         0         0         0         0;
   25.0000   37.5000   50.0000   62.5000;   
];

% connect.power.time =  [24];  % number of hours for each period of time sum()=24
% connect.power.demands.pd = [
%     0
%     0
%     0
%     0
%     45.0000
%     0
%     50.0000
%     0
%     62.5000
%     ];
% connect.power.demands.qd = [
%     0
%     0
%     0
%     0
%     15.0000
%     0
%     17.5000
%     0
%     25.0000
%     ];
%%
connect.power.cost = 5000; % power for non-suplied power demand

% matrix, nt -> colums, nz -> rows
connect.power.sr = [                            
    220 40  50  70;
    50 10  60  10;
%     300.0000  200.0000  150.0000  120.0000;
    ];
% connect.power.sr = [                            
%     220 40  50  70;
%     50 10  60  10;
%     50 10  60  10;
%     ];  
%   genid       maxenergy
connect.power.energy = ...
    [
    1   24*100;
    2   24*100;
    ];%


%% ---------------------- Matpower gas case version ----------------------
%                   idcomp      bus
connect.interc.comp =  [2   5];     % compressor connected to bus
connect.interc.comp =  [];     
% connect.interc.comp =  [1   5;      % compressor connected to bus 
%                         2   7]; 

%                   idgen	node    eff    
connect.interc.term =  [3   7   10e-3];% gas-fired unit connected to node
% connect.interc.term =  [];

% connect.interc.term =  [
%     1   40	10e-3;
%     3   7   8e-3;
%     ];


             
                  