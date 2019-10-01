function mpc_out = multi_period(mpc_in,time,pd,qd)
% MULTI_PERIOD constructs a matpower case for multiple periods of analysis 
%   through the creation of islands for every period of time.
% 
%   mpc_out = multi_period_mpc(mpc_in,factors)
%   mpc_out = multi_period_mpc(mpc_in,time,pd,qd)
% 
%   Inputs:     mpc_in  -> Matpower original case
%               factors -> Factor which multiplies the demands in every period 
%                          of time. Lenght of "factors" determines the
%                          number of periods. 
%               times   -> A vector with the factors for every period of analysis.          
%               pd      -> A matrix with nb rows and nt columns. Every column
%                          represents the active power demand in all buses for 
%                          each period of time.
%               qd      -> A matrix with nb rows and nt columns. Every column
%                          represents the reactive power demand in all buses for 
%                          each period of time.

%   MPNG: Matpower - Natural Gas
%   Copyright (c) 2019 - v0.99alpha
%   Sergio García Marín - Universidad Nacional de Colombia - Sede Manizales
%   Wilson González Vanegas - Universidad Tecnológica de Pereira
%   Carlos E. Murillo Sánchez - Universidad Nacional de Colombia - Sede Manizales

%   3-clause bsd license 

%%
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

mpc = mpc_in;

[nb,cb]= size(mpc.bus);
[ng,cg]= size(mpc.gen);
[nl,cl] = size(mpc.branch);
[ngen,cgen] = size(mpc.gencost);
nt = length(time);

if nargin == 1
     error('multi_period: Not enough input arguments');
elseif nargin == 2
    factors = time; %number of periods
    pd = mpc.bus(:,PD);
    pd = factors.*pd;
    qd = mpc.bus(:,QD);
    qd = factors.*qd;
elseif nargin > 2
    if nb ~= size(pd,1) || nb ~= size(qd,1)
        error('multi_period: Input demands must agree with number of buses.');
    elseif  size(pd) ~= size(qd)
        error('multi_period: Active and reactive power demands disagree.');
    end
    if length(time) ~= size(pd,2)
        error('multi_period: Demands don´t match with number of periods.');
    end
end
if sum(time) ~= 24
    error('multi_period: Total number of hours don´t correspond to a day.');
end

bus = zeros(nb*nt,cb);
gen = zeros(ng*nt,cg);
branch = zeros(nl*nt,cl);
gencost = zeros(ngen*nt,cgen);
if isfield(mpc, 'busnames') && iscell(mpc.busnames)
    busnames = cell(nb*nt, 1);
end
if isfield(mpc, 'genfuel') && iscell(mpc.genfuel)
    genfuel = cell(ng*nt, 1);
end
if isfield(mpc, 'gennames') && iscell(mpc.gennames)
    gennames = cell(ng*nt, 1);
end

for i = 0:(nt-1)
    
    j = i+1;
    bus(((i*nb)+1):j*nb,:) = mpc.bus;
    bus(((i*nb)+1):j*nb,PD) = pd(:,j);
    bus(((i*nb)+1):j*nb,QD) = qd(:,j);
    
    gen(((i*ng)+1):j*ng,:) = mpc.gen;
    gb_temp = mpc.gen(:,GEN_BUS);
    gen(((i*ng)+1):j*ng,GEN_BUS) = gb_temp + nb*i;
    
    branch(((i*nl)+1):j*nl,:) = mpc.branch;
    fb_temp = mpc.branch(:,F_BUS);
    tb_temp = mpc.branch(:,T_BUS);
    branch(((i*nl)+1):j*nl,F_BUS) = fb_temp + nb*i;
    branch(((i*nl)+1):j*nl,T_BUS) = tb_temp + nb*i;
    
    gencost(((i*ngen)+1):j*ngen,:) = mpc.gencost;
    gc_temp = time(j)*gencost(((i*ngen)+1):j*ngen,COST:end);
    gencost(((i*ngen)+1):j*ngen,COST:end) = gc_temp;
    
    if isfield(mpc, 'busnames') && iscell(mpc.busnames)
        busnames((nb*i + 1):(nb*j)) = mpc.busnames;
    end
    if isfield(mpc, 'genfuel') && iscell(mpc.genfuel)
        genfuel((ng*i + 1):(ng*j)) = mpc.genfuel;
    end
    if isfield(mpc, 'gennames') && iscell(mpc.gennames)
        gennames((ng*i + 1):(ng*j)) = mpc.gennames;
    end
   
end
%% 
bus(:,BUS_I) = 1:nb*nt;

mpc_out = mpc;
mpc_out.bus = bus;
mpc_out.gen = gen;
mpc_out.branch = branch;
mpc_out.gencost = gencost;
if isfield(mpc, 'busnames') && iscell(mpc.busnames)
    mpc_out.busnames = busnames;
end
if isfield(mpc, 'genfuel') && iscell(mpc.genfuel)
    mpc_out.genfuel = genfuel;
end
if isfield(mpc, 'gennames') && iscell(mpc.gennames)
    mpc_out.gennames = gennames;
end
mpc_out.multi_period.time = time;
mpc_out.multi_period.status = 1;
% mpc_out.multi_period.original.nb = nb;
% mpc_out.multi_period.original.ng = ng;
% mpc_out.multi_period.original.nl = nl;
% mpc_out = ext2int(mpc_out);







