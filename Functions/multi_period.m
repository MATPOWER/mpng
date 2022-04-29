function mpc_out = multi_period(mpc0,time,pd,qd)
% MULTI_PERIOD constructs a matpower case for multiple periods of analysis 
%   through the creation of islands for every period of time inside the 
%   same case.
% 
%   MPC_OUT = MULTI_PERIOD(MPC0,TIME,PD,QD)
%   MPC_OUT = MULTI_PERIOD(MPC0,FACTORS)
% 
%   This function takes a MATPOWER case file or struct, a vector time which
%   defines the number of periods and their length, and the active and
%   reactive demand, to constructs multiple islands inside the case. It is
%   also possible to not define explicitly the demands, but defining a
%   FACTORS vector which multiplies the original demands and divides the
%   day in same intervals of time. Inputs are as follows:   
% 
%   MPC0 - File name or struct with original MATPOWER case.
%   
%   TIMES - Vector which specifies the number of hours for each period of
%       time taking into considaration. Sum of elements in the vector must 
%       equal 24.
%   
%   PD - Active power demand matrix for all buses and periods of time.
% 
%   QD - Rective power demand matrix for all buses and periods of time.
% 
%   FACTORS - Vector factor which multiplies the demands in every periods 
%       of time. Lenght of FACTORS determines the number of periods in the 
%       day, all with the same length. 
%
%   See also MPC2GAS_PREP


%   MPNG: MATPOWER - Natural Gas
%   Copyright (c) 2019-2022 - v0.99beta
%   Sergio García-Marín - Universidad Nacional de Colombia - Sede Manizales
%   Wilson González-Vanegas - Universidad Tecnológica de Pereira
%   Carlos E. Murillo-Sánchez - Universidad Nacional de Colombia - Sede Manizales
% 
%   This file is part of MPNG.
%   Covered by the 3-clause BSD License (see LICENSE file for details).

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

mpc = mpc0;

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
    nt = length(factors);
    time = ones(1,nt)*(24/nt);
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
    gc_temp = modcost(gencost(((i*ngen)+1):j*ngen,:),time(j));
    gencost(((i*ngen)+1):j*ngen,:) = gc_temp;
    
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
%% Create new ids for generators
if isfield(mpc, 'genid')
    idgen_original_old = mpc_out.genid.original;
    n_original_old = length(idgen_original_old);
    idgen_original = [];
    for i = 0:(nt-1)
        idgen_original((n_original_old*i + 1):(n_original_old*i + n_original_old),1) = (idgen_original_old) + (ng*i);
    end    
    mpc_out.genid.original = idgen_original;
    
    if isfield(mpc_out.genid,'comp')
        idgen_comp_old = mpc_out.genid.comp;
        n_comp_old = length(idgen_comp_old);
        idgen_comp = [];
        for i = 0:(nt-1)
            idgen_comp((n_comp_old*i + 1):(n_comp_old*i + n_comp_old),1) = (idgen_comp_old) + (ng*i);
        end
        mpc_out.genid.comp = idgen_comp;
    end
end
%% 
mpc_out.multi_period.time = time;
mpc_out.multi_period.status = 1;
end 


% mpc_out.multi_period.original.nb = nb;
% mpc_out.multi_period.original.ng = ng;
% mpc_out.multi_period.original.nl = nl;
% mpc_out = ext2int(mpc_out);