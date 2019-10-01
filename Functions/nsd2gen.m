function mpc = nsd2gen(mpc,ndc)
% nsd2gen is used to consider non-supplied demand through dispatchable loads
%   (negative generation). 
% 
%   There are some things taken from MATPOWER fcn 'load2disp'.
%
%% 
if nargin == 1
    ndcost = 1000;
elseif nargin == 2
    ndcost = ndc;
end

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

[~, cg] = size(mpc.gen);
mbase = mpc.baseMVA;

pd = mpc.bus(:,PD);
qd = mpc.bus(:,QD);

id_dem = find(pd > 0);
ndl = size(id_dem,1);
bus_dem = mpc.bus(id_dem,BUS_I);
% id_dem = find(id_dem==1);

gen_dem = zeros(ndl,cg);
gen_dem(:,GEN_BUS) = bus_dem;
gen_dem(:,PG) = -pd(id_dem);
gen_dem(:,PMIN) = -pd(id_dem);
gen_dem(:,VG) = 1;
gen_dem(:,GEN_STATUS) = 1;
gen_dem(:,MBASE) = mbase;

gen_dem(:,QG) = -qd(id_dem);


id_qpos = find(qd(id_dem) > 0);
id_qneg = find(qd(id_dem) < 0);

if ~isempty(id_qpos)
    gen_dem(id_qpos,QMIN) = -qd(id_dem(id_qpos));
end
if ~isempty(id_qneg)   
    gen_dem(id_qneg,QMAX) = -qd(id_dem(id_qneg));
end

[~,ccost] = size(mpc.gencost);

gencost_dem = zeros(length(bus_dem),ccost);
gencost_dem(:,1) = 2;
gencost_dem(:,NCOST) = 3;
gencost_dem(:,COST+1) = ndcost;
gencost_dem(:,COST+2) = ndcost*pd(id_dem);

mpc.gen = [mpc.gen; gen_dem];
mpc.gencost = [mpc.gencost; gencost_dem];

%% (optional) generator fuel types (taken from load2disp)
if isfield(mpc, 'genfuel') && iscell(mpc.genfuel)
    genfuel = cell(ndl, 1);
    for k = 1:ndl
        genfuel{k} = 'dl';
    end
    mpc.genfuel = [mpc.genfuel; genfuel];
end
%%
mpc.nsd.id_dem = id_dem;
mpc.nsd.original.PD = pd;
mpc.nsd.original.QD = qd;
mpc.nsd.N = ndl;
mpc.bus(:,PD) = 0;
mpc.bus(:,QD) = 0;
end