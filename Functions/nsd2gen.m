function mpc = nsd2gen(mpc0,nsdcost)
% nsd2gen Converts all loads to dispatchable in order to model non-supplied
%   demands through negegative generators.
%   MPC = NSD2GEN(MPC0);
%   MPC = NSD2GEN(MPC0,NDC);
%
%   Takes a MATPOWER case file or struct and converts fixed loads to
%   dispatchable loads and returns the resulting case struct. Inputs
%   are as follows:
%   
%   MPC0 - File name or struct with initial MATPOWER case.
%
%   NSDCOST (optional) - Scalar specifying the value of non-supplied demand 
%       to use as the value for the dispatchable loads. Default is $5000 
%       per MWh.
% 
%   There are some things taken from MATPOWER function LOAD2DISP.
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

%% default values
if nargin == 1
    ndcost = 5000;
elseif nargin == 2
    ndcost = nsdcost;
end
%% define constants
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%% initialize some things
mpc = mpc0;

[ng, cg] = size(mpc0.gen);
mbase = mpc0.baseMVA;

pd = mpc0.bus(:,PD);
qd = mpc0.bus(:,QD);

id_dem = find(pd > 0);
ndl = size(id_dem,1);
bus_dem = mpc0.bus(id_dem,BUS_I);
% id_dem = find(id_dem==1);
%% gen matrix
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
mpc.gen = [mpc0.gen; gen_dem];
%% gencost matrix
[~,ccost] = size(mpc0.gencost);

gencost_dem = zeros(length(bus_dem),ccost);
gencost_dem(:,1) = 2;
gencost_dem(:,NCOST) = 3;
gencost_dem(:,COST+1) = ndcost;
gencost_dem(:,COST+2) = ndcost*pd(id_dem);

mpc.gencost = [mpc0.gencost; gencost_dem];

%% (optional) generator fuel types (taken from load2disp)
if isfield(mpc0, 'genfuel') && iscell(mpc0.genfuel)
    genfuel = cell(ndl, 1);
    for k = 1:ndl
        genfuel{k} = 'dl';
    end
    mpc.genfuel = [mpc0.genfuel; genfuel];
end
%% organize final info
mpc.nsd.id_dem = id_dem;
mpc.nsd.original.PD = pd;
mpc.nsd.original.QD = qd;
mpc.nsd.N = ndl;
mpc.bus(:,PD) = 0;
mpc.bus(:,QD) = 0;
%% Create new ids for generators
mpc.genid.dl = ((ng+1):(ng+ndl))';
end