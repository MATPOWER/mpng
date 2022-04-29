function results = mpng(mpgc,mpopt)
% MPNG (MATPOWER - Natural Gas) Solves an optimal power and natural gas flow. 
%   RESULTS = mpng(MPGC,MPOPT)
% 
%   Runs an optimal power and natural gas flow and returns a RESULTS
%   struct. MPGC corresponds to the input struct case while MPOPT is the       
%   default MATPOWER options struct.
%
%   The input case (MPGC) must contain two additional fields, where a gas 
%   case (MPGC.mgc) and an interconnection case (MPGC.connect) must be 
%   included as follows: 
% 
%       MPGC.mgc: This structure must contain the information related to
%       the natural gas case, which has to have the following fields:
%                   
%           .version    current version of gas cases.
%           .fbase      base for the natural gas volumes, (MMSCFD)    
%           .pbase      base for the natural gas pressures,(PSI)
%           .wbase      base for the active power in the gas system,
%                       usually the same as in the power system, (MW) 
%           .node
%               .info   information related to the nodes of the natural 
%                       gas case, excluding costs.
%               .dem    matrix where the gas demand is specified for every
%                       node and every type.
%               .demcost matrix that specifies the cost of the non-supplied
%                        gas for every node and every type of demand.
%           .well       information related to the wells of the natural gas
%                       case, including extraction costs.
%           .pipe       information ralated to the pipelines of the natural
%                       gas case, including transport costs.
%           .comp       information related to the compressors of the
%                       natural gas case, including compression costs.
%           .sto        information related to the storage units of the
%                       natural gas system, including storage costs. 
% 
%       MPGC.connect: This struct contains the information related to the 
%       interconection case, which specifies how the power and natural gas
%       systems are related. It has the following fields:
% 
%           .version    current version for connect struct.
%           .power
%               .time   vector which specifies the number of hours for each
%                       period of time taking into considaration. Sum of  
%                       elements in the vector must equal 24.
%               .demands
%                   .pd     active power demand matrix for all buses and
%                           periods of time.
%                   .qd     reactive power demand matrix for all buses and
%                           periods of time.
%               .cost   general cost of non-supplied active power ($/MW)
%               .sr     spinning reserve matrix for all zones and periods
%                       of time. 
%               .energy a two-colum matrix for specifying the daily energy
%                       available for certain generators.
%           .interc
%                   .comp   a two-colum index matrix for locating compressors
%                           in specific buses. 
%                   .term   a three-column matrix for locating gas-fired 
%                           generators connected to specific gas nodes and
%                           fixing their corresponding efficiencies.
%
%   For a better undestanding on all contained information in this two
%   structs please see the manual and the following functions:
%
%   See also DEFINE_CONSTANTS_GAS MPC2GAS_PREP

%   MPNG: MATPOWER - Natural Gas
%   Copyright (c) 2019-2022 - v0.99beta
%   Sergio García-Marín - Universidad Nacional de Colombia - Sede Manizales
%   Wilson González-Vanegas - Universidad Tecnológica de Pereira
%   Carlos E. Murillo-Sánchez - Universidad Nacional de Colombia - Sede Manizales
% 
%   This file is part of MPNG.
%   Covered by the 3-clause BSD License (see LICENSE file for details).                                                   

%% check for Matpower options or set default
if nargin == 1
    mpopt = mpoption;
    mpopt.opf.ac.solver = 'IPOPT';     % set the default solver
    mpopt.ipopt.opts.max_iter = 1e5;   % default max number of iterations
end
nb = size(mpgc.bus,1);
if nb == 2
    mpopt.out.sys_sum = 0;
    mpopt.out.area_sum = 0;
    mpopt.out.bus = 0;
    mpopt.out.branch = 0;
    mpopt.out.lim.pg = 0;
    mpopt.out.lim.qg = 0;
%     mpopt.out.suppress_detail = 0;
end
verbose = mpopt.verbose;
%% check for proper gas inputs
if verbose
    fprintf(1,'MPNG: MATPOWER-Natural Gas - Version 0.99beta  \n')
    fprintf(1,'     - Thanks to UNAL-Mz & UTP \n\n')
end
if ~isfield(mpgc, 'connect') || ...
        ~isfield(mpgc, 'mgc') || ...
        ~isstruct(mpgc.connect) || ...
        ~isfield(mpgc.connect, 'power') || ...
        ~isfield(mpgc.connect, 'interc') 
    error('matpowergas: case must contain a ''connect'' and a ''matgas case'' field');
else %% add user functions
    if isfield(mpgc.mgc,'wey_percent')
        args.percent = mpgc.mgc.wey_percent;
    else
        args.percent = 0.2;
    end
    mpgc = add_userfcn(mpgc, 'ext2int', @userfcn_mpng_ext2int);
    mpgc = add_userfcn(mpgc, 'formulation', @userfcn_mpng_formulation,args);
    mpgc = add_userfcn(mpgc, 'int2ext', @userfcn_mpng_int2ext);
    mpgc = add_userfcn(mpgc, 'printpf', @userfcn_mpng_printpf);
    mpgc = add_userfcn(mpgc, 'savecase', @userfcn_mpng_savecase);
end

%% modify the power case for special features
if verbose
    fprintf(1,'Preparing case...\n')
end
mpgc = mpc2gas_prep(mpgc,mpopt);
%% run program
results = runopf(mpgc,mpopt);
end
function mpgc = userfcn_mpng_ext2int(mpgc, mpopt, args)
%
%   mpc = userfcn_mpng_ext2int(mpc, mpopt, args)
%
%   This is the 'ext2int' stage userfcn callback that prepares the input
%   data for the formulation stage. It expects to find a 'MGC' and a
%   'CONNECT' fields in mpc as described above. The optional args are
%	not currently used. 
%
%   Currently there is no change between external and internal indexing,
%   and status is not a consideration.

%% Define power and gas constants
% power
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
% [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
%     TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
%     ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
% [GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
%     MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
%     QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
% [PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
% [CT_LABEL, CT_PROB, CT_TABLE, CT_TBUS, CT_TGEN, CT_TBRCH, CT_TAREABUS, ...
%     CT_TAREAGEN, CT_TAREABRCH, CT_ROW, CT_COL, CT_CHGTYPE, CT_REP, ...
%     CT_REL, CT_ADD, CT_NEWVAL, CT_TLOAD, CT_TAREALOAD, CT_LOAD_ALL_PQ, ...
%     CT_LOAD_FIX_PQ, CT_LOAD_DIS_PQ, CT_LOAD_ALL_P, CT_LOAD_FIX_P, ...
%     CT_LOAD_DIS_P, CT_TGENCOST, CT_TAREAGENCOST, CT_MODCOST_F, ...
%     CT_MODCOST_X] = idx_ct;
% gas 
[DEM, WELL, NODE_I, NODE_TYPE, PR, PRMAX, PRMIN, OVP, UNP, COST_OVP, ...
    COST_UNP, GD, NGD] = idx_node;
[STO_NODE, STO, STO_0, STOMAX, STOMIN, FSTO, FSTO_OUT, FSTO_IN,...
    FSTOMAX, FSTOMIN, S_STATUS, COST_STO, COST_OUT, COST_IN] = idx_sto;
[F_NODE, T_NODE, FG_O, K_O, DIAM, LNG, FMAX_O, FMIN_O, COST_O] = idx_pipe;
[COMP_G, COMP_P, F_NODE, T_NODE, TYPE_C, RATIO_C, B_C, Z_C, ALPHA, BETA, ...
    GAMMA, FMAX_C, COST_C, FG_C, GC_C, PC_C] = idx_comp;
[WELL_NODE, G, PW, GMAX, GMIN, WELL_STATUS, COST_G] = idx_well;
%% initialize some things
mpc = mpgc;
mgc = mpgc.mgc;
connect = mpgc.connect;
% 
nt = size(connect.power.time,2);% number of periods of time
time = connect.power.time;
% power
nb = size(mpc.bus,1);
nl = size(mpc.branch,1);
ng = size(mpc.gen,1);
ndl = mpc.nsd.N;
nc_p = length(find(mgc.comp(:,TYPE_C) == COMP_P)); 
ng_pg_ori = length(mpc.genid.original)/nt;
areas = unique(mpc.bus(:,BUS_AREA));
nareas = length(areas);
% gas
nn = size(mgc.node.info,1);     % number of nodes
nw = size(mgc.well,1);          % number of wells
no = size(mgc.pipe,1);          % number of pipelines
nc = size(mgc.comp,1);          % number of compressors
ns = size(mgc.sto,1);           % number of storage units
ngd = mgc.node.info(1,NGD);     % number of gas demand types

%% check data for consistent dimensions
%% cas case
% total gas demand must be equal than the sum of the diversified gas in
% every node
dem_tot = mgc.node.info(:,GD);
dem_g = mgc.node.dem;
dem_diff = abs(dem_tot - sum(dem_g,2));
tol = 1e-3;
if any(dem_diff > tol) 
    error('mpng: Total gas demand does not correspont to the sum of the diversified demand.');
end
% matrix for gas demand must be nn x ngd
if any(size(dem_g) ~= [nn ngd])
    error('mpng: Diversified demand gas matrix has wrong dimenssions.');
end
% Ratio of compressors must be bigger than 1 and lower than 2
ratio_c =  mgc.comp(:,RATIO_C);
% if any(1 > ratio_c)
if any(0.25 > ratio_c)
    error('mpng: Compressors ratio must be greater than 1.');
% elseif any(ratio_c > 2)
elseif any(ratio_c > 4)
    error('mpng: Compressors ratio must be lower than 2.');
end
% Check if 'wey_percent' size is equal to the number of pipelines
if isfield(mgc,'wey_percent')
    [nweyp,cweyp] = size(mgc.wey_percent);
    if nweyp == 1 
        if cweyp ~= 1
             error('mpng: Weymouth percent vector dimension mismatch.');
        end
    else 
        if (cweyp ~= 1) || (nweyp ~= no)
             error('mpng: Weymouth percent vector dimension mismatch.');
        end
    end
end
% Check if initial storage is in the range of max and min storage
is_sto = ~isempty(mgc.sto);
if is_sto
    sto0 = mgc.sto(:,STO_0);
    sto_max = mgc.sto(:,STOMAX);
    sto_min = mgc.sto(:,STOMIN);
    sto_lim_l = (sto0 - sto_min) < -1e-4;
    sto_lim_u = (sto_max - sto0) < -1e-4;
    sto_lim_l_set_bound = (sto0 - sto_min) >= -1e-4 & (sto0 - sto_min) < 0;
    sto_lim_u_set_bound = (sto_max - sto0) >= -1e-4 & (sto_max - sto0) < 0;
    if any(sto_lim_l) || any(sto_lim_u)
        error('mpng: Initial storage is not inside the storage limits.');
    end
    if any(sto_lim_l_set_bound)
        mgc.sto(sto_lim_l_set_bound,STO_0) = mgc.sto(sto_lim_l_set_bound,STOMIN);
    end
    if any(sto_lim_u_set_bound)
        mgc.sto(sto_lim_u_set_bound,STO_0) = mgc.sto(sto_lim_u_set_bound,STOMAX);
    end
else
    mgc.sto = [ 1 	0   0   0   0	0   0   0   0   0   0   0   0   0];
end
%% connection case

% Spinning reserve asked for every period of time must be lower than total
%	demand.
sr = connect.power.sr;
if ~isempty(sr)
    [nsr, csr] = size(sr);
    if ~(nsr == nareas) || ~(csr == nt)
        error('mpng: Spinning reserves zones matrix has wrong dimensions.');
    end
end

% Unit commitment matrix dimensions
UC = mpgc.connect.power.UC;
if ~isempty(UC)
    [nuc, cuc] = size(UC);
    % UC matrix dimensions must correspont to ng and nt
    if ~(nuc == ng_pg_ori) || ~(cuc == nt)
        error('mpng: Unit commitment matrix has wrong dimensions.');
    end
    % UC matrix-values must be 0 or 1
    if any(~ismember(UC(:),[0 1]))
        error('mpng: Unit commitment matrix has wrong values.');
    end
end

% Ramp time matrix dimenssions
ramp_time = connect.power.ramp_time;
if ~isempty(ramp_time) && length(time)>1
    [nr_t, cr_t] = size(ramp_time);
    % ramp_time matrix dimensions must correspont to ng and nt
    if ~(nr_t == ng_pg_ori) || ~(cr_t == nt)
        error('mpng: Ramp time matrix has wrong dimensions.');
    end
    max_rt = max(connect.power.ramp_time,[],1);
    max_rt = max_rt(2:end);     % At the moment the firts column of ramp_time is unused       
    % ramp_time matrix values must be lower than the corresponding period of time
    if any(time(1:end-1)<max_rt)
        error('mpng: Ramp time matrix values are larger than times vector.');
    end
end
% Check max energy for power plants and matrix dimenssions 
%%
mgc = mgc_PU(mgc); 
mpgc.mgc = mgc;
end 
function om = userfcn_mpng_formulation(om, mpopt, args)
%
%   om = userfcn_mpng_formulation(om, mpopt, args)
%
%   This is the 'formulation' stage userfcn callback that defines the
%   user costs and constraints for the optimal power and natural gas flow. 
%   It expects to find a 'MGC' and a 'CONNECT' fields in mpc as described 
%   above. The optional args are not currently used.

%% Define power and gas constants
% power
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
[CT_LABEL, CT_PROB, CT_TABLE, CT_TBUS, CT_TGEN, CT_TBRCH, CT_TAREABUS, ...
    CT_TAREAGEN, CT_TAREABRCH, CT_ROW, CT_COL, CT_CHGTYPE, CT_REP, ...
    CT_REL, CT_ADD, CT_NEWVAL, CT_TLOAD, CT_TAREALOAD, CT_LOAD_ALL_PQ, ...
    CT_LOAD_FIX_PQ, CT_LOAD_DIS_PQ, CT_LOAD_ALL_P, CT_LOAD_FIX_P, ...
    CT_LOAD_DIS_P, CT_TGENCOST, CT_TAREAGENCOST, CT_MODCOST_F, ...
    CT_MODCOST_X] = idx_ct;
% gas 
[DEM, WELL, NODE_I, NODE_TYPE, PR, PRMAX, PRMIN, OVP, UNP, COST_OVP, ...
    COST_UNP, GD, NGD] = idx_node;
[STO_NODE, STO, STO_0, STOMAX, STOMIN, FSTO, FSTO_OUT, FSTO_IN,...
    FSTOMAX, FSTOMIN, S_STATUS, COST_STO, COST_OUT, COST_IN] = idx_sto;
[F_NODE, T_NODE, FG_O, K_O, DIAM, LNG, FMAX_O, FMIN_O, COST_O] = idx_pipe;
[COMP_G, COMP_P, F_NODE, T_NODE, TYPE_C, RATIO_C, B_C, Z_C, ALPHA, BETA, ...
    GAMMA, FMAX_C, COST_C, FG_C, GC_C, PC_C] = idx_comp;
[WELL_NODE, G, PW, GMAX, GMIN, WELL_STATUS, COST_G] = idx_well;
[GEN_ID, MAX_ENER, COMP_ID, BUS_ID, NODE_ID, EFF] = idx_connect;
%% 

mpc = om.get_mpc();
mgc = mpc.mgc;
connect = mpc.connect;

baseMVA = mpc.baseMVA;
fbase = mgc.fbase;
wbase = mgc.wbase;
pbase = mgc.pbase;

%% ------------------------------------------------------------------------
% Here we add the gas realated variables, costs and constraints
%  ------------------------------------------------------------------------
%% initialize some things
nt = size(connect.power.time,2);% number of periods of time
time = connect.power.time;
% power
nb = size(mpc.bus,1);
nl = size(mpc.branch,1);
ng = size(mpc.gen,1);
ndl = mpc.nsd.N;
is_power_tool = nb == 2;
% gas
nn = size(mgc.node.info,1);     % number of nodes
nw = size(mgc.well,1);          % number of wells
no = size(mgc.pipe,1);          % number of pipelines
nc = size(mgc.comp,1);          % number of compressors
ns = size(mgc.sto,1);           % number of storage units
ngd = mgc.node.info(1,NGD);     % number of gas demand types

from_o = mgc.pipe(:,F_NODE);
to_o = mgc.pipe(:,T_NODE);

from_c = mgc.comp(:,F_NODE);
to_c = mgc.comp(:,T_NODE);  


%% Gas production
g0 = mgc.well(:,G);
gmin = mgc.well(:,GMIN);
gmax = mgc.well(:,GMAX);
gcost = mgc.well(:,COST_G);

%% Pressure
p0 = mgc.node.info(:,PR);
% p0 = p0.^2;
pmin0 = zeros(length(p0),1);
pmin = mgc.node.info(:,PRMIN);
% pmin = pmin.^2;
pmax = mgc.node.info(:,PRMAX);
% pmax = pmax.^2;

ovp0 = zeros(length(p0),1);
% ovp0 = mgc.node.info(:,OVP);
ovpmin = zeros(length(p0),1);
ovpcost = mgc.node.info(:,COST_OVP);

unp0 = zeros(length(p0),1);
% unp0 = mgc.node.info(:,UNP);
unpmin = zeros(length(p0),1);
unpcost = mgc.node.info(:,COST_UNP);

Apr = eye(length(p0));
Aovp = -Apr;
Aunp = Apr;

Appos = [Apr Aovp];
Apneg = [Apr Aunp];

% % Pressures in compressors
% % (1) pj - pi >= 0
% ip = [(1:nc)'; (1:nc)'];
% jp = [from_c; to_c;];
% val1 = [-ones(nc,1); ones(nc,1)];
% Apcomp1 = sparse(ip,jp,val1,nc,nn);
% Lpcomp1 = zeros(nc,1);
% % (2) pj - ratio*pi <= 0
% ratio_c =  mgc.comp(:,RATIO_C);
% % ratio_c =  ratio_c;
% val2 = [-ratio_c; ones(nc,1)];
% Apcomp2 = sparse(ip,jp,val2,nc,nn);
% Upcomp2 = zeros(nc,1);

% Pressures in compressors
ratio_c =  mgc.comp(:,RATIO_C);
id_turbo = find(ratio_c < 1);

% (2) pj - ratio*pi <= 0
ip = [(1:nc)'; (1:nc)'];
jp = [from_c; to_c;];
val2 = [-ratio_c; ones(nc,1)];
Apcomp2 = sparse(ip,jp,val2,nc,nn);
Upcomp2 = zeros(nc,1);

% (1) pj - pi >= 0 (compressors) and  pi - pj >= 0 (turboexpanders)
jp([id_turbo id_turbo+nc]) = jp([id_turbo+nc id_turbo]);
val1 = [-ones(nc,1); ones(nc,1)];
Apcomp1 = sparse(ip,jp,val1,nc,nn);
Lpcomp1 = zeros(nc,1);

%% Non-supplied gas demand
dem = mgc.node.dem;                     % gas demand 
totaldem = mgc.node.info(:,GD);
isdem = (mgc.node.info(:,GD) ~= 0);     
idgamma = find(isdem);                  % node where demand is located
ndn = length(idgamma);                  % number of nodes with demand
ngamma = ngd*ndn;
gamma0 = zeros(ngamma,1);
gammamin = zeros(ngamma,1);
gammamax = (dem(idgamma,:))';
gammamax = gammamax(:);
gammacost = mgc.node.demcost;
gammacost = (gammacost(idgamma,:))';
gammacost = gammacost(:);

%% Storage

sto0 = mgc.sto(:,STO_0);                % inital storage level
stomax = mgc.sto(:,STOMAX);             % max storage capacity
stomin = mgc.sto(:,STOMIN);             % min storage capacity

sto_diff0 = mgc.sto(:,FSTO);
sto_diffmin = mgc.sto(:,FSTOMIN);
sto_diffmax = mgc.sto(:,FSTOMAX);

sto_out0 = mgc.sto(:,FSTO_OUT);
sto_outmin = zeros(ns,1);
sto_outmax = sto0 - stomin;        % abs is added to avoid lims problems

sto_in0 = mgc.sto(:,FSTO_IN);
sto_inmin = zeros(ns,1);
sto_inmax = stomax - sto0;         % abs is added to avoid lims problems

Asto_diff = eye(ns);
Asto_out = -eye(ns);
Asto_in = eye(ns);

Asto = [Asto_diff Asto_out Asto_in];
Lsto = zeros(ns,1);
Usto = zeros(ns,1);

sto_cost = mgc.sto(:,COST_STO); 
sto_out_cost = mgc.sto(:,COST_OUT);
sto_in_cost = mgc.sto(:,COST_IN);

k_sto_cost = sto_cost.*sto0;

%% Gas flow in pipeline

k_oij = mgc.pipe(:,K_O); % <------ [WGV]: See notes for K_O. Is it ok?

fgo0 = mgc.pipe(:,FG_O);
fgomax = mgc.pipe(:,FMAX_O);
fgomin = mgc.pipe(:,FMIN_O);

isfgopos = find(fgo0>=0);
fgopos0 = zeros(no,1);
fgopos0(isfgopos) = fgo0(isfgopos);
fgoposmax = 1.1*fgomax;
fgoposmin = zeros(no,1);

isfgoneg = find(fgo0<0);
fgoneg0 = zeros(no,1);
fgoneg0(isfgoneg) = fgo0(isfgoneg);
fgonegmax = zeros(no,1);
fgonegmin = 1.1*fgomin;

Afgo = [eye(no) -eye(no) -eye(no)];
Lfgo = zeros(no,1);
Ufgo = Lfgo;

fgocost = mgc.pipe(:,COST_O);

if numel(args.percent) == 1
    pi_star = ((args.percent*fgomax).^2)./(k_oij.^2);
elseif size(args.percent,2) == 1
    pi_star = args.percent;
end

parpipe.npipe = no;                 % \ 
parpipe.idfrom = from_o;            % | 
parpipe.idto = to_o;                % | <------ [WGV]: Create structure with pipeline parameters
parpipe.k_oij = k_oij;              % |
parpipe.other.nn = nn;              % |
parpipe.other.pi_star = pi_star;    % /

%% Gas flow in compressors
% 38.6e6/35/24/3600 * 100000  * (50/14.7)

% Gas-driven compressors
iscomp_g = (mgc.comp(:,TYPE_C) == COMP_G);      % compressors working with gas  
idcomp_g = find(iscomp_g);                      % id for gas-driven compressors
nc_g = length(idcomp_g);                        % number of compressors-gas

fgc0_g = mgc.comp(idcomp_g,FG_C);
fgcmax_g = mgc.comp(idcomp_g,FMAX_C);
fgcmin_g = zeros(nc_g,1);
fgccost_g = mgc.comp(idcomp_g,COST_C);

psi0_g = mgc.comp(idcomp_g,PC_C);               % must be treated apart
psimin_g = zeros(nc_g,1);

phi0  = mgc.comp(idcomp_g,GC_C);                % phi only applies in compressors-gas
phimin = zeros(nc_g,1);
 
% Power-driven compressors
if any(mgc.comp(:,TYPE_C) == COMP_P) &&...
        ~isempty(connect.interc.comp)
    iscomp_p = (mgc.comp(:,TYPE_C) == COMP_P);  % compressors working with power
    idcomp_p = find(iscomp_p);                  % id for power-driven compressors
    nc_p = length(idcomp_p);                    % number of compressors-power
    
    fgc0_p = mgc.comp(idcomp_p,FG_C);
    fgcmax_p = mgc.comp(idcomp_p,FMAX_C);
    fgcmin_p = zeros(nc_p,1);
    fgccost_p = mgc.comp(idcomp_p,COST_C);
    
    psi0_p = mgc.comp(idcomp_p,PC_C);               % must be treated apart
    psimin_p = zeros(nc_p,1);
    if any(mgc.comp(idcomp_p,RATIO_C) < 1)
        psimin_p(mgc.comp(idcomp_p,RATIO_C) < 1) = -inf;
    end
else
    iscomp_p = false(1);
    idcomp_p = find(iscomp_p);
    nc_p = length(idcomp_p);
end                                   

% parameters for compressors
parcomp.ncompg = nc_g;                % \ 
parcomp.ncompp = nc_p;                % |
parcomp.idfrom = from_c;              % | 
parcomp.idto = to_c;                  % | 
parcomp.x = mgc.comp(idcomp_g,ALPHA); % | 
parcomp.y = mgc.comp(idcomp_g,BETA);  % | 
parcomp.z = mgc.comp(idcomp_g,GAMMA); % | <------ Create structure with compressors parameters
parcomp.Bc = mgc.comp(:,B_C);         % | 
parcomp.Zc = mgc.comp(:,Z_C);         % | 
parcomp.other.nn = nn;                % |
parcomp.other.idcompg = idcomp_g;     % |
parcomp.other.idcompp = idcomp_p;     % /

%% Create new variables
om.add_var('g',nw,g0,gmin,gmax);
om.add_var('p',nn,p0,pmin0,[]);
om.add_var('ovp',nn,ovp0,ovpmin,[]);
om.add_var('unp',nn,unp0,unpmin,[]);
om.add_var('sto_diff',ns,sto_diff0,sto_diffmin,sto_diffmax);
om.add_var('sto_out',ns,sto_out0,sto_outmin,sto_outmax);
om.add_var('sto_in',ns,sto_in0,sto_inmin,sto_inmax);
om.add_var('fgo',no,fgo0,fgomin,fgomax);
om.add_var('fgopos',no,fgopos0,fgoposmin,fgoposmax);
om.add_var('fgoneg',no,fgoneg0,fgonegmin,fgonegmax);
om.add_var('gamma',ngamma,gamma0,gammamin,gammamax);
if any(iscomp_g)     % variables for gas-driven compressors
    om.add_var('fgc_g',nc_g,fgc0_g,fgcmin_g,fgcmax_g);
    om.add_var('psi_g',nc_g,psi0_g,psimin_g,[]);
    om.add_var('phi',nc_g,phi0,phimin,[]);
end
if any(iscomp_p)     % variables for power-driven compressors
    om.add_var('fgc_p',nc_p,fgc0_p,fgcmin_p,fgcmax_p);
    om.add_var('psi_p',nc_p,psi0_p,psimin_p,[]);
end

%% linear constraints
om.add_lin_constraint('p_pos', Appos, [], pmax, {'p', 'ovp'});    % overpressure
om.add_lin_constraint('p_neg', Apneg, pmin, [], {'p', 'unp'});    % underpressure
om.add_lin_constraint('sto', Asto, Lsto, Usto, {'sto_diff', 'sto_out','sto_in'}); % storage
om.add_lin_constraint('p_comp1', Apcomp1, Lpcomp1, [], {'p'});    % 
om.add_lin_constraint('p_comp2', Apcomp2, [], Upcomp2, {'p'});    % 
om.add_lin_constraint('fgo', Afgo, Lfgo, Ufgo, {'fgo','fgoneg','fgopos'});    % 

%% non-linear constraints
% -------------- [WGV] Nonlinear constraint: Equation (9) --------------
fcn_fpipe = @(x)fpipe_fcn(x, parpipe);
hess_fpipe = @(x, lambda)fpipe_hess(x, lambda, parpipe);
om.add_nln_constraint('f_pipe', no, 'true', fcn_fpipe, hess_fpipe, {'fgo','p'});

if any(iscomp_g)     % nln_constraint for gas-driven compressors
    % ------------ [WGV] Nonlinear constraint: Equation (11) --------------
    fcn_wcomp_g = @(x)wcompgas_fcn(x, parcomp);
    hess_wcomp_g = @(x, lambda)wcompgas_hess(x, lambda, parcomp);
    om.add_nln_constraint('w_comp_g', nc_g, '=', fcn_wcomp_g, hess_wcomp_g, {'psi_g','fgc_g','p'});

    % ------------ Gas consuption in gas-driven compressors Eq. (12) ------
    fcn_compgas = @(x)compgas_fcn(x, parcomp);
    hess_compgas = @(x, lam)compgas_hess(x, lam, parcomp);
    om.add_nln_constraint('g_comp_g', nc_g, '=', fcn_compgas, hess_compgas,{'phi','psi_g'});
end
% -------------------------------------------------------------------------
if any(iscomp_p)     % nln_constraint for power-driven compressors
    % ------------ [WGV] Nonlinear constraint: Equation (11) --------------
    fcn_wcomp_p = @(x)wcomppower_fcn(x, parcomp);
    hess_wcomp_p = @(x, lambda)wcomppower_hess(x, lambda, parcomp);
    om.add_nln_constraint('w_comp_p', nc_p, '=', fcn_wcomp_p, hess_wcomp_p, {'psi_p','fgc_p','p'});
end
%% cost
om.add_quad_cost('gcost',[],gcost,0,{'g'});
om.add_quad_cost('ovpcost',[],ovpcost,0,{'ovp'});
om.add_quad_cost('unpcost',[],unpcost,0,{'unp'});
om.add_quad_cost('gammacost',[],gammacost,0,{'gamma'});
om.add_quad_cost('sto_cost',[],-sto_cost,k_sto_cost,{'sto_diff'});
om.add_quad_cost('sto_out_cost',[],sto_out_cost,0,{'sto_out'});
om.add_quad_cost('sto_in_cost',[],-sto_in_cost,0,{'sto_in'});
om.add_quad_cost('fgoposcost',[],fgocost,0,{'fgopos'});
om.add_quad_cost('fgonegcost',[],-fgocost,0,{'fgoneg'});
if any(iscomp_g)  
    om.add_quad_cost('fgccost_g',[],fgccost_g,0,{'fgc_g'});
end
if any(iscomp_p)  
    om.add_quad_cost('fgccost_p',[],fgccost_p,0,{'fgc_p'});
end

%% Nodal balance equations

% Gas flow in pipelines
f_fgo = mgc.pipe(:,F_NODE);
t_fgo = mgc.pipe(:,T_NODE);
row_fgo = [f_fgo; t_fgo];
col_fgo = [(1:no)'; (1:no)'];
val_fgo = [-ones(no,1); ones(no,1)];
A_fgo = sparse(row_fgo,col_fgo,val_fgo,nn,no);

% Gas flow in compressors-power
f_fgc_p = mgc.comp(:,F_NODE); f_fgc_p = f_fgc_p(iscomp_p);
t_fgc_p = mgc.comp(:,T_NODE); t_fgc_p = t_fgc_p(iscomp_p);
row_fgc_p = [f_fgc_p; t_fgc_p];
col_fgc_p = [(1:nc_p)'; (1:nc_p)'];
val_fgc_p = [-ones(nc_p,1); ones(nc_p,1)];
A_fgc_p = sparse(row_fgc_p,col_fgc_p,val_fgc_p,nn,nc_p);

% Gas flow in compressors-gas
f_fgc_g = mgc.comp(:,F_NODE); f_fgc_g = f_fgc_g(iscomp_g);
t_fgc_g = mgc.comp(:,T_NODE); t_fgc_g = t_fgc_g(iscomp_g);
row_fgc_g = [f_fgc_g; t_fgc_g];
col_fgc_g = [(1:nc_g)'; (1:nc_g)'];
val_fgc_g = [-ones(nc_g,1); ones(nc_g,1)];
A_fgc_g = sparse(row_fgc_g,col_fgc_g,val_fgc_g,nn,nc_g);

% Gas consumed by compressors-gas
row_phi = f_fgc_g;
col_phi = (1:nc_g)';
val_phi = -ones(nc_g,1);
A_phi = sparse(row_phi,col_phi,val_phi,nn,nc_g);

% Gas wells
id_g = mgc.well(:,WELL_NODE);
row_g = id_g;
col_g = (1:nw)';
val_g = ones(nw,1);
A_g = sparse(row_g,col_g,val_g,nn,nw);

% Storage outflow
id_sto = mgc.sto(:,STO_NODE);            % nodes where storage is conected
row_sto = id_sto;
col_sto = (1:ns)';
val_sto = ones(ns,1);
A_sto = sparse(row_sto,col_sto,val_sto,nn,ns);

% Gas demanded by gas-fired units
ngtemp = size(mpc.gen,1) - ndl;
if ~isempty(connect.interc.term)    
    ng_ori = ngtemp/nt;
    gens = (1:ngtemp)';
    gens = reshape(gens,ng_ori,nt);
    id_pggas = connect.interc.term(:,GEN_ID);
    id_pgnode = connect.interc.term(:,NODE_ID);
    eff_pg = connect.interc.term(:,EFF);
    eff_pg = eff_pg/(fbase/baseMVA);
    
    row_pggas = id_pgnode*ones(1,nt);
    row_pggas = vec2mat(row_pggas,1);
    col_pggas = gens(id_pggas,:)';
    col_pggas = col_pggas(:);
    val_pggas = -(time'*eff_pg');
    val_pggas = val_pggas(:);
    A_pggas = sparse(row_pggas,col_pggas,val_pggas,nn,ng);
else
    A_pggas = zeros(nn,ng);
end 

% Non-supplied gas demands

row_gamma = (idgamma*ones(1,ngd))';
row_gamma = row_gamma(:);
col_gamma = (1:ngamma)';
val_gamma = ones(ngamma,1);
A_gamma = sparse(row_gamma,col_gamma,val_gamma,nn,ngamma);

% Bounds -> Demand in every node
dem_gas = mgc.node.info(:,GD);
L_node = dem_gas;
U_node = L_node;

if ~any(iscomp_g) && ~any(iscomp_p)
    A_node = [A_fgo A_g A_sto A_pggas A_gamma];
    om.add_lin_constraint('Nodes', A_node, L_node, U_node, {'fgo','g','sto_diff','Pg','gamma'});
elseif ~any(iscomp_g) && any(iscomp_p)
    A_node = [A_fgo A_fgc_p A_g A_sto A_pggas A_gamma];
    om.add_lin_constraint('Nodes', A_node, L_node, U_node, {'fgo','fgc_p','g','sto_diff','Pg','gamma'});
elseif any(iscomp_g) && ~any(iscomp_p)
    A_node = [A_fgo A_fgc_g A_phi A_g A_sto A_pggas A_gamma];
    om.add_lin_constraint('Nodes', A_node, L_node, U_node, {'fgo','fgc_g','phi','g','sto_diff','Pg','gamma'});
elseif any(iscomp_g) && any(iscomp_p)
    A_node = [A_fgo A_fgc_p A_fgc_g A_phi A_g A_sto A_pggas A_gamma];
    om.add_lin_constraint('Nodes', A_node, L_node, U_node, {'fgo','fgc_p','fgc_g','phi','g','sto_diff','Pg','gamma'});
end

%% ------------------------------------------------------------------------
% Power network additional constraints and cost
%  ------------------------------------------------------------------------
%% Power demanded by compressors must be iqual for all periods of time and
%   the first must be equal to psi_p
if ~is_power_tool
    ng = size(mpc.gen,1);
    if any(iscomp_p)
        idgencomp = mpc.genid.comp;
        % Power consumption must equal for all periods of time
        row = [(1:nc_p*(nt-1))'; (1:nc_p*(nt-1))'];
        col = reshape(idgencomp,[nc_p,nt]);
        col1 = (col(:,(1:end-1)))';
        col1 = reshape(col1,[numel(col1),1]);
        col2 = (col(:,(2:end)))';
        col2 = reshape(col2,[numel(col2),1]);
        col = [col1;col2];
        val = [ones(nc_p*(nt-1),1); -ones(nc_p*(nt-1),1)];
        
        Agencomp = sparse(row,col,val,nc_p*(nt-1),ng);
        Lgencomp = zeros(nc_p*(nt-1),1);
        Ugencomp = Lgencomp;
        
        om.add_lin_constraint('power_gencomp', Agencomp, Lgencomp, Ugencomp, {'Pg'});
        
        idgencomp_1 = reshape(idgencomp,[nc_p,nt]);
        idgencomp_1 = idgencomp_1(:,1);
        
        row = [(1:nc_p)';(1:nc_p)'];
        col = [(idgencomp_1); (ng+1:ng+nc_p)'];
        val = ones(nc_p*2,1);
        Agencomp_1 = sparse(row,col,val,nc_p,(ng+nc_p));
        Lgencomp_1 = zeros(numel(idgencomp_1),1);
        Ugencomp_1 = Lgencomp_1;
        
        om.add_lin_constraint('power_gencomp_1', Agencomp_1, Lgencomp_1, Ugencomp_1, {'Pg','psi_p'});
    end
    
    %% Spinning reserves related to zones
    if ~isempty(connect.power.sr)
        areas = mpc.bus(:,BUS_AREA);
        areas = sort(unique(areas));    % areas in the power network
        nareas = length(areas);         % number of areas in the power network
        
        sr = connect.power.sr;
        sr = sr/baseMVA;
        
        idxpgen = mpc.genid.original; % find real generators (non dl non comp)
        idxpgen = reshape(idxpgen,[length(idxpgen)/nt,nt]);
        busgen = mpc.gen(idxpgen(:,1),GEN_BUS);
        areasgen = mpc.bus(busgen,BUS_AREA);
        pgmax = mpc.gen(idxpgen(:,1),PMAX);
        pgmax = pgmax/baseMVA;
        
        for j = 1:nareas
            idx_area = areasgen==areas(j);
            pgmax_area(j,1) = sum(pgmax(idx_area));
        end
        
        U_sr = zeros(nt*nareas,1);
        L_sr = U_sr;
        A_sr = zeros(nt*nareas,ng);
        
        row_sr = areasgen;
        val_sr = ones(length(areasgen),1);
        for j = 1:nt
            U_sr(((j-1)*nareas)+1:j*nareas) = pgmax_area - sr(:,j);
            col_sr = idxpgen(:,j);
            A_sr(((j-1)*nareas)+1:j*nareas,:) = sparse(row_sr,col_sr,val_sr,nareas,ng);
        end
        om.add_lin_constraint('SR', A_sr, L_sr, U_sr, {'Pg'});
    end
    %% Energy available for hydro power generators
    if ~isempty(connect.power.energy)
        idg_ori = mpc.genid.original;
        ng_ori = length(idg_ori)/nt;        
        idg_ori = reshape(idg_ori,[ng_ori,nt]);
        ngenr = size(connect.power.energy,1);
        genr = connect.power.energy(:,GEN_ID);
        row = [];
        val = [];
        col = [];
        for i = 1:ngenr
            row = [row; i*ones(nt,1)];
            col = [col; idg_ori(genr(i),:)'];
            val = [val; time'];
        end
        A_ener = sparse(row,col,val,ngenr,ng);
        maxenergy = connect.power.energy(:,MAX_ENER);
        maxenergy = maxenergy/baseMVA;
        U_ener = maxenergy;
        L_ener = zeros(size(U_ener));
        om.add_lin_constraint('hydro_energy', A_ener, L_ener, U_ener, {'Pg'});    % Max hydro energy
    end
    %% Unit commitment of power generators
    if ~isempty(connect.power.UC)
        idg_ori = mpc.genid.original;
        UC = connect.power.UC;
        nUC = length(find(~UC));
        UC = UC(:);        
        row = (1:nUC)';
        col = idg_ori(~UC);
        val = ones(nUC,1);
        A_UC = sparse(row,col,val,nUC,ng);
        L_UC = zeros(nUC,1);
        U_UC = L_UC;
        om.add_lin_constraint('UCpg', A_UC, L_UC, U_UC, {'Pg'});
        om.add_lin_constraint('UCqg', A_UC, L_UC, U_UC, {'Qg'});
    end
    %% Generators ramp
    if ~isempty(connect.power.ramp_time) && length(time)>1
        idg_ori = mpc.genid.original;
        ng_ori = length(idg_ori);
        ng_pg_ori = ng_ori/nt;
        ramp_time = connect.power.ramp_time;
        ramp_time = ramp_time(:,2:end);
        ramp_time = ramp_time(:);
        ramp30 = mpc.gen(idg_ori,RAMP_30);
        ramp30 = ramp30(1:end-ng_pg_ori)./mpc.baseMVA;
        row = [(1:(ng_pg_ori*(nt-1)))'; (1:(ng_pg_ori*(nt-1)))'];
        col1 = idg_ori(ng_pg_ori+1:end);
        col2 = idg_ori(1:end-ng_pg_ori);
        col = [col1; col2];
        val = [ones(ng_pg_ori*(nt-1),1); -ones(ng_pg_ori*(nt-1),1)];        
        A_r_t = sparse(row,col,val,ng_pg_ori*(nt-1),ng);
        L_r_t = -2*ramp30.*ramp_time;
        U_r_t = 2*ramp30.*ramp_time;
        % Units with UC=0 do not have incoming ramp
        UC = connect.power.UC;
        if ~isempty(connect.power.UC) && size(UC,2)>1
            UC = UC(:,2:end);
            UC = UC(:);
            A_r_t = A_r_t(logical(UC),:);
            L_r_t = L_r_t(logical(UC),:);
            U_r_t = U_r_t(logical(UC),:);
        end
        % Add linear constraint
        om.add_lin_constraint('ramp', A_r_t, L_r_t, U_r_t, {'Pg'});
    end
end
end
function results = userfcn_mpng_int2ext(results, mpopt,args)
%
%   results = userfcn_mpng_int2ext(results, mpopt, args)
%
%   This is the 'int2ext' stage userfcn callback that packages up the 
%   results. It expects to find a 'MGC' and a 'CONNECT' fields in mpc 
%   as described above. The optional args are not currently used.
%
%   Currently there is no change between external and internal indexing,
%   and status is not a consideration.

%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[DEM, WELL, NODE_I, NODE_TYPE, PR, PRMAX, PRMIN, OVP, UNP, COST_OVP, ...
    COST_UNP, GD, NGD] = idx_node;
[STO_NODE, STO, STO_0, STOMAX, STOMIN, FSTO, FSTO_OUT, FSTO_IN,...
    FSTOMAX, FSTOMIN, S_STATUS, COST_STO, COST_OUT, COST_IN] = idx_sto;
[F_NODE, T_NODE, FG_O, K_O, DIAM, LNG, FMAX_O, FMIN_O, COST_O] = idx_pipe;
[COMP_G, COMP_P, F_NODE, T_NODE, TYPE_C, RATIO_C, B_C, Z_C, ALPHA, BETA, ...
    GAMMA, FMAX_C, COST_C, FG_C, GC_C, PC_C] = idx_comp;
[WELL_NODE, G, PW, GMAX, GMIN, WELL_STATUS, COST_G] = idx_well;

%% organize info
mvabase = results.baseMVA;

mgc = results.mgc;
mgc = mgc_REAL(mgc);
connect = results.connect;
fbase = mgc.fbase;
wbase = mgc.wbase;
pbase = mgc.pbase;

nn = size(mgc.node.info,1);     % number of nodes
nw = size(mgc.well,1);          % number of wells
no = size(mgc.pipe,1);          % number of pipelines
nc = size(mgc.comp,1);          % number of compressors
ns = size(mgc.sto,1);           % number of storage units
nt = size(connect.power.time,2);% number of periods of time
ngamma = mgc.node.info(1,NGD);
time = connect.power.time;

from_o = mgc.pipe(:,F_NODE);
to_o = mgc.pipe(:,T_NODE);

from_c = mgc.comp(:,F_NODE);
to_c = mgc.comp(:,T_NODE);  

iscomp_g = (mgc.comp(:,TYPE_C) == COMP_G); 
idcomp_g = find(iscomp_g);                      % id for gas-driven compressors
nc_g = length(idcomp_g); 

if any(mgc.comp(:,TYPE_C) == COMP_P) &&...
        ~isempty(connect.interc.comp)
    iscomp_p = (mgc.comp(:,TYPE_C) == COMP_P);
    idcomp_p = find(iscomp_p);                  % id for gas-driven compressors
    nc_p = length(idcomp_p);
else
    iscomp_p = false(1);
    idcomp_p = find(iscomp_p);
    nc_p = length(idcomp_p);
end
%% Nodes
p = sqrt(results.var.val.p)*pbase;
ovp = p -(mgc.node.info(:,PRMAX));
unp = (mgc.node.info(:,PRMIN)) - p;
mgc.node.info(:,PR) = p;
mgc.node.info(:,OVP) = ovp;
mgc.node.info(:,UNP) = unp;

%% Well
g = results.var.val.g*fbase;
mgc.well(:,G) = g;

%% Pipeline
fgo = (results.var.val.fgo)*fbase;
mgc.pipe(:,FG_O) = fgo;

%% Compressor

if ~any(iscomp_g) && ~any(iscomp_p)
    fgc_p = 0;
    psi_p = 0;
    
    fgc_g = 0;
    psi_g = 0;  
elseif ~any(iscomp_g) && any(iscomp_p)
    fgc_p = (results.var.val.fgc_p)*fbase;
    psi_p = (results.var.val.psi_p)*wbase;
    
    fgc_g = 0;
    psi_g = 0;
elseif any(iscomp_g) && ~any(iscomp_p)
    fgc_p = 0;
    psi_p = 0;
    
    fgc_g = (results.var.val.fgc_g)*fbase;
    psi_g = (results.var.val.psi_g)*wbase;
    phi = (results.var.val.phi)*fbase;
elseif any(iscomp_g) && any(iscomp_p)
    fgc_p = (results.var.val.fgc_p)*fbase;
    psi_p = (results.var.val.psi_p)*wbase;
    
    fgc_g = (results.var.val.fgc_g)*fbase;
    psi_g = (results.var.val.psi_g)*wbase;
    phi = (results.var.val.phi)*fbase;
end

if any(iscomp_p)
    mgc.comp(iscomp_p,FG_C) = fgc_p;
    mgc.comp(iscomp_p,PC_C) = psi_p;
    mgc.comp(iscomp_p,GC_C) = 0;
end 

if any(iscomp_g)
    mgc.comp(iscomp_g,FG_C) = fgc_g;
    mgc.comp(iscomp_g,PC_C) = psi_g;
    mgc.comp(iscomp_g,GC_C) = phi;
end
%% Storage 
sto_diff = results.var.val.sto_diff*fbase;
sto_out = results.var.val.sto_out*fbase;
sto_in = results.var.val.sto_in*fbase;

mgc.sto(:,STO) = mgc.sto(:,STO_0) - sto_diff;
mgc.sto(:,FSTO) = sto_diff;
mgc.sto(:,FSTO_OUT) = sto_out;
mgc.sto(:,FSTO_IN) = sto_in;
%% extract individual periods of time for power network
mpc.version = results.version;
mpc.baseMVA = results.baseMVA;
mpc.bus = results.bus;
mpc.gen = results.gen;
mpc.branch = results.branch;
mpc.gencost = results.gencost;
mpc.f = results.f;
mpc.success = results.success;
results.multi_period.periods = extract_islands(mpc);
for i = 1:length(results.multi_period.periods)
    results.multi_period.periods{i}.et = 0;
end
%%
results.mgc = mgc;
end
function results = userfcn_mpng_printpf(results, fd, mpopt, args)
%
%   results = userfcn_mpng_printpf(results, fd, mpopt, args)
%
%   This is the 'printpf' stage userfcn callback that pretty-prints the
%   results. It expects a results struct, a file descriptor and a MATPOWER
%   options struct. The optional args are not currently used.

%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[DEM, WELL, NODE_I, NODE_TYPE, PR, PRMAX, PRMIN, OVP, UNP, COST_OVP, ...
    COST_UNP, GD, NGD] = idx_node;
[STO_NODE, STO, STO_0, STOMAX, STOMIN, FSTO, FSTO_OUT, FSTO_IN,...
    FSTOMAX, FSTOMIN, S_STATUS, COST_STO, COST_OUT, COST_IN] = idx_sto;
[F_NODE, T_NODE, FG_O, K_O, DIAM, LNG, FMAX_O, FMIN_O, COST_O] = idx_pipe;
[COMP_G, COMP_P, F_NODE, T_NODE, TYPE_C, RATIO_C, B_C, Z_C, ALPHA, BETA, ...
    GAMMA, FMAX_C, COST_C, FG_C, GC_C, PC_C] = idx_comp;
[WELL_NODE, G, PW, GMAX, GMIN, WELL_STATUS, COST_G] = idx_well;

%% organize info
baseMVA = results.baseMVA;

mgc = results.mgc;
% mgc = mgc_REAL(mgc);
connect = results.connect;
fbase = mgc.fbase;
wbase = mgc.wbase;
pbase = mgc.pbase;

nn = size(mgc.node.info,1);     % number of nodes
nw = size(mgc.well,1);          % number of wells
no = size(mgc.pipe,1);          % number of pipelines
nc = size(mgc.comp,1);          % number of compressors
ns = size(mgc.sto,1);           % number of storage units
nt = size(connect.power.time,2);% number of periods of time
ngamma = mgc.node.info(1,NGD);
time = connect.power.time;

from_o = mgc.pipe(:,F_NODE);
to_o = mgc.pipe(:,T_NODE);

from_c = mgc.comp(:,F_NODE);
to_c = mgc.comp(:,T_NODE);  

iscomp_g = (mgc.comp(:,TYPE_C) == COMP_G); 
idcomp_g = find(iscomp_g);                      % id for gas-driven compressors
nc_g = length(idcomp_g); 

if any(mgc.comp(:,TYPE_C) == COMP_P) &&...
        ~isempty(connect.interc.comp)
    iscomp_p = (mgc.comp(:,TYPE_C) == COMP_P);
    idcomp_p = find(iscomp_p);                  % id for gas-driven compressors
    nc_p = length(idcomp_p);
else
    iscomp_p = false(nc,1);
    idcomp_p = find(iscomp_p);
    nc_p = length(idcomp_p);
end

%% data
% wells
idwell = mgc.well(:,WELL_NODE);
g = mgc.well(:,G);

% nodes
nodes_id = mgc.node.info(:,1);
p = mgc.node.info(:,PR);
ovp = mgc.node.info(:,OVP);
unp = mgc.node.info(:,UNP);
demg = mgc.node.info(:,GD);


isdem = (mgc.node.info(:,GD) ~= 0);     
idgamma = find(isdem); 
gamma = results.var.val.gamma*fbase;
gamma = vec2mat(gamma,ngamma);
gamma_nodes = sum(gamma,2);
nondem_tot = sum(results.var.val.gamma*fbase); % add demands from gas-fired units and compressors

mu_nodes = (results.lin.mu.l.Nodes-results.lin.mu.u.Nodes)/fbase; % lambda nodal

% pipelines
k_o = mgc.pipe(:,K_O);
fgo_max = mgc.pipe(:,FMAX_O);
fgo = mgc.pipe(:,FG_O);


% compressors
comp_ratio = p(mgc.comp(:,T_NODE))./p(mgc.comp(:,F_NODE));
if ~any(iscomp_g) && ~any(iscomp_p)
    fgc_p = 0;
    psi_p = 0;
    fgccost_p = 0;
    
    fgc_g = 0;
    psi_g = 0;  
    fgccost_g = 0;
elseif ~any(iscomp_g) && any(iscomp_p)
    fgc_p = mgc.comp(iscomp_p,FG_C);
    psi_p = mgc.comp(iscomp_p,PC_C);
    fgccost_p = results.qdc.fgccost_p;
    
    fgc_g = 0;
    psi_g = 0;
    fgccost_g = 0;
elseif any(iscomp_g) && ~any(iscomp_p)
    fgc_p = 0;
    psi_p = 0;
    fgccost_p = 0;
    
    fgc_g = mgc.comp(iscomp_g,FG_C);
    psi_g = mgc.comp(iscomp_g,PC_C);
    phi = mgc.comp(iscomp_g,GC_C);
    fgccost_g = results.qdc.fgccost_g;
elseif any(iscomp_g) && any(iscomp_p)
    fgc_p = mgc.comp(iscomp_p,FG_C);
    psi_p = mgc.comp(iscomp_p,PC_C);
    fgccost_p = results.qdc.fgccost_p;
    
    fgc_g = mgc.comp(iscomp_g,FG_C);
    psi_g = mgc.comp(iscomp_g,PC_C);
    phi = (results.var.val.phi)*fbase;
    fgccost_g = results.qdc.fgccost_g;
end

% storage
sto_node    = mgc.sto(:,STO_NODE);
sto         = mgc.sto(:,STO);
sto_0       = mgc.sto(:,STO_0);
sto_max     = mgc.sto(:,STOMAX);
fsto        = mgc.sto(:,FSTO);
fsto_out	= mgc.sto(:,FSTO_OUT);
fsto_in     = mgc.sto(:,FSTO_IN);
sto_tot     = sum(sto);

sto_cost = results.qdc.sto_cost;
outflow_cost = results.qdc.sto_out_cost;
inflow_cost = results.qdc.sto_in_cost;

% cost
gcost = results.qdc.gcost;
gcost_tot = sum(gcost);
ovpcost = results.qdc.ovpcost;
unpcost = results.qdc.unpcost;
gammacost = results.qdc.gammacost;
sto_cost = results.qdc.sto_cost;
sto_out_cost = results.qdc.sto_out_cost;
sto_in_cost = results.qdc.sto_in_cost;
fgocost = (results.qdc.fgonegcost + results.qdc.fgoposcost);

% Gas-fired units
if ~isempty(connect.interc.term)
    id_gfu = connect.interc.term(:,1);
    node_gfu = connect.interc.term(:,2);
    eff_gfu = connect.interc.term(:,3);
    pg = results.var.val.Pg*baseMVA;
    ndl = results.nsd.N;
    pg_gfu_day = pg(1:end-ndl);
    pg_gfu_day = reshape(pg_gfu_day,length(pg_gfu_day)/nt,nt);
    pg_gfu_day = pg_gfu_day(id_gfu,:);
    pg_gfu_day = pg_gfu_day*time';
    g_gfu = pg_gfu_day.*eff_gfu;
end

%% Power system
% Non-supplied demand
ndl = results.nsd.N; 
pd =  results.nsd.original.PD;
pd = reshape(pd,length(pd)/nt,nt);
id_dem = find(pd(:,1) ~= 0);
pd = pd(id_dem,:);

e_demanded = pd*time';
d_supplied = -results.gen(end-ndl+1:end,PG);
d_supplied = reshape(d_supplied,ndl/nt,nt);
e_supplied = d_supplied*time';

e_non = e_demanded - e_supplied;

%% print options
isOPF           = isfield(results, 'f') && ~isempty(results.f);
SUPPRESS        = mpopt.out.suppress_detail;
OUT_ALL         = mpopt.out.all;
OUT_FORCE       = mpopt.out.force;
OUT_RES         = OUT_ALL == 1 || (OUT_ALL == -1 && ~SUPPRESS && (mpopt.out.bus || mpopt.out.gen));

%%
nb = size(results.bus(:,1));
is_power_tool = nb == 2; 
if ~is_power_tool
fprintf(fd, '\n================================================================================');
fprintf(fd, '\n|     Non-supplied power demand (MWh/d)                                        |');
fprintf(fd, '\n================================================================================');
fprintf(fd, '\n  Bus     Total    Supplied     Non-S.                                     ');
fprintf(fd, '\n   #      Energy    Energy      Energy                               ');
fprintf(fd, '\n -----   --------  ---------   --------                              ');
for i = 1:length(id_dem)
    fprintf(fd,'\n  %3d',id_dem(i));
    fprintf(fd,'   %9.2f',e_demanded(i));
    fprintf(fd,'  %9.2f',e_supplied(i));
    fprintf(fd,' %9.2f',e_non(i));
end
fprintf(fd, '\n -----  ---------  ---------   --------                              ');
fprintf(fd, '\n Total  %9.2f  %9.2f %9.2f',sum(e_demanded),sum(e_supplied),sum(e_non));
end
%%
fprintf(fd, '\n================================================================================');
fprintf(fd, '\n|     Gas system summary                                                       |');
fprintf(fd, '\n================================================================================');
fprintf(fd, '\n\nHow many?                  How much?              ');
fprintf(fd, '\n---------------------    -------------------------------- ');
fprintf(fd, '\nNodes          %3d         Total Well Capacity    %2.2f ', nn, sum(mgc.well(:, GMAX)));
fprintf(fd, '\nWells          %3d         On-line Capacity       %2.2f', nw, sum(mgc.well(:, GMAX)));
fprintf(fd, '\nPipelines      %3d         Gas Production         %2.2f', no,  sum(g));
fprintf(fd, '\nCompressors    %3d         Total Demand           %2.2f', nc, sum(mgc.node.info(:, GD)));
fprintf(fd, '\n  Gas Comp.    %3d           Supplied Demand      %2.2f', nc_g, sum(mgc.node.info(:, GD))-nondem_tot);
fprintf(fd, '\n  Power Comp.  %3d           Non-Supplied Demand  %2.2f', nc_p,nondem_tot);
fprintf(fd, '\nStorage Units  %3d         Gas Stored             %2.2f', ns, sto_tot);
fprintf(fd, '\n')

fprintf(fd, '\nGas total extraction cost =   %5.2f ($/day)',gcost_tot);

fprintf(fd, '\n=================================================================================');
fprintf(fd, '\n|     Nodes Data                                                                |');
fprintf(fd, '\n=================================================================================');
fprintf(fd, '\n Node   Pressure   Over      Under     Demand    Non-S       Gas        Nodal    ');
fprintf(fd, '\n  #               Pressure  Pressure             Demand   Extraction    Lambda   ');
fprintf(fd, '\n         (psi)     (psi)     (psi)    (MMSCFD)  (MMSCFD)   (MMSCFD)    ($/MMSCFD)');
fprintf(fd, '\n ----   --------  --------  --------  --------  --------  -----------  ----------');
k = 1;
for i = 1:nn
    fprintf(fd,'\n  %2d',nodes_id(i));
    fprintf(fd,'   %9.3f',p(i));    
    
    if ovp(i) <= 0
        fprintf(fd,'     ---  ');
    else
        fprintf(fd,' %8.2f ',ovp(i));
    end
    
    if unp(i) <= 0
        fprintf(fd,'     ---  ');
    else
        fprintf(fd,' %8.2f ',unp(i));        
    end    
    isdem = (mgc.node.info(:,GD) ~= 0);    
    if ~isdem(i)  
        fprintf(fd,'     ---  ');
        fprintf(fd,'     ---  ');
    else
        fprintf(fd,'  %7.2f ',demg(i));
        fprintf(fd,'  %7.2f ',gamma_nodes(k));
        k = k+1;
    end
    
    iswell = find(idwell == i);
    gnode = sum(g(iswell)); 
    if ~any(iswell)
        fprintf(fd,'       ---  ');
    else
        fprintf(fd,'    %7.2f ',gnode);
    end
    fprintf(fd,'   %9.2f ',mu_nodes(i));
end
fprintf(fd, '\n');
fprintf(fd, '\n================================================================================');
fprintf(fd, '\n|     Pipeline Data                                                            |');
fprintf(fd, '\n================================================================================');
fprintf(fd, '\n Pipeline    From      To       Weymouth      Max Gas       Gas       Transport');
fprintf(fd, '\n    #        Node     Node      Constant        Flow        Flow        Cost   ');
fprintf(fd, '\n                              (MMSCFD/psia)   (MMSCFD)    (MMSCFD)     ($/day) ');
fprintf(fd, '\n --------    -----    -----   -------------   --------    --------    ---------');
for i = 1:no
    fprintf(fd,'\n   %2d',i);
    fprintf(fd,'         %2d',from_o(i));
    fprintf(fd,'       %2d',to_o(i));
    fprintf(fd,'        %7.3f',k_o(i));
    fprintf(fd,'      %7.2f',fgo_max(i));
    fprintf(fd,'     %7.2f',fgo(i));
    fprintf(fd,'     %8.2f ',fgocost(i));
end
fprintf(fd, '\n');
fprintf(fd, '\n================================================================================');
fprintf(fd, '\n|     Compressor Data                                                          |');
fprintf(fd, '\n================================================================================');
fprintf(fd, '\n Comp.   Comp.   From    To     Comp.      Power       Gas      Comp.    Comp.  ');
fprintf(fd, '\n   #     Type    Node   Node    Flow      Consumed   Consumed   Ratio    Cost   ');
fprintf(fd, '\n                               (MMSCFD)   (MW/day)   (MMSCFD)           ($/day) ');
fprintf(fd, '\n -----   -----   ----   ----   --------   ---------  --------   -----   --------');
j = 1; k = 1;
for i = 1:nc
    fprintf(fd,'\n  %2d ',i);
    if iscomp_p(i)
        fprintf(fd,'      P');
        fprintf(fd,'      %2d',from_c(i));
        fprintf(fd,'     %2d',to_c(i));
        fprintf(fd,'    %8.3f',fgc_p(j));
        fprintf(fd,'   %7.2f',psi_p(j));
        fprintf(fd,'       ---  ');
        fprintf(fd,'   %5.3f',comp_ratio(i));
        fprintf(fd,'    %5.2f',fgccost_p(j));
        j = j + 1;
    end
    if iscomp_g(i)
        fprintf(fd,'      G');
        fprintf(fd,'      %2d',from_c(i));
        fprintf(fd,'     %2d',to_c(i));
        fprintf(fd,'    %8.3f',fgc_g(k));
        fprintf(fd,'   %7.2f',psi_g(k));
        fprintf(fd,'     %6.4f',phi(k));
        fprintf(fd,'    %5.3f',comp_ratio(i));
        fprintf(fd,'    %5.2f',fgccost_g(k));
        k = k + 1;
    end    

end
fprintf(fd, '\n');

sto_max_tot = sum(sto_max);
exist_sto = sto_max_tot ~= 0;
if exist_sto
    fprintf(fd, '\n====================================================================================');
    fprintf(fd, '\n|     Storage Units Data                                                           |');
    fprintf(fd, '\n====================================================================================');
    fprintf(fd, '\n Storage  Initial   Final    Sto.      Sto.      Sto.      Sto.    Outflow   Inflow ');
    fprintf(fd, '\n   Node     Sto.     Sto.    Diff.     Out       In        Cost     Cost      Cost  ');
    fprintf(fd, '\n          (MMSCF)  (MMSCF)  (MMSCFD)  (MMSCFD)  (MMSCFD)  ($/day)  ($/day)   ($/day)');
    fprintf(fd, '\n -------  -------  -------  --------  --------  --------  -------  -------  --------');
    for i = 1:ns
        fprintf(fd,'\n    %2d ',sto_node(i));
        fprintf(fd,'    %5.1f',sto_0(i));
        fprintf(fd,'    %5.1f',sto(i));
        fprintf(fd,'    %6.2f',fsto(i));
        fprintf(fd,'    %6.2f',fsto_out(i));
        fprintf(fd,'    %6.2f',fsto_in(i));
        fprintf(fd,'   %7i',round(sto_cost(i)));
        fprintf(fd,'   %6i',round(outflow_cost(i)));
        fprintf(fd,'    %6i',round(inflow_cost(i)));
    end
    fprintf(fd, '\n');
end 
if ~isempty(connect.interc.term)
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n|     Gas-fired Generators Data                                                |');
    fprintf(fd, '\n================================================================================');
    fprintf(fd, '\n Unit.   Node   Plant         Daily       Gas     ');
    fprintf(fd, '\n   #      #      Eff.         Energy     Consumed ');
    fprintf(fd, '\n               (MMSCFD/MW)   (MWh/day)   (MMSCFD) ');
    fprintf(fd, '\n -----   ----  -----------   ---------   -------- ');
    for i = 1:length(id_gfu)
        fprintf(fd,'\n  %2d',id_gfu(i));
        fprintf(fd,'      %2d',node_gfu(i));
        fprintf(fd,'     %4.2e',eff_gfu(i));
        fprintf(fd,'   %9.2f',pg_gfu_day(i));
        fprintf(fd,'    %7.2f',g_gfu(i));
    end
    fprintf(fd, '\n -----   ----  -----------   ---------   --------                           ');
    fprintf(fd, '\n                     Total  %9.2f    %7.2f   ',sum(pg_gfu_day),sum(g_gfu));
    fprintf(fd, '\n');
end
end
function mpgc = userfcn_mpng_savecase(mpgc, fd, prefix, args)
%
%   This is the 'savecase' stage userfcn callback that prints the M-file
%   code to save the results of the optminal power and natural gas flow.
%
%   Not currently used.
%
end
%% -------------- [WGV] Jacobian and Hessian for equation (12) --------------
function [g, dg] = compgas_fcn(x, parcomp)
%
%
%
%%
alpha = parcomp.x;
beta = parcomp.y;
gamma = parcomp.z;
n = parcomp.ncompg;
[phi, psi] = deal(x{:});
% psi = x(1:n);
% phi = x(n+1:2*n);

if nargout == 1
    g = phi - alpha - diag(beta)*psi - diag(gamma)*(psi.^2);
else
    g = phi - alpha - diag(beta)*psi - diag(gamma)*(psi.^2);
    m1 = eye(n);
    m2 = -diag(beta) - 2*diag(gamma)*diag(psi);
    dg = [m1 m2];
end
end
function d2H = compgas_hess(x, lambda, parcomp)
%
%  
gamma = parcomp.z;
n = parcomp.ncompg;
% [phi, psi] = deal(x{:});

m0 = zeros(n);
m1 = -2*eye(n) * diag(gamma) * diag(lambda);

d2H = [ m0 m0;
        m0 m1];
end 

%% -------------- [WGV] Jacobian and Hessian for equation (9) --------------
function [g, dg] = fpipe_fcn(x, parpipe)
npipe = parpipe.npipe;   % Extract pipeline parameters
from_o = parpipe.idfrom;
to_o = parpipe.idto;     
k_oij = parpipe.k_oij;
k_oij = k_oij.^2;
nn = parpipe.other.nn;
pi_star = parpipe.other.pi_star;

[fgo,p] = deal(x{:});   % Extract variables 

[f,df] = wey_approx(k_oij,p(from_o),p(to_o),pi_star); % Eval Weymouth equation approximation

if nargout < 2
    g = fgo - f; % Eval nonlinear constraints
else
    g = fgo - f; % Eval nonlinear constraints
    jac1 = speye(npipe); % Eval dg/d(fgo)
    
    id_jac2   = [(1:npipe)' from_o    % Indexes for dg/d(pi)
                 (1:npipe)'  to_o  ]; % Indexes for dg/d(pj)           
    vals_jac2 = [-df                  % Values for dg/d(pi)
                  df];                % Values for dg/d(pj)
        jac2 = sparse(id_jac2(:,1),id_jac2(:,2),vals_jac2,npipe,nn); % Eval dg/dp
    dg = [jac1 jac2]; % Compute Jacobian of the Lagrangian   
end
end
function d2L = fpipe_hess(x, lambda, parpipe)
npipe = parpipe.npipe;   % Extract pipeline parameters
from_o = parpipe.idfrom;
to_o = parpipe.idto;     
k_oij = parpipe.k_oij;
k_oij = k_oij.^2;
nn = parpipe.other.nn;
pi_star = parpipe.other.pi_star;

[~,p] = deal(x{:});   % Extract variables

hess1 = spalloc(npipe,npipe,npipe*npipe); % Eval d^2(L)/d^2(fgo)

[~,~,d2f] = wey_approx(k_oij,p(from_o),p(to_o),pi_star); % Eval Weymouth equation approximation

id_hess2 = [from_o  from_o   % Indexes for d^2(L)/d^2(pi)
            to_o    to_o     % Indexes for d^2(L)/d^2(pj)
            from_o  to_o     % Indexes for d^2(L)/d(pi)d(pj)
            to_o    from_o]; % Indexes for d^2(L)/d(pj)d(pi)        
vals_hess2 = [-lambda.*d2f   % Values for d^2(L)/d^2(pi) 
               lambda.*d2f   % Values for d^2(L)/d^2(pj) 
               lambda.*d2f   % Values for d^2(L)/d(pi)d(pj) 
               lambda.*d2f]; % Values for d^2(L)/d(pj)d(pi) 

% Option 1: using sparse Matlab function
hess2 = sparse(id_hess2(:,1),id_hess2(:,2),vals_hess2,nn,nn);
d2L = [    hess1        spalloc(npipe,nn,npipe*nn) % Compute Hessian of the Lagrangian
       spalloc(nn,npipe,nn*npipe)       hess2    ];
end

%% -------------- [WGV] Jacobian and Hessian for equation (11) --------------
function [g, dg] = wcompgas_fcn(x, parcomp)
nc_g = parcomp.ncompg;  % Extract compressor parameters
from_c = parcomp.idfrom;
to_c = parcomp.idto;
nn = parcomp.other.nn;
id_comp_gas = parcomp.other.idcompg;
Bc = parcomp.Bc; Bc = Bc(id_comp_gas);
Zc = parcomp.Zc; Zc = Zc(id_comp_gas);

[psi,fgc,p] = deal(x{:}); % Extract variables

if nargout < 2
    pi = p(from_c(id_comp_gas));
    pj = p(to_c(id_comp_gas));
    g = psi - Bc.*fgc.*((pj./pi).^(0.5*Zc) - 1); % Eval nonlinear constraint
else
    pi = p(from_c(id_comp_gas));
    pj = p(to_c(id_comp_gas));
    g = psi - Bc.*fgc.*((pj./pi).^(0.5*Zc) - 1); % Eval nonlinear constraint
    
    jac1 = speye(nc_g); % Eval dg/d(psi)
    
    id_jac2   = [(1:nc_g)' (1:nc_g)'];         % Indexes for dg/d(fgc)                            
    vals_jac2 = -Bc.*((pj./pi).^(0.5*Zc) - 1); % Values for dg/d(fgc)
    jac2 = sparse(id_jac2(:,1),id_jac2(:,2),vals_jac2,nc_g,nc_g); % Eval dg/d(fgc)
    
    id_jac3   = [(1:nc_g)' from_c(id_comp_gas)    % Indexes for dg/d(pi)
                 (1:nc_g)'  to_c(id_comp_gas)  ]; % Indexes for dg/d(pj)           
    vals_jac3 = [0.5*Zc.*Bc.*fgc.*((pj.^(0.5*Zc))./(pi.^(0.5*Zc+1)))     % Values for dg/d(pi)
                 -0.5*Zc.*Bc.*fgc.*((pj.^(0.5*Zc-1))./(pi.^(0.5*Zc)))];  % Values for dg/d(pj)
    jac3 = sparse(id_jac3(:,1),id_jac3(:,2),vals_jac3,nc_g,nn); % Eval dg/dp
    
    dg = [jac1 jac2 jac3]; % Compute Jacobian of the Lagrangian               
end
end
function d2L = wcompgas_hess(x, lambda, parcomp)
nc_g = parcomp.ncompg;  % Extract compressor parameters
from_c = parcomp.idfrom;
to_c = parcomp.idto;
nn = parcomp.other.nn;
id_comp_gas = parcomp.other.idcompg;
Bc = parcomp.Bc; Bc = Bc(id_comp_gas);
Zc = parcomp.Zc; Zc = Zc(id_comp_gas);

[~,fgc,p] = deal(x{:}); % Extract variables

pi = p(from_c(id_comp_gas));
pj = p(to_c(id_comp_gas));

id_hess1 = [from_c(id_comp_gas)  from_c(id_comp_gas)   % Indexes for d^2(L)/d^2(pi)
            to_c(id_comp_gas)    to_c(id_comp_gas)     % Indexes for d^2(L)/d^2(pj)
            from_c(id_comp_gas)  to_c(id_comp_gas)     % Indexes for d^2(L)/d(pi)d(pj)
            to_c(id_comp_gas)    from_c(id_comp_gas)]; % Indexes for d^2(L)/d(pj)d(pi)
cte1 = -0.25*lambda.*Bc.*Zc.*(Zc + 2);    
cte2 = -0.25*lambda.*Bc.*Zc.*(Zc - 2);    
cte3 = 0.25*lambda.*Bc.*(Zc.^2);
vals_hess1 = [cte1.*fgc.*((pj.^(0.5*Zc))./(pi.^(0.5*Zc+2)))     % Values for d^2(L)/d^2(pi)
              cte2.*fgc.*((pj.^(0.5*Zc-2))./(pi.^(0.5*Zc)))     % Values for d^2(L)/d^2(pj)
              cte3.*fgc.*((pj.^(0.5*Zc-1))./(pi.^(0.5*Zc+1)))   % Values for d^2(L)/d(pi)d(pj)
              cte3.*fgc.*((pj.^(0.5*Zc-1))./(pi.^(0.5*Zc+1)))]; % Values for d^2(L)/d(pj)d(pi)
hess1 = sparse(id_hess1(:,1),id_hess1(:,2),vals_hess1,nn,nn);

id_hess2 = [(1:nc_g)'   from_c(id_comp_gas)     % Indexes for d^2(L)/d(fgc)d(pi)
            (1:nc_g)'   to_c(id_comp_gas)];     % Indexes for d^2(L)/d(fgc)d(pj)
cte4 = 0.5*lambda.*Bc.*Zc;
cte5 = - cte4;          
vals_hess2 = [cte4.*((pj.^(0.5*Zc))./(pi.^(0.5*Zc+1)))    % Values for d^2(L)/d(fgc)d(pi)
              cte5.*((pj.^(0.5*Zc-1))./(pi.^(0.5*Zc)))];  % Values for d^2(L)/d(fgc)d(pj)
hess2 = sparse(id_hess2(:,1),id_hess2(:,2),vals_hess2,nc_g,nn);

L1 = spalloc(nc_g,nc_g,nc_g*nc_g);
L2 = spalloc(nc_g,nn,nc_g*nn);

d2L = [ L1    L1      L2       % Compute Hessian of the Lagrangian
        L1    L1     hess2   
        L2'  hess2'  hess1];   
end
%% -------------- [WGV] Jacobian and Hessian for equation (11) --------------
function [g, dg] = wcomppower_fcn(x, parcomp)
nc_p = parcomp.ncompp;  % Extract compressor parameters
from_c = parcomp.idfrom;
to_c = parcomp.idto;
nn = parcomp.other.nn;
id_comp_power = parcomp.other.idcompp;
Bc = parcomp.Bc; Bc = Bc(id_comp_power);
Zc = parcomp.Zc; Zc = Zc(id_comp_power);

[psi,fgc,p] = deal(x{:}); % Extract variables

if nargout < 2
    pi = p(from_c(id_comp_power));
    pj = p(to_c(id_comp_power));
    g = psi - Bc.*fgc.*((pj./pi).^(0.5*Zc) - 1); % Eval nonlinear constraint
else
    pi = p(from_c(id_comp_power));
    pj = p(to_c(id_comp_power));
    g = psi - Bc.*fgc.*((pj./pi).^(0.5*Zc) - 1); % Eval nonlinear constraint
    
    jac1 = speye(nc_p); % Eval dg/d(psi)
    
    id_jac2   = [(1:nc_p)' (1:nc_p)'];         % Indexes for dg/d(fgc)                            
    vals_jac2 = -Bc.*((pj./pi).^(0.5*Zc) - 1); % Values for dg/d(fgc)
    jac2 = sparse(id_jac2(:,1),id_jac2(:,2),vals_jac2,nc_p,nc_p); % Eval dg/d(fgc)
    
    id_jac3   = [(1:nc_p)' from_c(id_comp_power)    % Indexes for dg/d(pi)
                 (1:nc_p)'  to_c(id_comp_power)  ]; % Indexes for dg/d(pj)           
    vals_jac3 = [0.5*Zc.*Bc.*fgc.*((pj.^(0.5*Zc))./(pi.^(0.5*Zc+1)))     % Values for dg/d(pi)
                 -0.5*Zc.*Bc.*fgc.*((pj.^(0.5*Zc-1))./(pi.^(0.5*Zc)))];  % Values for dg/d(pj)
    jac3 = sparse(id_jac3(:,1),id_jac3(:,2),vals_jac3,nc_p,nn); % Eval dg/dp
    
    dg = [jac1 jac2 jac3]; % Compute Jacobian of the Lagrangian               
end
end

function d2L = wcomppower_hess(x, lambda, parcomp)
nc_p = parcomp.ncompp;  % Extract compressor paramet ers
from_c = parcomp.idfrom;
to_c = parcomp.idto;
nn = parcomp.other.nn;
id_comp_power = parcomp.other.idcompp;
Bc = parcomp.Bc; Bc = Bc(id_comp_power);
Zc = parcomp.Zc; Zc = Zc(id_comp_power);

[~,fgc,p] = deal(x{:}); % Extract variables

pi = p(from_c(id_comp_power));
pj = p(to_c(id_comp_power));

id_hess1 = [from_c(id_comp_power)  from_c(id_comp_power)   % Indexes for d^2(L)/d^2(pi)
            to_c(id_comp_power)    to_c(id_comp_power)     % Indexes for d^2(L)/d^2(pj)
            from_c(id_comp_power)  to_c(id_comp_power)     % Indexes for d^2(L)/d(pi)d(pj)
            to_c(id_comp_power)    from_c(id_comp_power)]; % Indexes for d^2(L)/d(pj)d(pi)
cte1 = -0.25*lambda.*Bc.*Zc.*(Zc + 2);    
cte2 = -0.25*lambda.*Bc.*Zc.*(Zc - 2);    
cte3 = 0.25*lambda.*Bc.*(Zc.^2);
vals_hess1 = [cte1.*fgc.*((pj.^(0.5*Zc))./(pi.^(0.5*Zc+2)))     % Values for d^2(L)/d^2(pi)
              cte2.*fgc.*((pj.^(0.5*Zc-2))./(pi.^(0.5*Zc)))     % Values for d^2(L)/d^2(pj)
              cte3.*fgc.*((pj.^(0.5*Zc-1))./(pi.^(0.5*Zc+1)))   % Values for d^2(L)/d(pi)d(pj)
              cte3.*fgc.*((pj.^(0.5*Zc-1))./(pi.^(0.5*Zc+1)))]; % Values for d^2(L)/d(pj)d(pi)
hess1 = sparse(id_hess1(:,1),id_hess1(:,2),vals_hess1,nn,nn);

id_hess2 = [(1:nc_p)'   from_c(id_comp_power)     % Indexes for d^2(L)/d(fgc)d(pi)
            (1:nc_p)'   to_c(id_comp_power)];     % Indexes for d^2(L)/d(fgc)d(pj)
cte4 = 0.5*lambda.*Bc.*Zc;
cte5 = - cte4;          
vals_hess2 = [cte4.*((pj.^(0.5*Zc))./(pi.^(0.5*Zc+1)))    % Values for d^2(L)/d(fgc)d(pi)
              cte5.*((pj.^(0.5*Zc-1))./(pi.^(0.5*Zc)))];  % Values for d^2(L)/d(fgc)d(pj)
hess2 = sparse(id_hess2(:,1),id_hess2(:,2),vals_hess2,nc_p,nn);

L1 = spalloc(nc_p,nc_p,nc_p*nc_p);
L2 = spalloc(nc_p,nn,nc_p*nn);

d2L = [ L1    L1      L2       % Compute Hessian of the Lagrangian
        L1    L1     hess2   
        L2'  hess2'  hess1];   
end
