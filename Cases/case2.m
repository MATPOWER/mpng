function mpc = case2
%CASE2   Power flow data for 2 bus, 1 generator case.
%   This case is used to run an independent Optimal Natural Gas Flow in MPNG

%   MPNG
%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
    1	3	0	0	0	0	1	1	0	345	1	1.1	0.9;
    2	1	5	1	0	0	1	1	0	345	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	5	0   300     -300	1.0	100     1	15	0   0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	2	0	0.05	0	250	250	250	0	0	1	-360	360
    ];

%%-----  OPF Data  -----%%
%% generator cost data

mpc.gencost = [
	2	1500	0	3	0	5	0;
];

%% gen fuel
mpc.genfuel = {'hydro'};