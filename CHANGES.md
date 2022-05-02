Change history for MATPOWER-Natural Gas (MPNG)
==============================================

#### 28/04/2022 (Changes included in released version 0.99beta)

- Addition of unit commitment as an input parameter. Extra field `connect.power.UC`

- Addition of inter-temporal power ramp constraints. Extra field `connect.power.ramp_time` for setting generation ramp-times.

- Improved nargout checking in `wey_approx.m` for passing only needed data based on Jacobian and Hessian evaluation calls from `mpng.m`

- Redefinition of the default region of approximation of Weymouth's equation (Pi*), according to a percentage of the pipeline's maximum flow.

- Addition of a new optional field `mgc.wey_percent` for allowing the user to specify a different approximation region for Weymouth's equation (Pi*).

- Temporal workaround for allowing compressor ratios to be lower or equal to 1.0 in order to model turboexpanders (i.e compressors with negative power consumption, that is, positive generation). *Thanks to a question from Roman Korab*.

- Minor fix to ensure that actual nodal gas demand in `mgc.node.dem` matches total gas demand in `mgc.node.info(:,GD)`.

- Fixes in per unit conversion of some quantities.

- Release of `Example3.m`. It shows MPNG's performance taking turboexpanders into account.

- Release of `Example4.m` and `connection_pg_case118.m`. This example simulates the 118-bus & 48-node power&gas system over a 24-hour horizon and tests the new MPNG's features regarding the assessment of given unit commitment schedules and start-up/shut-down generation ramps.

- A new option has been enabled in all `Examples` for allowing passing the initial point contained in power&gas cases to the solver.
