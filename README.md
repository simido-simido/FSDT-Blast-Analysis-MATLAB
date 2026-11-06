# FSDT blast FEM adapted from Ferreira's first-edition MATLAB code


It predicts **linear elastic, small-deflection** response to a **prescribed
pressure history**, using five degrees of freedom per node and first-order
shear-deformation theory (FSDT). It includes general ply data, modal analysis,
full-system Newmark integration, probe stress recovery, and numerical checks.

The Friedlander/cubic pressure definitions follow the same mathematical pulse
forms used in the supplied DYNAblast 2.0 code. This is a separate FEM model;
it does not reproduce DYNAblast's GUI, all calibration options, or its nonlinear
plate formulations. There is no dependency on the DYNAblast installation.

## Run in MATLAB

Extract the ZIP, open the `FSDT_Blast_Ferreira_v1` folder in MATLAB, then run:

```matlab
result = run_fsdt_blast;
```

The inputs are in `fsdt_config.m`. No Symbolic Math Toolbox or FEM toolbox is
needed. The original encrypted `.p` plotting/meshing files are not needed.
Modern base MATLAB is the intended environment; the package was executed in
GNU Octave 8.4.0 during preparation, not in proprietary MATLAB.

Other entry points:

```matlab
report = run_validation;       % numerical verification; errors on a failed check
summary = run_study;           % 8 cases: layup x support x negative-phase option
checks = run_convergence;      % actual blast mesh/time refinement tables
```

For a custom case, retain all fields supplied by the default configuration:

```matlab
cfg = fsdt_config;
cfg.name = 'my_laminate';
cfg.plies(1).angle = 45;       % degrees; bottom ply
cfg.plies(2).angle = -45;
cfg.plies(3).angle = -45;
cfg.plies(4).angle = 45;
cfg.pulse.pmax = 1000;         % Pa
cfg.pulse.includeNegative = false;
cfg.output.folder = fullfile(pwd,'results',cfg.name);
result = run_fsdt_blast(cfg);
```

Set `cfg.output.makePlots=false` for a machine without graphics. Set
`saveFullHistory=true` to retain every nodal displacement, free velocity, and
free acceleration history in the MAT result. This increases memory and disk use.

## Demonstration inputs and meaning

- Plate: 0.50 m x 0.40 m, four 1 mm plies, bottom-to-top [0/90/90/0].
- Orthotropic ply: E1=135 GPa, E2=10 GPa, nu12=0.30, G12=G13=5 GPa,
  G23=3.8 GPa, rho=1600 kg/m3.
- Q4 structured mesh: 12 x 10 elements, 143 nodes, 715 total DOFs.
- Simply supported bending conditions with explicitly defined Ferreira
  membrane constraints; see below.
- Uniform pulse: pmax=1000 Pa, td=1 ms, decay=1.5, followed by cubic
  underpressure of magnitude 100 Pa and duration 4 ms.
- tEnd=40 ms; dt=20 microseconds; zero damping; initially at rest.

These are demonstrative engineering inputs, not a tested commercial material
or an experimentally calibrated explosive event. All quantities are in
consistent **m, kg, s, Pa** units. Positive gauge pressure acts along negative
z by default; `cfg.load.direction` controls that sign. Underpressure is negative
gauge pressure and reverses the applied force.

## What comes out

Results are written under `results/<case>/`:

| File | Contents |
|---|---|
| response.csv | Time, pressure, probe displacement/velocity/acceleration, w/h, centre displacement |
| frequencies.csv | Angular frequency, frequency in Hz, transverse mass fraction, eigen residual |
| probe_ply_stress.csv | Global and local in-plane stress on both faces of every ply, at the probe |
| energy.csv | Kinetic energy, strain energy, external work, damping loss, balance error |
| mesh_nodes.csv / mesh_elements.csv | Explicit mesh and connectivity |
| peak_probe_time_field.csv | All nodal fields at the time of maximum absolute probe displacement |
| result.mat | Configuration, laminate matrices, modes, histories, peak field, 

CSV stresses are in **Pa**, while plot stresses are explicitly converted to
**MPa**. Stress is evaluated at the configured probe point; it is not a search
for the maximum stress over the whole plate. At an element boundary, the stress
is the one-sided value from the element selected by `fsdt_probe`.

## Supports and DOF convention

The global vector follows Ferreira's field blocks:

`[w(all nodes); theta_x(all); theta_y(all); u(all); v(all)]`.

At an element, the corresponding field blocks each contain its four nodes.
The rotations are defined by u=u0+z*theta_x and v=v0+z*theta_y; use these
definitions rather than guessing a shell rotational sign from a GUI label.

- `transverse='ssss'`: w=0 at every edge, theta_x=0 on y edges,
  theta_y=0 on x edges, matching the source example.
- `transverse='cccc'`: w=theta_x=theta_y=0 at every edge.
- `membrane='ferreira'`: u=0 at x=0 and v=0 at y=0 for ssss; u=v=0
  on every edge for cccc, matching the original routine.
- `membrane='immovable'`: u=v=0 on every edge, independent of bending support.
- `membrane='minimal'`: u=v=0 at the bottom-left corner, v=0 at the
  bottom-right corner; removes in-plane rigid motions with otherwise free edges.

The study uses `membrane='immovable'` for both support types so that its
in-plane restraint does not change along with the bending support. Labels such
as "simply supported" alone are insufficient when matching a coupled laminate.

## Pulse choices

- `friedlander_cubic`: Friedlander positive phase and optional finite cubic
  negative phase.
- `friedlander`: positive phase only.
- `extended`: Friedlander expression continued beyond td into its negative
  exponential tail when `includeNegative=true`; pmin/tm do not define that tail.
- `triangular` and `rectangular`: finite positive pulses; the negative option
  has no effect for these types.

Arrival time is t=0; this version does not apply a travel delay or infer a pulse
from charge mass and distance. A rectangular step is sampled with p=0 at td;
time-step refinement is necessary at discontinuities.

`cfg.load.patch=[]` applies pressure to the full plate. A patch is
`[xmin xmax ymin ymax]` in metres and must coincide with element boundaries.
For the default mesh, `[0.125 0.375 0.08 0.32]` is an aligned example.
Cut-element pressure integration is intentionally rejected rather than
approximated by an undocumented quadrature rule.

## Verification and practical limits

`run_validation` verifies source assembly equivalence, physical mass and load
resultants, rigid motion, section inertia, pulse integrals, time integration,
and an independent analytical FSDT plate solution. Its data and definitions are
in `verification_results/` and `VERIFICATION_REPORT.md`.

The independent plate transient benchmark is a **spatial sine load with a
temporal step**, not an experimental blast history. The demonstration pulse
has separate mesh and time sensitivity tables. A default input file is not a
claim that every requested output has converged, especially ply stresses.

This model does not include geometric membrane stiffening, plasticity,
progressive damage, delamination, fracture, fluid-structure interaction, or
explosive chemistry. In particular, large w/h should trigger assessment of
the small-deflection assumption, not a claim of failure prediction. The model
uses constant user-selected shear correction factors; it does not reproduce
Ferreira's example-specific energy-equivalent correction calculation.

The Q4 element retains selective reduced integration. Thin/slender meshes,
distorted grids, arbitrary laminates, and local stress fields require their own
locking, mechanism, and convergence assessments. No universal accuracy is
claimed from the included cases.
