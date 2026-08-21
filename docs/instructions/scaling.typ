#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 0.7em)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 17pt, weight: "bold")[SURFcua MCWF Runtime and Moment-Storage Scaling]
]

= Purpose and Scope

Use this file to estimate SURFcua runtime for balanced two-group MCWF runs
`Ni=(m, m)` from `m` and `dN`, where $N=2m$. The scaling assumes $Gamma=1$,
$Omega,delta prop N Gamma$, and $"dt" prop 1/(N Gamma)$. Keep the protocol
durations, solver configuration, and multiprocessing efficiency comparable.
The stored-output estimate uses the trajectory-resolved generator-moment
contract described in
`docs/instructions/solvers/mcwf/ensemble_simulation.typ`.

= SURFcua Baseline

The reference run used

$
Omega_0 = 0.1 N Gamma,
quad delta_0 = 0.02 N Gamma,
quad Gamma=1.
$

No explicit `dt` was supplied. With `drive_factor=0.01` and
`decay_factor=0.1`, `mcwf_dt_from_scales(...)` gives

$
"dt" = min(0.1/(N Gamma), 0.01/abs(Omega_0), 0.01/abs(delta_0))
     = 0.1/(N Gamma).
$

#align(center)[
  #table(
    columns: (auto, auto, auto, auto, auto, auto),
    align: (center, center, center, center, center, center),
    inset: 6pt,
    stroke: 0.5pt,
    [*`Ni`*], [*`dN`*], [*Elapsed*], [*CPUs*], [*Trajectories*], [*Baseline*],
    [`(10, 10)`], [`0`], [`03:21` = 201 s], [64], [200], [64.32 CPU-s/traj.],
  )
]

Here CPU-seconds per trajectory is the requested seconds-per-trajectory-per-CPU
metric:

$
b_0 = (t_"wall" n_"CPU") / n_"traj"
    = (201 " s" times 64) / 200
    = 64.32 " CPU-s/trajectory".
$

= Scaling Factors

== Vary `dN` at `Ni=(10, 10)`

#align(center)[
  #table(
    columns: (auto, auto, auto),
    align: (center, center, center),
    inset: 6pt,
    stroke: 0.5pt,
    [*`dN`*], [*Retained sectors*], [*$F_d$ runtime*],
    [`0`], [11], [$1.0$],
    [`1`], [31], [$2.9$],
    [`2`], [49], [$4.7$],
    [`3`], [65], [$6.5$],
    [`4`], [79], [$8.2$],
  )
]

For balanced groups, the number of represented sectors is

$
S(m,d) = (2d+1)(m+1)-d(d+1), quad d = "dN".
$

== Vary `Ni=(m, m)` at `dN=0`

#align(center)[
  #table(
    columns: (auto, auto, auto),
    align: (center, center, center),
    inset: 6pt,
    stroke: 0.5pt,
    [*`Ni`*], [*Retained sectors*], [*$F_N$ runtime*],
    [`(10, 10)`], [11], [$1.0$],
    [`(20, 20)`], [21], [$13$],
    [`(30, 30)`], [31], [$110$],
    [`(40, 40)`], [41], [$540$],
    [`(50, 50)`], [51], [$1900$],
  )
]

At `dN=0`, there are $m+1$ represented sectors and

$
sum_(j=0)^m (j+1)(m-j+1) = binom(m+3, 3)
$

sector-block amplitudes in the propagated state. This dimension affects
propagation cost and the working memory used while extracting moments.

= Stored Moment Size

For $G$ physical groups, there are

$
P=8G
$

collective generators and

$
K=frac(P(P+1), 2)
$

packed symmetric second moments. Each trajectory and saved time therefore
stores

$
F=P+K=32G^2+12G
$

real floating-point values. With `float64` arrays, the first- and second-moment
payload is

$
B_"moments" = 8 R n_t (32G^2+12G) " bytes",
$

where $R$ is the trajectory count and $n_t$ is the number of saved times. For
two groups, $F=152$, or 1216 bytes per trajectory and saved time. The stored
array size depends on $R$, $n_t$, and $G$; the retained sector count and sector
dimensions instead affect simulation runtime and transient working memory.

= Future Runtime Estimate

For inputs `Ni=(m, m)` and `dN=d`, estimate the required CPU-seconds per
trajectory as

$
b(m,d) = 64.32 F_N(m) F_d(d)
quad "CPU-s/trajectory".
$

Use the factors from the two tables. This product treats the `Ni` and `dN`
effects as approximately separable. For $n_"traj"$ trajectories on
$n_"CPU"$ CPUs, convert it to wall time with

$
t_"wall" = (b(m,d) n_"traj") / n_"CPU".
$

The factors are planning estimates. Recalibrate them when more SURFcua timing
data or materially different solver settings become available.

= Small-Run Validation

Local serial runs used two trajectories, 500 saved times, and
$T_1=T_2=T_3=50$ with the baseline parameter scaling.

#align(center)[
  #text(size: 9pt)[
    #table(
      columns: (auto, auto, auto),
      align: (center, center, center),
      inset: 5pt,
      stroke: 0.5pt,
      [*Case*],
      [*Runtime / trajectory*],
      [*Runtime factor (measured / predicted)*],
      [`Ni=(10,10)`, `dN=0`], [1.43 s], [1.00 / 1.00],
      [`Ni=(20,20)`, `dN=0`], [7.49 s], [5.24 / 13],
      [`Ni=(10,10)`, `dN=2`], [6.55 s], [4.58 / 4.7],
    )
  ]
]

The `dN` runtime factor agreed within about 3%. The `Ni` runtime factor
overestimated the measured increase by about 2.5 times. These timings predate a
dedicated calibration of the current per-saved-time generator-moment extraction,
whose cost also grows with the number of generators and generator pairs. Treat
both runtime-factor tables and the combined runtime formula as provisional
until current moment-producing runs have been benchmarked.
