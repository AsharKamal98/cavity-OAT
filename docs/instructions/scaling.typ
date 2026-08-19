#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 0.7em)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 17pt, weight: "bold")[SURFcua MCWF Runtime and Memory Scaling]
]

= Purpose and Scope

Use this file to estimate SURFcua runtime for balanced two-group MCWF runs
`Ni=(m, m)` from `m` and `dN`, where $N=2m$. The scaling assumes $Gamma=1$,
$Omega,delta prop N Gamma$, and $"dt" prop 1/(N Gamma)$. Keep the protocol
durations, solver configuration, and multiprocessing efficiency comparable.

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
    columns: (auto, auto, auto, auto),
    align: (center, center, center, center),
    inset: 6pt,
    stroke: 0.5pt,
    [*`dN`*], [*Stored sectors*], [*$F_d$ runtime*], [*Memory factor*],
    [`0`], [11], [$1.0$], [$1.0$],
    [`1`], [31], [$2.9$], [$3.0$],
    [`2`], [49], [$4.7$], [$4.9$],
    [`3`], [65], [$6.5$], [$6.7$],
    [`4`], [79], [$8.2$], [$8.5$],
  )
]

For balanced groups, the number of stored sectors is

$
S(m,d) = (2d+1)(m+1)-d(d+1), quad d = "dN".
$

== Vary `Ni=(m, m)` at `dN=0`

#align(center)[
  #table(
    columns: (auto, auto, auto, auto),
    align: (center, center, center, center),
    inset: 6pt,
    stroke: 0.5pt,
    [*`Ni`*], [*Stored sectors*], [*$F_N$ runtime*], [*Memory factor*],
    [`(10, 10)`], [11], [$1.0$], [$1.0$],
    [`(20, 20)`], [21], [$13$], [$6.2$],
    [`(30, 30)`], [31], [$110$], [$19$],
    [`(40, 40)`], [41], [$540$], [$43$],
    [`(50, 50)`], [51], [$1900$], [$82$],
  )
]

At `dN=0`, there are $m+1$ stored sectors and

$
sum_(j=0)^m (j+1)(m-j+1) = binom(m+3, 3)
$

complex snapshot amplitudes. Memory is proportional to this sum when `ntraj`
and `num_snapshots` are fixed.

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

Local serial runs used two trajectories, 500 snapshots, and
$T_1=T_2=T_3=50$ with the baseline parameter scaling.

#align(center)[
  #text(size: 9pt)[
    #table(
      columns: (auto, auto, auto, auto, auto),
      align: (center, center, center, center, center),
      inset: 5pt,
      stroke: 0.5pt,
      [*Case*],
      [*Runtime / trajectory*],
      [*Runtime factor (measured / predicted)*],
      [*Pickle size*],
      [*Memory factor (measured / predicted)*],
      [`Ni=(10,10)`, `dN=0`], [1.43 s], [1.00 / 1.00], [5.03 MB], [1.00 / 1.00],
      [`Ni=(20,20)`, `dN=0`], [7.49 s], [5.24 / 13], [29.20 MB], [5.81 / 6.2],
      [`Ni=(10,10)`, `dN=2`], [6.55 s], [4.58 / 4.7], [24.15 MB], [4.80 / 4.9],
    )
  ]
]

The `dN` runtime factor agreed within about 3%, and both pickle-memory factors
agreed within about 7%. The `Ni` runtime factor overestimated the measured
increase by about 2.5 times. Treat the $F_N$ runtime table and combined runtime
formula as provisional until an additional SURFcua `Ni`-scaling run is
available; the `dN` runtime and memory factors are suitable planning estimates.
