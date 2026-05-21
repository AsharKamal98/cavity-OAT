# Half-width analysis

The half-width analysis studies how the dynamics change when the initial wavefunction is supported on a larger set of $N_J$ sectors. For a chosen half-width $d_N$, the initial state includes a symmetric set of sectors around the central value, and the number of included sectors is reported on the $x$-axis. The two-panel plot then shows:

\begin{itemize}
\item the runtime of the custom ensemble simulation as a function of the number of $N_J$ sectors,
\item the average phase-2 jump count as a function of the number of $N_J$ sectors.
\end{itemize}

The jump-count panel compares three different quantities. The first is the simulated jump count extracted directly from Monte Carlo trajectories. The second is a theory-informed jump count evaluated using the full simulated sector distribution at each saved time. The third is a simpler theoretical approximation in which the sector dependence is reduced to $\hat S_z \approx \frac{N}{2} - \hat N_J$.

## Simulated jump count

The simulated curve corresponds to the shifted jump operator used in the trajectory simulation,
\[
\hat l = \hat J_- + i\frac{\Omega}{\Gamma}.
\]

For each trajectory, all jumps whose jump times fall inside phase 2 are counted. If the phase-2 time window is $[t_1,t_2)$, then the jump count for trajectory $k$ is
\[
n_k^{(\mathrm{sim})}
=
\sum_{m} \mathbf{1}\!\left(t_{k,m}\in [t_1,t_2)\right),
\]
where $t_{k,m}$ denotes the time of the $m$-th jump in trajectory $k$.

The plotted simulated value is the sample mean over all trajectories,
\[
\bar n^{(\mathrm{sim})}

=
\frac{1}{N_{\mathrm{traj}}}
\sum_{k=1}^{N_{\mathrm{traj}}}
n_k^{(\mathrm{sim})}.
\]

## Theoretical jump count

The second curve uses the effective theoretical jump operator
\[
\hat l
\approx
\frac{\delta}{\Gamma}\tan\tilde{\theta}_J
+
\frac{2\delta\sin\tilde{\theta}_J}
{N\Gamma\cos^3\tilde{\theta}_J}
\hat S_z.
\]

This curve is not computed from realized jump events. Instead, for each trajectory one evaluates the corresponding theoretical jump rate from the full time-dependent sector distribution, integrates that rate over phase 2, and then averages over trajectories. In other words, each trajectory produces a theoretical expected jump count
\[
n_k^{(\mathrm{th})}
=
\int_{t_1}^{t_2} R_k^{(\mathrm{th})}(t)\,dt,
\]
and the plotted curve is
\[
\bar n^{(\mathrm{th})}
=
\frac{1}{N_{\mathrm{traj}}}
\sum_{k=1}^{N_{\mathrm{traj}}}
n_k^{(\mathrm{th})}.
\]

The important point is that this theoretical curve still uses the simulated trajectories as input, because the required moments of $\hat S_z$ are extracted from the simulated sector weights at each saved time.

## Theoretical approximate jump count

The third curve uses the same effective jump-operator formula,
\[
\hat l
\approx
\frac{\delta}{\Gamma}\tan\tilde{\theta}_J
+
\frac{2\delta\sin\tilde{\theta}_J}
{N\Gamma\cos^3\tilde{\theta}_J}
\hat S_z,
\qquad
\hat S_z \approx \frac{N}{2} - \hat N_J.
\]

Here the dependence on the actual simulated sector distribution is removed. Instead, the jump count is estimated directly from the half-width description alone, using the simplified relation between $\hat S_z$ and $\hat N_J$. This produces a deterministic theoretical approximation for each half-width value,
\[
n^{(\mathrm{th,approx})}(d_N),
\]
without requiring a trajectory-by-trajectory evaluation.

## Error bars

Error bars are shown for the simulated curve and for the theoretical curve, but not for the theoretical approximate curve.

### Simulated curve

For the simulated curve, the error bars are standard errors of the mean computed from the sample of trajectory jump counts:
\[
\mathrm{SEM}_{\mathrm{sim}}
=
\frac{s_{\mathrm{sim}}}{\sqrt{N_{\mathrm{traj}}}},
\]
where
\[
s_{\mathrm{sim}}
=
\sqrt{
\frac{1}{N_{\mathrm{traj}}-1}
\sum_{k=1}^{N_{\mathrm{traj}}}
\left(n_k^{(\mathrm{sim})}-\bar n^{(\mathrm{sim})}\right)^2
}.
\]

These error bars quantify the trajectory-to-trajectory fluctuations of the realized number of jumps.

### Theoretical curve

For the theoretical curve, the error bars are again standard errors of the mean, but now applied to the per-trajectory theoretical counts:
\[
\mathrm{SEM}_{\mathrm{th}}
=
\frac{s_{\mathrm{th}}}{\sqrt{N_{\mathrm{traj}}}},
\]
with
\[
s_{\mathrm{th}}
=
\sqrt{
\frac{1}{N_{\mathrm{traj}}-1}
\sum_{k=1}^{N_{\mathrm{traj}}}
\left(n_k^{(\mathrm{th})}-\bar n^{(\mathrm{th})}\right)^2
}.
\]

These quantify the spread induced by trajectory-to-trajectory differences in the simulated sector distributions that enter the theoretical rate evaluation.

### Theoretical approximate curve

No error bars are shown for the theoretical approximate curve because it is not obtained from a sample of trajectories. For each half-width, it is a single deterministic prediction from the approximation
\[
\hat S_z \approx \frac{N}{2} - \hat N_J.
\]

Since there is no ensemble of independent estimates behind this curve, there is no sample variance from which to construct a standard error bar in the same way as for the other two curves.

## Interpretation

The three curves therefore answer slightly different questions:

\begin{itemize}
\item the simulated curve measures the actual stochastic jump counts generated by the Monte Carlo wave-function dynamics,
\item the theoretical curve measures a theory-based prediction evaluated on the full simulated sector statistics,
\item the theoretical approximate curve measures a simpler closed-form prediction based only on the half-width description.
\end{itemize}

This makes the comparison useful for separating three effects: stochastic jump noise, deviations between exact simulated sector statistics and the theoretical effective description, and the additional approximation introduced by replacing $\hat S_z$ with $\frac{N}{2}-\hat N_J$.
