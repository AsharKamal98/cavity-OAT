\section{Inhomogeneous Coupling (two-group split)}
We split the atoms into two groups according to the drive they feel. The Hamiltonian changes as
\begin{equation}
\hat{H}_{\delta}
=
\Omega \hat{J}_x
-
\delta \hat{N}_e
\quad
\rightarrow
\quad
\hat{H}_{\delta}
=
\Omega
(
\omega_1 \hat{J}_{1x}
+
\omega_2 \hat{J}_{2x}
)
-
\delta
(
\hat{N}_{e1}
+
\hat{N}_{e2}
),
\end{equation}
with
\begin{equation}
\hat{J}_{1x}
=
\frac{1}{2}
(
\hat{J}_{1+}
+
\hat{J}_{1-}
),
\qquad
\hat{J}_{2x}
=
\frac{1}{2}
(
\hat{J}_{2+}
+
\hat{J}_{2-}
).
\end{equation}

The natural Dicke/occupation basis is now described by the number of excitations in each subgroup,
\begin{equation}
|n_e\rangle
\quad
\rightarrow
\quad
|n_{e1},n_{e2}\rangle
=
|n_{e1}\rangle_1
\otimes
|n_{e2}\rangle_2.
\end{equation}

Previously, the wavefunction could be thought of as split into different sectors labelled by 
\begin{equation}
    N_J,
    \quad 
    0 \leq N_J \leq N,
\end{equation}
with each $N_J$ sector being decoupled from the others. Now, the $N_J$ sectors are further resolved into sectors labelled by
\begin{equation}
(N_{J1}, N_{J2}),
\quad
0 \leq N_{J1} \leq N_1,
\quad
0 \leq N_{J2} \leq N_2,
\end{equation}
where $N_1$ and $N_2$ denote the number of atoms with respective drives $\omega_1 \Omega$ and $\omega_2 \Omega$, with $N_1 + N_2 = N$. Here, $N_{J1}$ and $N_{J2}$ count how many atoms in each group belong to the active $\{|\downarrow\rangle, |e\rangle\}$ manifold.
\\

Analogously to before, different $(N_{J1}, N_{J2})$ sectors do not mix. However, within a fixed $(N_{J1}, N_{J2})$ sector, the excitation states
\begin{equation}
(n_{e1}, n_{e2})
\end{equation}
do mix under the Hamiltonian and collective dissipation, with
\begin{equation}
0 \leq n_{e1} \leq N_{J1},
\qquad
0 \leq n_{e2} \leq N_{J2}.
\end{equation}

Below, we highlight the steady-state phase equation. Its detailed mean-field derivation can be found in \texttt{mean\_field\_equations\_inh\_couplings}.
\begin{equation}
\boxed{
\begin{aligned}
\frac{\Omega\omega_a}{2}
e^{-i\phi_a}
\sin\theta_a
=
&\frac{\delta}{2}
\sin \theta_a \tan \theta_a
-
i\frac{\Gamma \omega_a}{4}
e^{-i\phi_a}
\sin\theta_a
\left[
\omega_1 N_{J,1} e^{i\phi_1} \sin\theta_1
+
\omega_2 N_{J,2}e^{i\phi_2} \sin\theta_2
\right].
\end{aligned}
}
\end{equation}
This is the steady-state phase equation for group \(a=1,2\).
\\

As a sanity check, in the homogeneous limit,
\begin{equation}
\omega_1=\omega_2=1,
\qquad
\theta_1=\theta_2=\theta_J,
\qquad
\phi_1=\phi_2=\phi_J,
\end{equation}
the group equation reduces to
\begin{equation}
\begin{aligned}
\frac{\Omega}{2}
e^{-i\phi_J}
\sin\theta_J
=
&\frac{\delta}{2}
\sin \theta_J \tan \theta_J
-
i\frac{\Gamma N_J}{4}
\sin^2\theta_J.
\end{aligned}
\end{equation}
This matches the homogeneous steady-state phase equation derived in the original paper.

\subsection{Simulation check: homogeneous limit}

As a first numerical check, we compare the original homogeneous implementation with the inhomogeneous implementation in the homogeneous limit
\begin{equation}
\omega_1=\omega_2=1.
\end{equation}
Code-wise, the inhomogeneous run still uses two physical subgroups and the product Dicke basis
\begin{equation}
|n_{e1},n_{e2}\rangle
=
|n_{e1}\rangle_1
\otimes
|n_{e2}\rangle_2,
\end{equation}
but both groups feel the same drive. The parameters used are
\begin{equation}
N=50,
\qquad
N_1=N_2=25,
\qquad
dN=2,
\qquad
N_{\rm traj}=100,
\qquad
\Gamma=1.
\end{equation}

Figure~\ref{fig:inh_homogeneous_limit_observables} compares the main observables from the homogeneous and inhomogeneous codes. The agreement of the polar and azimuthal angles, in particular \(\theta(t)\) and \(\phi(t)\), shows that resolving the Hilbert space into two subgroup Dicke bases does not change the dynamics when \(\omega_1=\omega_2\).

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{output/inhomogeneous_nj1_zero_vs_homogeneous.png}
\caption{Homogeneous-limit validation of the inhomogeneous code. The inhomogeneous simulation still uses two subgroup Dicke bases, but with \(\omega_1=\omega_2=1\). The relevant observables agree with the original homogeneous implementation, including the angular dynamics \(\theta(t)\) and \(\phi(t)\).}
\label{fig:inh_homogeneous_limit_observables}
\end{figure}

Figure~\ref{fig:inh_homogeneous_limit_residuals} shows the two complex residuals \(R_1\) and \(R_2\) associated with the group steady-state equations above, plotted separately as \(|{\rm Re}\,R_a|\) and \(|{\rm Im}\,R_a|\). The residuals are nonzero at the beginning of a phase, but relax close to zero as the state approaches the corresponding phase steady state.

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{output/inhomogeneous_mfe_residuals.png}
\caption{Mean-field steady-state residuals for the inhomogeneous run in the homogeneous limit. The real and imaginary parts of \(R_1\) and \(R_2\) are plotted separately. The residuals relax after each phase change, indicating convergence toward the steady state of the current phase.}
\label{fig:inh_homogeneous_limit_residuals}
\end{figure}

The residual magnitudes at the ends of the three phases were
\begin{equation}
\begin{array}{c|c}
\text{phase} & |R_1|+|R_2| \\
\hline
1 & 9.898346\times 10^{-2} \\
2 & 4.838606\times 10^{-3} \\
3 & 0
\end{array}
\end{equation}
The remaining phase-1 residual is larger than in phases 2 and 3, but the overall trend confirms that the inhomogeneous implementation reaches the appropriate steady states in this validation case.

The cost of the inhomogeneous representation is substantially higher even when the results match the homogeneous code. This is expected: instead of one Dicke block per total \(N_J\), the inhomogeneous code propagates many \((N_{J1},N_{J2})\) blocks, each with dimension
\begin{equation}
(N_{J1}+1)(N_{J2}+1).
\end{equation}
For this run, the measured runtimes were
\begin{equation}
t_{\rm hom}=4.977\,{\rm s},
\qquad
t_{\rm inh}=118.926\,{\rm s},
\qquad
\frac{t_{\rm inh}}{t_{\rm hom}}=23.894.
\end{equation}
Thus, the homogeneous-limit test validates the implementation, while also showing the expected runtime penalty from resolving the two subgroup sectors explicitly.

\subsection{Simulation with unequal couplings}

We next run a genuinely inhomogeneous case with
\begin{equation}
\omega_1=0.75.
\end{equation}
The second coupling is chosen so that the atom-number weighted average coupling is unchanged relative to the homogeneous run,
\begin{equation}
N_1\omega_1+N_2\omega_2=N,
\qquad
\omega_2=\frac{N-N_1\omega_1}{N_2}.
\end{equation}
For \(N_1=N_2=25\), this gives
\begin{equation}
\omega_2=1.25.
\end{equation}
This normalization keeps the average coupling equal to one, so the inhomogeneous run can be compared directly to the homogeneous reference.

Figure~\ref{fig:inh_omega1_small_observables} compares the homogeneous and inhomogeneous dynamics for this unequal-coupling case. The azimuthal angle \(\phi(t)\) remains essentially the same, while the polar angle \(\theta(t)\) differs: the inhomogeneous curve lies below the homogeneous one.

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{output/inhomogeneous_nj1_zero_vs_homogeneous_omega1_small.png}
\caption{Comparison of homogeneous and inhomogeneous dynamics for \(\omega_1=0.75\), with \(\omega_2=1.25\) chosen to preserve the average coupling. The azimuthal angle agrees well, while the inhomogeneous polar angle is shifted below the homogeneous result.}
\label{fig:inh_omega1_small_observables}
\end{figure}

Figure~\ref{fig:inh_omega1_small_residuals} shows the residuals of the two group steady-state equations for the same run. \textbf{Placeholder: discuss what the residual plot shows here.}

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{output/inhomogeneous_mfe_residuals_omega1_small.png}
\caption{Mean-field steady-state residuals for the unequal-coupling run with \(\omega_1=0.75\). The real and imaginary parts of \(R_1\) and \(R_2\) are plotted separately.}
\label{fig:inh_omega1_small_residuals}
\end{figure}

The residual magnitudes at the ends of the three phases were
\begin{equation}
\begin{array}{c|c}
\text{phase} & |R_1|+|R_2| \\
\hline
1 & 9.366626\times 10^{-2} \\
2 & 1.675788\times 10^{-2} \\
3 & 4.184092\times 10^{-4}
\end{array}
\end{equation}

To better understand the unequal-coupling dynamics, we also inspect the angles of the two subgroup collective Bloch vectors separately, shown in Fig.~\ref{fig:inh_omega1_small_group_angles}. The two groups develop different polar angles, \(\theta_1(t)\neq\theta_2(t)\). The smaller-coupling group has the smaller \(\theta\), consistent with weaker excitation. The subgroup angles show oscillatory behavior in both \(\theta_a\) and \(\phi_a\), but these oscillations largely cancel in the combined observables, producing approximately constant average \(\theta(t)\) and \(\phi(t)\).

\begin{figure}[H]
\centering
\includegraphics[width=\textwidth]{output/inhomogeneous_group_angles_omega1_small.png}
\caption{Group-resolved active-manifold angles for the unequal-coupling run. The two subgroups have different polar angles because they experience different effective drives. The subgroup oscillations largely cancel in the averaged collective angles.}
\label{fig:inh_omega1_small_group_angles}
\end{figure}
