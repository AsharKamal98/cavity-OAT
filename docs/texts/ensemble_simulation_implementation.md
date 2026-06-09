\documentclass[11pt]{article}

\usepackage[a4paper,margin=1in]{geometry}
\usepackage{amsmath,amssymb}
\usepackage{graphicx}
\usepackage{xcolor}
\usepackage{booktabs}
\usepackage{float}
\usepackage{caption}
\usepackage{hyperref}

\newcommand{\up}{\uparrow}
\newcommand{\down}{\downarrow}
\newcommand{\ee}{e}
\newcommand{\Nj}{N_J}
\newcommand{\Ne}{N_e}
\newcommand{\GammaC}{\Gamma}
\newcommand{\OmegaDrive}{\Omega}

\setlength{\parindent}{0pt}

\title{Quantum-Trajectory Code Documentation}
\author{}
\date{}

\begin{document}
\maketitle


% \section*{\texttt{run\_trajectory\_ensemble(...)}: implementation summary}

This note explains how the custom MCWF ensemble code is executed, following the code path in \texttt{quantum\_trajectories/ensamble\_sim.py} and \texttt{quantum\_trajectories/sim.py}. The focus here is not the full physics derivation, but how the strong-symmetry structure and precomputation are used to make the simulation efficient.

\section*{1. Top-level call: \texttt{run\_trajectory\_ensemble(...)}}

The ensemble entry point is \texttt{run\_trajectory\_ensemble(...)}. Its main jobs are:

\begin{enumerate}
\item build all reusable trajectory data once with \texttt{build\_precomputed\_trajectory\_data(...)};
\item run \texttt{simulate\_single\_trajectory(...)} either serially or in parallel;
\item collect all \texttt{TrajectoryResult}s into a \texttt{TrajectoryEnsemble}.
\end{enumerate}

The important point is that the expensive objects that are identical for all trajectories are built once at ensemble level, not separately inside each trajectory.

\section*{2. Strong-symmetry reduction}

The core simplification is that the many-body state is not propagated in the full $3^N$ Hilbert space. Instead, the code uses the strong-symmetry label
\[
N_J = N_d + N_e,
\]
and represents the full wavefunction as a direct sum over populated $N_J$ sectors:
\[
|\psi(t)\rangle = \bigoplus_{N_J} c_{N_J}\,|\psi_{N_J}(t)\rangle.
\]

Inside each sector, the basis is the Dicke-like active-manifold basis labeled by $n_e=0,\dots,N_J$, so the sector dimension is only
\[
\dim(\mathcal H_{N_J}) = N_J + 1.
\]

This is the main structural reason the simulation is feasible at large $N$: we evolve a list of small sector blocks instead of one exponentially large state vector.

\section*{3. Ensemble-level precomputation: \texttt{build\_precomputed\_trajectory\_data(...)}}

Before any trajectory is run, \texttt{run\_trajectory\_ensemble(...)} calls \texttt{build\_precomputed\_trajectory\_data(...)}. This function builds all objects that depend only on
\begin{itemize}
\item \texttt{N},
\item \texttt{Gamma},
\item the protocol \texttt{phases},
\item the populated initial sectors,
\item and the fixed base timestep \texttt{dt}.
\end{itemize}

It returns a dictionary containing:

\begin{enumerate}
\item \texttt{sector\_list}: the sorted list of populated $N_J$ sectors.

\item \texttt{ops\_list}: one \texttt{SectorOperators} object per sector, created by \texttt{build\_sector\_ops(Nj)}. This packages the reduced operators such as $J_x$, $J_-$, $N_e$, and $J_+J_-$ for that sector.

\item \texttt{multiplicities}: the combinatorial multiplicity of each sector.

\item \texttt{dims}: the sector dimension $N_J+1$.

\item \texttt{phase\_jump\_operators}: for each phase and sector, the jump operator
\[
l =
\begin{cases}
J_- & \text{regular picture},\\[4pt]
J_- + i\,\Omega/\Gamma & \text{shifted picture}.
\end{cases}
\]

\item \texttt{phase\_generators}: for each phase and sector, the non-Hermitian effective generator
\[
H_{\mathrm{eff}} = H - \frac{i\Gamma}{2}\,l^\dagger l.
\]
For the regular jump operator this becomes
\[
H_{\mathrm{eff}} = \Omega J_x - \delta N_e - \frac{i\Gamma}{2}J_+J_-.
\]

\item \texttt{phase\_propagators}: for each phase and sector, the full-step propagator
\[
U_{\mathrm{eff}}(dt) = e^{-i H_{\mathrm{eff}} dt}.
\]
\end{enumerate}

This is one of the most important optimizations: as long as the protocol is piecewise constant and the base step is exactly \texttt{dt}, we do not need to recompute the matrix exponential at every timestep. \\

These operators are represented in sparse reduced-basis form, so both the basis reduction and the matrix representation are exploited.

\section*{4. Why the precompute is valid}

The precomputed data can be reused because within one phase:
\begin{itemize}
\item $\Omega$ is constant,
\item $\delta$ is constant,
\item $\Gamma$ is constant,
\item the sector basis is fixed,
\item and the nominal internal timestep is fixed to \texttt{dt}.
\end{itemize}

Therefore the same $H_{\mathrm{eff}}$ and the same full-step propagator $U_{\mathrm{eff}}(dt)$ apply every time an ordinary step of length \texttt{dt} is taken inside that phase.

\section*{5. Parallel execution: \texttt{\_init\_trajectory\_worker(...)}}

If \texttt{n\_processes > 1}, \texttt{run\_trajectory\_ensemble(...)} creates a process pool and uses \texttt{\_init\_trajectory\_worker(...)} to place the large read-only data in a process-local global dictionary \texttt{\_WORKER\_STATE}. Each worker receives the precomputed object once, then only a child \texttt{SeedSequence} is passed per trajectory. \\

The worker-side wrapper \texttt{\_simulate\_single\_trajectory\_worker(seed\_sequence)} simply reads the common data from \texttt{\_WORKER\_STATE} and calls \texttt{simulate\_single\_trajectory(...)}.

\section*{6. Single-trajectory setup: \texttt{simulate\_single\_trajectory(...)}}

Each trajectory then proceeds in \texttt{simulate\_single\_trajectory(...)}.

\subsection*{6.1 Output grid}

First, the common \texttt{t\_eval} grid is constructed by \texttt{build\_t\_eval\_from\_phases(phases, num\_snapshots)}. All trajectories in the ensemble use the same \texttt{t\_eval}, which means observables can later be averaged at identical physical times. Individual trajectories may still use additional internal timesteps to update the wavefunction, but only states saved on the common \texttt{t\_eval} grid are used when computing observables.

\subsection*{6.2 Initial state}

The initial reduced-basis state is built by \texttt{build\_initial\_sector\_state(...)} and starts from a list
\[
\{\psi_{N_J}(0)\}_{N_J \in \text{sector\_list}}.
\]

\section*{7. Main propagation loop}

The actual MCWF time evolution is then organized as
\begin{enumerate}
\item loop over protocol phases;
\item inside each phase, loop until the phase end is reached;
\item at each attempted step, decide whether the step is a full \texttt{dt} step or a shorter partial step.
\end{enumerate}

For each phase, \texttt{simulate\_single\_trajectory(...)} loads from the precompute:
\begin{itemize}
\item \texttt{jump\_operators\_list},
\item \texttt{generators\_list},
\item \texttt{full\_step\_propagators}.
\end{itemize}

\section*{8. Fast path: ordinary full \texttt{dt} steps}

If the next step length is exactly \texttt{dt}, the code uses \texttt{propagate\_blocks\_with\_propagators(...)} which applies the already computed matrices
\[
\psi_{N_J}(t+dt) = U_{\mathrm{eff},N_J}(dt)\,\psi_{N_J}(t)
\]
sector by sector. This is the main cheap path in the simulation. In particular, the code avoids recomputing
\[
e^{-iH_{\mathrm{eff}}dt}
\]
at every timestep.

\section*{9. When the precomputed path does \textbf{not} work}

The precomputed propagators are only valid for one exact step size: \texttt{dt}. Whenever the required step is shorter than \texttt{dt}, the code cannot use the precomputed matrix, because it would need
\[
U_{\mathrm{eff}}(\Delta t) = e^{-iH_{\mathrm{eff}}\Delta t},
\qquad
\Delta t \neq dt.
\]

This happens in three important places.

\subsection*{9.1 Phase boundaries}

If the next full step would cross the end of the current phase, the step is shortened to land exactly on the phase boundary. Since the duration is no longer \texttt{dt}, the precomputed full-step propagator is not valid.

\subsection*{9.2 Saved output times \texttt{t\_eval}}

If the next full step would cross a requested output time, the step is split so the state lands exactly on the next \texttt{t\_eval} point. Again this is a partial step, so the precomputed \texttt{dt} propagator cannot be used. This is the implementation cost of having a common exact output grid across trajectories.

\subsection*{9.3 Jump refinement}

If a jump occurs during a step, the jump time is located by bisection inside that step. The midpoint propagations, the propagation to the refined jump time, and any remainder after the jump are all variable-step propagations. These also cannot use the precomputed full-step matrices.

\section*{10. Variable-step path: \texttt{propagate\_blocks(...)}}

Whenever a partial step is needed, the code falls back to \texttt{propagate\_blocks(...)} which uses \texttt{expm\_multiply} sector by sector:
\[
\psi_{N_J}(t+\Delta t)
=
e^{-iH_{\mathrm{eff},N_J}\Delta t}\,\psi_{N_J}(t).
\]

This is more flexible, because it works for arbitrary $\Delta t$, but it is more expensive than applying a precomputed full-step propagator. This is why the code keeps counters for total propagation calls and for calls that could not use the precomputed path.

\section*{11. Jump detection and jump application}

After a trial propagation, the trajectory norm is checked against the MCWF threshold. If the norm stays above the threshold, the step is accepted.

If the norm falls below threshold, a jump occurred within that step. The code then:
\begin{enumerate}
\item bisects the step a fixed number of times to localize the jump time;
\item propagates to the refined jump time;
\item renormalizes the state;
\item applies the jump operator with \texttt{apply\_jump(...)};
\item records the jump time;
\item draws a new random threshold;
\item optionally propagates the remainder of the original step.
\end{enumerate}

The jump itself is applied sector by sector as
\[
\psi_{N_J} \mapsto l_{N_J}\psi_{N_J},
\]
followed by renormalization of the direct-sum state.

\section*{12. Saving snapshots}

\texttt{maybe\_save\_snapshot()} saves the state exactly when the trajectory lands on the next requested \texttt{t\_eval} point. This means: all trajectories save the same number of snapshots; snapshots occur at identical physical times across the ensemble.


\section*{13. Final output}

Finally, \texttt{simulate\_single\_trajectory(...)} returns a \texttt{TrajectoryResult} containing:
\begin{itemize}
\item the saved snapshots;
\item final sector blocks;
\item the jump times and jump count;
\item sector metadata;
\item the common \texttt{t\_eval} grid;
\item diagnostic counters for total propagation calls and calls that did not use precompute.
\end{itemize}

\texttt{run\_trajectory\_ensemble(...)} collects these results, prints the average step diagnostics per trajectory, and returns a \texttt{TrajectoryEnsemble}.

\section*{14. Summary of the main efficiency ideas}

The implementation is efficient mainly because of four choices:
\begin{enumerate}
\item \textbf{Strong-symmetry block decomposition.} The simulation evolves only the populated $N_J$ sectors, each with dimension $N_J+1$ (for homogeneous couplings), instead of the full $3^N$ Hilbert space.

\item \textbf{Precomputed sector operators.} Reduced operators such as $J_x$, $J_-$, $N_e$, and $J_+J_-$ are built once per sector.

\item \textbf{Precomputed phase-dependent generators and propagators.} For each phase and sector, the code builds $H_{\mathrm{eff}}$ and the full-step propagator $e^{-iH_{\mathrm{eff}}dt}$ once, and reuses them on ordinary steps.  

\item \textbf{Shared read-only ensemble data.} In multiprocessing mode, the expensive precomputed data are loaded once per worker and reused across many trajectories.
\end{enumerate}

The main situations where this efficiency breaks down are exactly the situations where the step length is no longer the fixed \texttt{dt}: phase boundaries, \texttt{t\_eval} boundary splits, and jump-time refinement.

\section*{15. What changes for inhomogeneous couplings}

The main implementation change is that the sector label is no longer a single integer \(N_J\). In the current two-group implementation, each sector is labeled by
\[
(N_{J,1},N_{J,2}),
\]
and the basis inside each block becomes
\[
|n_{e,1},n_{e,2}\rangle,
\qquad
n_{e,1}=0,\dots,N_{J,1},
\qquad
n_{e,2}=0,\dots,N_{J,2}.
\]
So the sector dimension changes from \(N_J+1\) in the homogeneous code to
\[
(N_{J,1}+1)(N_{J,2}+1).
\]

As a result, several stored objects change shape: \texttt{sector\_coeffs}, \texttt{sector\_list}, \texttt{sector\_blocks}, \texttt{sector\_dimensions}, and the saved snapshots all use tuple keys instead of scalar \(N_J\) keys.

The operator construction also changes. Instead of building one reduced operator bundle for a single Dicke sector, the code builds one reduced operator bundle for each tuple sector \((N_{J,1},N_{J,2})\). The one-group operators are constructed for each group and combined once through sparse Kronecker products. The weighted jump and drive operators stored in the precomputed data are
\[
A = \omega_1 J_{1,-} + \omega_2 J_{2,-},
\qquad
\omega_2 =
\frac{N-N_1\omega_1}{N_2}.
\]
Here \(\omega_2\) is fixed once from the physical group sizes, not recomputed separately for each \((N_{J,1},N_{J,2})\) sector.
\[
\Omega(\omega_1 J_{1,x} + \omega_2 J_{2,x}).
\]
In the shifted picture, the stored jump operator becomes
\[
l = A + i\frac{\Omega}{\Gamma}.
\]

The overall simulation strategy stays the same. The code still:
\begin{itemize}
\item builds the sector operators once per populated sector;
\item precomputes phase-dependent jump operators, effective generators, and fixed-\texttt{dt} propagators;
\item reuses those precomputed objects across all trajectories;
\item falls back to variable-step propagation only for partial steps.
\end{itemize}

So the main change is not the time-stepping logic, but the block structure that is propagated and precomputed: homogeneous couplings use one Dicke block per \(N_J\), while inhomogeneous couplings use one product-Dicke block per \((N_{J,1},N_{J,2})\).



\end{document}
