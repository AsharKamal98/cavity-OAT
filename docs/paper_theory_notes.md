# Paper Theory Notes

Source paper: **Engineering One Axis Twisting via a Dissipative Berry Phase Using Strong Symmetries**.

## Conceptual goal

The paper engineers one-axis twisting from a dissipative Berry phase. The system is driven and dissipative, but the useful squeezing is generated coherently through phase accumulation between strong-symmetry sectors.

The central idea is:

\[
\text{different } N_J \text{ sectors}
\quad\Rightarrow\quad
\text{different Berry phases}
\quad\Rightarrow\quad
\text{effective OAT}.
\]

The strong symmetry protects the inter-sector coherences needed for the Berry phases to matter.

## Strong symmetry

The conserved quantity is

\[
N_J=N_\downarrow+N_e.
\]

The state \(|\uparrow\rangle\) is not coupled to the \(|\downarrow\rangle,|e\rangle\) manifold, so the number of atoms in the active manifold is conserved.

A strong symmetry means the jump operators themselves respect the symmetry. Therefore each quantum trajectory remains in its initial symmetry sector.

This differs from a weak symmetry, where conservation may only hold after averaging over trajectories.

## Why strong symmetry matters

The initial state is a superposition over \(N_J\) sectors. For example, an equal superposition of \(|\uparrow\rangle\) and \(|\downarrow\rangle\) for each atom gives a binomial distribution over \(N_J\).

The protocol needs the relative phase between sectors to remain coherent:

\[
\sum_{N_J} c_{N_J}|N_J\rangle
\rightarrow
\sum_{N_J} c_{N_J}e^{i\phi_B(N_J)}|N_J\rangle.
\]

If the environment learned \(N_J\), the off-diagonal coherences between sectors would decay. Then the Berry phase would not generate useful OAT.

In the ideal spin-polarized phase at \(\delta=0\), the sectors are indistinguishable to the environment.

## Dark state picture

In the usual Hamiltonian picture,

\[
H=\Omega J_x,
\qquad
l=J_-.
\]

The steady state is reached by a balance between coherent drive and dissipation.

In the shifted-jump picture,

\[
H=0,
\qquad
l=J_-+i\frac{\Omega}{\Gamma}
\]

at \(\delta=0\). In this picture, the steady state is dark:

\[
l|\psi_{\rm ss}\rangle=0.
\]

This means the net output field vanishes. It does not mean \(J_-|\psi\rangle=0\). Instead,

\[
J_-|\psi_{\rm ss}\rangle
=
-i\frac{\Omega}{\Gamma}|\psi_{\rm ss}\rangle.
\]

The atomic radiation cancels the coherent drive contribution in the output channel.

## Physical meaning of no \(N_J\) leakage

The environment monitors the physical output field, represented by the shifted jump operator. If the shifted output field is zero or independent of \(N_J\), the environment cannot learn which \(N_J\) sector the system occupies.

Thus, coherence between sectors is preserved.

This is why the shifted-jump picture is useful: it makes environmental distinguishability explicit.

## Detuning and Berry phase

During phase 2,

\[
\delta\neq0.
\]

The steady-state Bloch vector in the \(|\downarrow\rangle,|e\rangle\) manifold changes phase. This generates a Berry phase.

The Berry phase depends on \(N_J\), because the steady-state angle depends on the collective spin length \(N_J/2\). Therefore different sectors acquire different phases.

Expanding this dependence gives

\[
H_{\rm eff}\sim S_z^2,
\]

with

\[
S_z\approx \frac{N}{2}-N_J.
\]

This is the emergent OAT.

## Collective dephasing

The same detuning that produces Berry phase also makes the shifted output field nonzero and \(N_J\)-dependent.

The effective shifted jump operator contains

\[
l_{\rm eff}\propto S_z.
\]

This means the environment weakly measures \(S_z\), equivalently \(N_J\). This causes collective dephasing.

For an effective jump

\[
L=\sqrt{\gamma_c}S_z,
\]

the off-diagonal density-matrix element between \(S_z\) eigenstates \(m,m'\) decays as

\[
\rho_{m,m'}(t)
=
\rho_{m,m'}(0)
\exp\left[
-\frac{\gamma_c}{2}(m-m')^2t
\right].
\]

Populations remain fixed, but coherences decay.

## Collective dephasing vs single-particle decoherence

Collective dephasing:
- Jump operator is proportional to a collective variable, approximately \(S_z\).
- Environment learns \(N_J\).
- \(N_J\) populations are preserved.
- Off-diagonal \(N_J\)-sector coherences decay.

Single-particle dephasing:
- Local noise acts on individual atoms.
- Environment learns local state information.
- Destroys local \(|0\rangle\)-\(|1\rangle\) phase coherence.
- More damaging to the symmetric collective picture.

Single-particle spin flips:
- Local decay changes \(|1\rangle\rightarrow |0\rangle\).
- Changes effective \(S_z\).
- Reduces squeezing by changing populations.

## Independent decay

Independent decay means local jump operators such as

\[
L_i=|\downarrow_i\rangle\langle e_i|
\]

instead of a collective jump

\[
L=J_-=\sum_i|\downarrow_i\rangle\langle e_i|.
\]

Independent decay gives which-atom information and appears as single-particle decoherence in the effective \(S\)-basis.

Depending on the final state:
- \(e\rightarrow\downarrow\) mainly gives effective single-particle dephasing.
- \(e\rightarrow\uparrow\) gives effective spin flips.

## Effective \(S\)-spin picture

The paper defines

\[
|0\rangle=|\uparrow\rangle,
\]

\[
|1\rangle=
\cos\frac{\tilde\theta_J}{2}|\downarrow\rangle
+
e^{-i\tilde\phi_J}
\sin\frac{\tilde\theta_J}{2}|e\rangle.
\]

The effective \(S\)-spin is the collective spin formed from \(|0\rangle,|1\rangle\).

The many-body state can be approximately viewed as a superposition of Dicke states of this effective \(S\)-spin. The \(N_J\) sector label corresponds to the \(S_z\) coordinate.

## Generalized squeezing

The ordinary squeezing parameter assumes a two-level system and one collective spin vector. During phase 2, the state has support in three levels, so the ordinary formula is insufficient.

The paper defines a generalized squeezing parameter using:
- \(|c\rangle\): single-particle state associated with the mean Bloch vector.
- \(|j\rangle\): orthogonal fluctuation direction associated with the \(J\)-sphere.
- \(|s\rangle\): orthogonal fluctuation direction associated with the \(S\)-sphere.

The covariance matrix of four transverse fluctuation quadratures gives the generalized squeezing.

This should be used during the driven three-level part of the protocol. After the drive is turned off and \(|e\rangle\) decays away, ordinary two-level squeezing in the \(|\uparrow\rangle,|\downarrow\rangle\) basis can also be computed.

## Numerical interpretation

For exact quantum trajectories, the paper does not reduce entirely to an \(N_J\)-only Hilbert space. It keeps the internal Dicke structure inside each \(N_J\) sector.

The strong symmetry allows the operators inside each \(N_J\) sector to be smaller and independent.

The initial \(N_J\) distribution can be truncated around \(N/2\), because the initial equal superposition gives a binomial distribution.

For Fig. 4:

\[
N=1000,\qquad dN=80.
\]

The paper uses the shifted jump operator

\[
l=J_-+i\frac{\Omega}{\Gamma}
\]

because it better represents actual photon loss and minimizes jump frequency at the steady state.

## Practical coding implications

For coding:
- Keep \(N_J\) sectors explicit.
- Preserve sector coherences when simulating the full state.
- Do not treat different \(N_J\) sectors as independently normalized trajectories if inter-sector coherence matters.
- Use ensemble-averaged density matrices, not a single pure trajectory, to diagnose decoherence.
- For generalized squeezing, build collective one-body operators in the reduced basis rather than full \(3^N\) matrices.
- For Fig. 4 analytics, use fixed \(\tilde\theta_J\) with \(\cos\tilde\theta_J\approx0.5\).
