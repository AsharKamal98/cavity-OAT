# Paper Key Equations

Source paper: **Engineering One Axis Twisting via a Dissipative Berry Phase Using Strong Symmetries**.

## Basic states

The physical single-particle basis is

\[
|\uparrow\rangle,\qquad |\downarrow\rangle,\qquad |e\rangle.
\]

The driven/cavity-active manifold is

\[
|\downarrow\rangle,\qquad |e\rangle.
\]

The spectator state is

\[
|\uparrow\rangle.
\]

## Collective operators

The collective lowering operator in the active manifold is

\[
J_- = \sum_i |\downarrow_i\rangle\langle e_i|.
\]

The excited-state number operator is

\[
N_e = \sum_i |e_i\rangle\langle e_i|.
\]

The conserved strong-symmetry sector label is

\[
N_J = N_\downarrow + N_e.
\]

In the \(|\downarrow\rangle,|e\rangle\) manifold,

\[
J_z = \frac{N_e-N_\downarrow}{2}.
\]

## Master equation

The effective spin master equation after eliminating the cavity is

\[
\partial_t\rho
=
-i[H_\delta,\rho]
+
\Gamma
\left(
l\rho l^\dagger
-\frac12\{l^\dagger l,\rho\}
\right).
\]

In the rotating frame of the drive,

\[
H_\delta = |\Omega|J_x-\delta N_e,
\qquad
l=J_-.
\]

Here \(\Omega\) is the effective Rabi frequency, \(\delta\) is the detuning, and

\[
\Gamma = \frac{4g_c^2}{\kappa}
\]

is the collective decay rate.

## Shifted jump operator

The same dynamics can be recast as

\[
H=-\delta N_e,
\qquad
l=J_-+i\frac{\Omega}{\Gamma}.
\]

This shifted jump operator is proportional to the output/cavity field.

At \(\delta=0\), in the spin-polarized phase, the steady state is dark in this shifted picture:

\[
\left(J_-+i\frac{\Omega}{\Gamma}\right)|\psi_{\rm ss}\rangle = 0.
\]

Equivalently,

\[
J_-|\psi_{\rm ss}\rangle
=
-i\frac{\Omega}{\Gamma}|\psi_{\rm ss}\rangle.
\]

## Critical drive and steady-state angle

For a fixed \(N_J\) sector,

\[
\Omega_c[N_J] = \frac{N_J\Gamma}{2}.
\]

At \(\delta=0\), in the spin-polarized phase \(\Omega<\Omega_c\),

\[
\cos\theta_J[N_J]
=
\sqrt{
1-\left|\frac{\Omega}{\Omega_c[N_J]}\right|^2
},
\]

and

\[
\phi_J[N_J]=\frac{\pi}{2}.
\]

The corresponding steady-state ansatz in the active manifold is

\[
|\psi_{N_J}\rangle_{\rm ss}
\approx
\left(
\cos\frac{\theta_J[N_J]}{2}|\downarrow\rangle
+
e^{-i\phi_J[N_J]}
\sin\frac{\theta_J[N_J]}{2}|e\rangle
\right)^{N_J}.
\]

## Effective qubit basis

Define

\[
|0\rangle = |\uparrow\rangle.
\]

Define

\[
|1\rangle =
\cos\frac{\tilde\theta_J}{2}|\downarrow\rangle
+
e^{-i\tilde\phi_J}
\sin\frac{\tilde\theta_J}{2}|e\rangle.
\]

The tilde means the mean-field value is evaluated around

\[
N_J = \langle N_J\rangle = N/2.
\]

In this effective basis,

\[
S_z \approx \frac{N}{2}-N_J.
\]

## Berry phase and OAT Hamiltonian

The Berry phase depends on \(N_J\). Expanding this dependence around \(N_J=N/2\) gives an effective OAT Hamiltonian:

\[
H_{\rm eff}
\approx
-\tilde\omega_B S_z
+
\check\chi S_z^2+\cdots.
\]

The OAT coefficient in the adiabatic small-\(\delta\) limit is

\[
\check\chi
=
-
\frac{\delta\sin^2\tilde\theta_J}
{2N\cos^3\tilde\theta_J}.
\]

The linear term \(-\tilde\omega_B S_z\) corresponds to a rotating frame on the effective \(S\)-Bloch sphere.

## Collective dephasing: small-\(\delta\) form

In the shifted-jump picture, the jump operator expands as

\[
\hat l
\approx
\frac{\delta}{\Gamma}\tan\tilde\theta_J
+
\frac{2\delta\sin\tilde\theta_J}
{N\Gamma\cos^3\tilde\theta_J}
S_z
+\cdots.
\]

The constant term does not affect the \(S\)-dynamics.

The linear \(S_z\) term causes collective dephasing because the environment can weakly learn

\[
S_z \approx \frac{N}{2}-N_J.
\]

## Beyond adiabatic \(\delta\)

The paper gives the generalized Holstein-Primakoff effective dynamics:

\[
H_{\rm eff}
=
-
\frac{\delta}{2N}
\frac{
N^2\Gamma^2\sin\tilde\theta_J\tan\tilde\theta_J
}{
N^2\Gamma^2\cos^2\tilde\theta_J
+
16\delta^2\sec^2\tilde\theta_J
}
S_z^2.
\]

The corresponding effective jump operator is

\[
\hat l
\approx
2e^{-i\tilde\phi_J}
\frac{
\delta\tan\tilde\theta_J
(iN\Gamma-4\delta\sec\tilde\theta_J)
}{
N^2\Gamma^2\cos^2\tilde\theta_J
+
16\delta^2\sec^2\tilde\theta_J
}
S_z.
\]

This is the effective theory benchmarked in Fig. 4.

## Single-particle decoherence

Spontaneous emission to \(|\uparrow\rangle\) acts as an effective spin flip:

\[
|1\rangle\rightarrow |0\rangle.
\]

The corresponding effective rate is

\[
\gamma_- =
\gamma_{e\uparrow}
\sin^2\frac{\tilde\theta_J}{2}.
\]

Spontaneous emission to \(|\downarrow\rangle\) acts as effective dephasing in the \(|0\rangle,|1\rangle\) basis:

\[
\gamma_d =
\gamma_{e\downarrow}
\sin^2\frac{\tilde\theta_J}{2}.
\]

## Wineland squeezing parameter

For an ordinary two-level collective spin,

\[
\xi^2 =
\frac{
N\min(\Delta S_\perp)^2
}{
|\langle S\rangle|^2
}.
\]

Here \(\min(\Delta S_\perp)^2\) is the minimum variance perpendicular to the mean spin.

## Generalized three-level squeezing parameter

For the three-level problem, use a generalized covariance matrix.

Let \(|c\rangle\) be the single-particle state associated with the Bloch vector. Let \(|j\rangle\) and \(|s\rangle\) be the two orthogonal single-particle states.

Define four collective fluctuation operators:

\[
\sum_i
\frac{|c_i\rangle\langle j_i|+|j_i\rangle\langle c_i|}{2},
\qquad
\sum_i
\frac{|c_i\rangle\langle j_i|-|j_i\rangle\langle c_i|}{2i},
\]

\[
\sum_i
\frac{|c_i\rangle\langle s_i|+|s_i\rangle\langle c_i|}{2},
\qquad
\sum_i
\frac{|c_i\rangle\langle s_i|-|s_i\rangle\langle c_i|}{2i}.
\]

Construct the associated covariance matrix \(C\). The generalized squeezing and anti-squeezing are eigenvalues of

\[
\frac{NC}{\langle N_c/2\rangle^2}.
\]

Here

\[
N_c =
\sum_i |c_i\rangle\langle c_i|.
\]

## Fig. 4 benchmark parameters

The paper uses:

\[
N=1000,
\qquad
dN=80,
\qquad
\delta=0.05N\Gamma,
\qquad
\Omega=0.465N\Gamma,
\]

corresponding to

\[
\cos\tilde\theta_J\approx0.5.
\]

They use 500 trajectories and turn off the drive at

\[
\Gamma t = 0.17.
\]

For the analytics in Fig. 4, they use a fixed steady-state value

\[
\cos\tilde\theta_J\approx0.5.
\]
