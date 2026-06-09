# Inhomogeneous Couplings

This note describes the physics and high-level implementation rules for adding inhomogeneous couplings to the custom MCWF trajectory code.

The homogeneous implementation uses one active-manifold sector label,

\[
N_J=N_\downarrow+N_e,
\]

and one Dicke-basis label inside each sector,

\[
|n_e\rangle,
\qquad
n_e=0,\dots,N_J.
\]

With inhomogeneous couplings, the atoms are split into coupling groups. Global permutation symmetry is broken, but permutation symmetry remains inside each group.

## Split the atoms into two coupling groups

Split the full atom ensemble into two coupling groups,

\[
N_1+N_2=N.
\]

Here \(N_1\) is the number of physical atoms with coupling \(\omega_1\), and \(N_2\) is the number of physical atoms with coupling \(\omega_2\).

The high-level user inputs should be

\[
N,
\qquad
dN,
\qquad
N_1,
\qquad
N_2,
\qquad
\omega_1.
\]

Set

\[
\omega_2
=
\frac{N-N_1\omega_1}{N_2},
\qquad
N=N_1+N_2.
\]

Do not require the user to pass \(\omega_2\) separately. It is computed once from the physical group sizes using the atom-number weighted normalization

\[
N_1\omega_1+N_2\omega_2=N.
\]

This keeps the average coupling fixed when comparing homogeneous and inhomogeneous runs. When \(N_2=0\), \(\omega_2\) multiplies an empty physical group and is irrelevant.

The homogeneous limit is

\[
\omega_1=\omega_2=1.
\]

## Group-resolved strong-symmetry sectors

The strong-symmetry sector is no longer only a scalar \(N_J\). In inhomogeneous mode, use the pair

\[
(N_{J,1},N_{J,2})
\]

as the sector label.

Here \(N_{J,1}\) is the number of active-manifold atoms in group 1, and \(N_{J,2}\) is the number of active-manifold atoms in group 2.

The total active-manifold number is

\[
N_J=N_{J,1}+N_{J,2}.
\]

This remains the physical number of atoms in the active manifold. The coupling weights \(\omega_1\) and \(\omega_2\) affect the drive and collective jump operators, but they do not redefine the sector label.

The allowed group-resolved sectors obey

\[
0\leq N_{J,1}\leq N_1,
\qquad
0\leq N_{J,2}\leq N_2.
\]

Inside each group-resolved sector, the basis is

\[
|n_{e,1},n_{e,2}\rangle,
\]

where

\[
n_{e,1}=0,\dots,N_{J,1},
\qquad
n_{e,2}=0,\dots,N_{J,2}.
\]

The sector dimension is therefore

\[
(N_{J,1}+1)(N_{J,2}+1).
\]

This replaces the homogeneous sector basis dimension \(N_J+1\).

## High-level sector construction

Do not require the user to manually pass all low-level group-resolved sector coefficients.

The low-level representation may remain

```python
{(Nj1, Nj2): coeff}
```

but user-facing inhomogeneous runs should be able to pass

```python
N, dN, N1, N2, omega1
```

and have the code construct all allowed group-resolved sectors automatically.

First choose the total active-manifold sectors exactly as in the homogeneous code:

\[
N_J\in\{N/2-dN,\dots,N/2+dN\}.
\]

For each target homogeneous-sector label \(N_J\), include every valid pair

\[
(N_{J,1},N_{J,2})
\]

satisfying

\[
N_{J,1}+N_{J,2}=N_J,
\]

with

\[
0\leq N_{J,1}\leq N_1,
\qquad
0\leq N_{J,2}\leq N_2.
\]

Example:

\[
N=10,
\qquad
dN=1,
\qquad
N_1=3,
\qquad
N_2=7.
\]

The selected total sectors are

\[
N_J=4,5,6.
\]

The generated group-resolved sectors are therefore the valid pairs whose physical active-manifold count matches one of these values.

For example, the generated group-resolved sectors should be

\[
N_J=4:
\qquad
(0,4),(1,3),(2,2),(3,1),
\]

\[
N_J=5:
\qquad
(0,5),(1,4),(2,3),(3,2),
\]

\[
N_J=6:
\qquad
(0,6),(1,5),(2,4),(3,3).
\]

So the helper should return sector coefficients for all 12 group-resolved sectors.

The default internal state in each group-resolved sector should be

\[
|n_{e,1}=0,n_{e,2}=0\rangle,
\]

meaning all active atoms in both groups start in \(|\downarrow\rangle\).

## Initial wavefunction coefficients

The high-level sector helper should support both the physical binomial distribution and the existing square-over-total-\(N_J\) convention.

### Binomial distribution

For the product state

\[
\left(\frac{|\uparrow\rangle+|\downarrow\rangle}{\sqrt{2}}\right)^N,
\]

the probability of a group-resolved sector is

\[
P(N_{J,1},N_{J,2})
=
\frac{
\binom{N_1}{N_{J,1}}
\binom{N_2}{N_{J,2}}
}{2^N}.
\]

Therefore the unnormalized amplitude should be

\[
c_{N_{J,1},N_{J,2}}
\propto
\sqrt{
\binom{N_1}{N_{J,1}}
\binom{N_2}{N_{J,2}}
}.
\]

After restricting to the selected total-\(N_J\) window, normalize all amplitudes so that

\[
\sum_{N_{J,1},N_{J,2}}
|c_{N_{J,1},N_{J,2}}|^2
=1.
\]

### Square-over-total-\(N_J\) distribution

If `sector_distribution="square"` is used, keep the old meaning: equal total probability for each selected total \(N_J\) sector.

Do not assign equal probability to every group-resolved pair directly. Instead:

1. Give each selected total \(N_J\) equal total probability.
2. Split that probability among the valid \((N_{J,1},N_{J,2})\) pairs using the conditional binomial distribution.

Thus

\[
P(N_J)=\frac{1}{N_{\rm sectors}},
\]

where \(N_{\rm sectors}\) is the number of valid selected total-\(N_J\) sectors.

For fixed \(N_J\), use

\[
P(N_{J,1},N_{J,2}\mid N_J)
=
\frac{
\binom{N_1}{N_{J,1}}
\binom{N_2}{N_{J,2}}
}{
\sum\limits_{\substack{a+b=N_J\\0\leq a\leq N_1\\0\leq b\leq N_2}}
\binom{N_1}{a}
\binom{N_2}{b}
}.
\]

Then

\[
P(N_{J,1},N_{J,2})
=
P(N_J)
P(N_{J,1},N_{J,2}\mid N_J),
\]

and the sector amplitude is

\[
c_{N_{J,1},N_{J,2}}
=
\sqrt{P(N_{J,1},N_{J,2})}.
\]

This preserves the old square distribution over total \(N_J\), while distributing each total sector physically across the two coupling groups.

### Flat group-resolved distribution

Do not use a flat distribution over all generated \((N_{J,1},N_{J,2})\) pairs by default. That changes the intended total-\(N_J\) distribution.

If needed later, add it as a separate explicit option such as

```python
sector_distribution="square_group_resolved"
```

but do not use it for the existing `"square"` behavior.

## Group operators

Build group operators using sparse tensor products.

For group 1,

\[
J_{1,-}=J_-(N_{J,1})\otimes I_2,
\qquad
J_{1,+}=J_+(N_{J,1})\otimes I_2.
\]

For group 2,

\[
J_{2,-}=I_1\otimes J_-(N_{J,2}),
\qquad
J_{2,+}=I_1\otimes J_+(N_{J,2}).
\]

The group \(x\)-operators are

\[
J_{1,x}=\frac{J_{1,+}+J_{1,-}}{2},
\qquad
J_{2,x}=\frac{J_{2,+}+J_{2,-}}{2}.
\]

The group excited-state number operators are

\[
N_{e,1}=N_e(N_{J,1})\otimes I_2,
\qquad
N_{e,2}=I_1\otimes N_e(N_{J,2}).
\]

The total excited-state number is

\[
N_e=N_{e,1}+N_{e,2}.
\]

## Inhomogeneous drive replacement

Replace the homogeneous drive

\[
\Omega J_x
\]

by

\[
H_{\Omega}^{\rm inh}
=
\Omega(\omega_1J_{1,x}+\omega_2J_{2,x}).
\]

Equivalently,

\[
H_{\Omega}^{\rm inh}
=
\frac{\Omega}{2}
\left[
\omega_1(J_{1,+}+J_{1,-})
+
\omega_2(J_{2,+}+J_{2,-})
\right].
\]

The detuning remains homogeneous:

\[
H_\delta=-\delta(N_{e,1}+N_{e,2}).
\]

Thus, in the unshifted picture,

\[
H^{\rm inh}
=
\Omega(\omega_1J_{1,x}+\omega_2J_{2,x})
-
\delta(N_{e,1}+N_{e,2}).
\]

## Collective jump operator convention

If only the drive is made inhomogeneous, the collective decay remains the same total collective channel:

\[
L=J_{1,-}+J_{2,-}.
\]

If the inhomogeneous coupling is also meant to modify the cavity/emission coupling, use the weighted collective channel

\[
A=\omega_1J_{1,-}+\omega_2J_{2,-}.
\]

This is still collective decay because there is one shared jump channel. It is not two independent jump channels.

Do not replace this by

\[
l_1=\omega_1J_{1,-},
\qquad
l_2=\omega_2J_{2,-},
\]

unless independent group-resolved decay is explicitly requested.

The dissipator is not linear in the jump operator:

\[
\mathcal{D}[\omega_1J_{1,-}+\omega_2J_{2,-}]
\neq
\omega_1^2\mathcal{D}[J_{1,-}]
+
\omega_2^2\mathcal{D}[J_{2,-}].
\]

The cross terms encode collective emission into a shared environment.

## Shifted-jump picture

For the weighted collective model, define

\[
A=\omega_1J_{1,-}+\omega_2J_{2,-}.
\]

The shifted jump operator is

\[
l=A+i\frac{\Omega}{\Gamma}.
\]

Expanding the dissipator gives

\[
\Gamma\mathcal{D}\left[A+i\frac{\Omega}{\Gamma}\right]\rho
=
\Gamma\mathcal{D}[A]\rho
-
i\left[\Omega\frac{A^\dagger+A}{2},\rho\right].
\]

Since

\[
\frac{A^\dagger+A}{2}
=
\omega_1J_{1,x}+\omega_2J_{2,x},
\]

the shifted jump operator absorbs the full inhomogeneous drive Hamiltonian for the weighted collective model.

Therefore,

\[
H=
\Omega(\omega_1J_{1,x}+\omega_2J_{2,x})
-
\delta(N_{e,1}+N_{e,2}),
\qquad
l=A
\]

is equivalent to

\[
H=-\delta(N_{e,1}+N_{e,2}),
\qquad
l=A+i\frac{\Omega}{\Gamma}.
\]

## Coding implications

- Add a high-level helper that expands `(N, dN, N1, N2, omega1)` into all valid group-resolved sectors `(Nj1, Nj2)`.
- Keep the low-level representation `{(Nj1, Nj2): coeff}` internally.
- Choose \(\omega_2=(N-N_1\omega_1)/N_2\) once from the physical group sizes, so \(N_1\omega_1+N_2\omega_2=N\).
- Use the physical count condition \(N_{J,1}+N_{J,2}=N_J\) when mapping a target homogeneous-sector label \(N_J\) onto inhomogeneous sector pairs.
- For `sector_distribution="square"`, keep equal total probability over selected total \(N_J\) sectors and split each total sector using the conditional binomial distribution.
- For `sector_distribution="binomial"`, use amplitudes proportional to \(\sqrt{\binom{N_1}{N_{J,1}}\binom{N_2}{N_{J,2}}}\), normalized over the selected window.
- Build separate group operators for each fixed pair \((N_{J,1},N_{J,2})\).
- The Hilbert basis should track both group excitation numbers, for example \((n_{e,1},n_{e,2})\).
- Use \(J_{1,\pm}\) only on group 1 and \(J_{2,\pm}\) only on group 2.
- Use the same detuning \(\delta\) for both groups through \(-\delta(N_{e,1}+N_{e,2})\).
- Use one collective jump operator, not two independent jumps.
- Precompute sector operators, jump operators, effective generators, and propagators as in the homogeneous implementation.
- Avoid rebuilding sparse Kronecker products or \(l^\dagger l\) inside the trajectory time loop.
- Recover the homogeneous implementation when one active group is empty, e.g. \(N_{J,1}=0\) or \(N_{J,2}=0\), within numerical tolerance.


## Mean-field steady-state residual check

Add a standalone diagnostic for checking the inhomogeneous mean-field steady-state equation from `cavity_OAT.pdf`.

Use the steady-state phase equation

Use the same fixed coupling weights as the simulation:

\[
\omega_2=\frac{N-N_1\omega_1}{N_2}.
\]

Do not recompute \(\omega_2\) from instantaneous or sector-averaged values of \(N_{J,1}\) and \(N_{J,2}\).

[
\frac{\Omega\omega_a}{2}
e^{-i\phi_a}
\sin\theta_a
============

\frac{\delta}{2}
\sin\theta_a
\tan\theta_a
------------

i
\frac{\Gamma\omega_a}{4}
e^{-i\phi_a}
\sin\theta_a
\left[
\omega_1N_{J,1}e^{i\phi_1}\sin\theta_1
+
\omega_2N_{J,2}e^{i\phi_2}\sin\theta_2
\right],
\qquad
a=1,2.
]

This is the same equation as in the PDF, but written with (\theta_a) instead of enforcing a common (\theta_J). Do not replace (\theta_a) by a common angle in this diagnostic.

Define the left-hand side

[
{\rm LHS}_a
===========

\frac{\Omega\omega_a}{2}
e^{-i\phi_a}
\sin\theta_a.
]

Define the right-hand side

[
{\rm RHS}_a
===========

\frac{\delta}{2}
\sin\theta_a
\tan\theta_a
------------

i
\frac{\Gamma\omega_a}{4}
e^{-i\phi_a}
\sin\theta_a
\left[
\omega_1N_{J,1}e^{i\phi_1}\sin\theta_1
+
\omega_2N_{J,2}e^{i\phi_2}\sin\theta_2
\right].
]

The residual should be

[
R_a
===

## {\rm LHS}_a

{\rm RHS}_a.
]

Explicitly,

[
R_a
===

\frac{\Omega\omega_a}{2}
e^{-i\phi_a}
\sin\theta_a
------------

\frac{\delta}{2}
\sin\theta_a
\tan\theta_a
+
i
\frac{\Gamma\omega_a}{4}
e^{-i\phi_a}
\sin\theta_a
\left[
\omega_1N_{J,1}e^{i\phi_1}\sin\theta_1
+
\omega_2N_{J,2}e^{i\phi_2}\sin\theta_2
\right].
]

For a steady state,

[
R_a\approx0,
\qquad
a=1,2.
]

Here:

* (a=1,2) labels the coupling group.
* (N_{J,1}) is the number of active-manifold atoms belonging to group 1.
* (N_{J,2}) is the number of active-manifold atoms belonging to group 2.
* (\omega_1) and (\omega_2) are the drive weights.
* (\Omega) is the drive amplitude in the current phase.
* (\delta) is the detuning in the current phase.
* (\Gamma) is the collective decay rate.
* (\theta_a,\phi_a) are the group-resolved active-manifold Bloch angles.

The residual plotting function should evaluate (R_1(t)) and (R_2(t)) from snapshot data using the phase-local values of (\Omega(t)) and (\delta(t)), then plot the absolute values of the real and imaginary parts of the two complex residuals in a (2\times2) grid:

```text
|Re R_1|    |Im R_1|
|Re R_2|    |Im R_2|
```


### Standalone residual plotting function

Add a standalone post-processing function, without changing the core trajectory implementation unless additional data must be stored.

Suggested function name:

```python
plot_inhomogeneous_mfe_residuals(result, *, omega1, N1, N2, axes=None, output_path=None)
```

The function should:

- read saved snapshot data from the trajectory or ensemble result;
- determine the current phase for each snapshot;
- use the phase-local values of \(\Omega(t)\) and \(\delta(t)\);
- use the global \(\Gamma\) stored in the result;
- compute group-resolved observables for each snapshot;
- compute \(\theta_1(t),\phi_1(t)\) from group 1;
- compute \(\theta_2(t),\phi_2(t)\) from group 2;
- compute \(N_{J,1}(t)\) and \(N_{J,2}(t)\), or use the fixed group-resolved sector labels when evaluating within a sector;
- compute \(R_1(t)\) and \(R_2(t)\);
- plot the absolute values of the real and imaginary parts of \(R_1(t)\) and \(R_2(t)\) in a \(2\times2\) grid.

Use the \(2\times2\) grid as

```text
|Re R_1|    |Im R_1|
|Re R_2|    |Im R_2|
```

The plotted curves should be close to zero when the trajectory is near the mean-field steady state.

### Group and average angle plotting function

Add another standalone post-processing function for plotting group-resolved and average angles.

Suggested function name:

```python
plot_inhomogeneous_group_angles(result, *, axes=None, output_path=None)
```

This function should create a \(1\times2\) grid.

First panel:

- plot \(\theta_1(t)\);
- plot \(\theta_2(t)\);
- plot the correctly computed average angle \(\theta_{\rm avg}(t)\).

Second panel:

- plot \(\phi_1(t)\);
- plot \(\phi_2(t)\);
- plot the correctly computed average phase \(\phi_{\rm avg}(t)\).

The average angle must not be computed as the arithmetic average of angles:

\[
\theta_{\rm avg}
\neq
\frac{\theta_1+\theta_2}{2}.
\]

Instead, compute average Bloch components first, then convert to angles. This should follow the same logic already used in the code for computing active-manifold angles from averaged \(J_x,J_y,J_z,N_e\).

For group \(a\), compute

\[
N_{{\rm active},a}
=
2(N_{e,a}-J_{z,a}),
\]

\[
s_{x,a}
=
\frac{2J_{x,a}}{N_{{\rm active},a}},
\qquad
s_{y,a}
=
\frac{2J_{y,a}}{N_{{\rm active},a}},
\qquad
s_{z,a}
=
-\frac{2J_{z,a}}{N_{{\rm active},a}}.
\]

For the average Bloch vector, first sum the group observables,

\[
J_x^{\rm avg}
=
J_{x,1}+J_{x,2},
\qquad
J_y^{\rm avg}
=
J_{y,1}+J_{y,2},
\qquad
J_z^{\rm avg}
=
J_{z,1}+J_{z,2},
\qquad
N_e^{\rm avg}
=
N_{e,1}+N_{e,2}.
\]

Then compute

\[
N_{\rm active}^{\rm avg}
=
2(N_e^{\rm avg}-J_z^{\rm avg}),
\]

\[
s_x^{\rm avg}
=
\frac{2J_x^{\rm avg}}{N_{\rm active}^{\rm avg}},
\qquad
s_y^{\rm avg}
=
\frac{2J_y^{\rm avg}}{N_{\rm active}^{\rm avg}},
\qquad
s_z^{\rm avg}
=
-\frac{2J_z^{\rm avg}}{N_{\rm active}^{\rm avg}}.
\]

Finally,

\[
\theta_{\rm avg}
=
\arccos(s_z^{\rm avg}),
\qquad
\phi_{\rm avg}
=
\operatorname{atan2}(s_y^{\rm avg},s_x^{\rm avg}).
\]

This avoids incorrect angle averaging.

### Data requirements

The diagnostic functions require group-resolved observables. If these are not already stored or computable from snapshots, add only the minimum needed support to compute

\[
J_{x,1},J_{y,1},J_{z,1},N_{e,1},
\qquad
J_{x,2},J_{y,2},J_{z,2},N_{e,2}.
\]

Do not modify unrelated plotting or trajectory code for this feature. The residual check should be implemented as post-processing of saved snapshots whenever possible.
