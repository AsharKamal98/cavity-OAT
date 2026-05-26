# Instruction: Implement Generalized Three-Level Squeezing Parameter

Construct the generalized three-level squeezing parameter for the custom MCWF code. The goal is to compute \\(\xi^2_{\rm gen}(t)\\) at each saved timestep after the simulation has run.

The squeezing calculation should preferably be implemented as a standalone post-processing function that takes the simulation result or observable output as input. If it is simpler, it may be integrated into the observable function, but it must be optional via a boolean flag, since this calculation may be expensive and is not always needed.

If additional per-timestep data must be saved for this calculation, add a boolean flag such as `save_squeezing_data` or `compute_generalized_squeezing`.

---

## Per-timestep calculation

At each timestep, do the following.

### 1. Define the effective dressed state \\(|1\rangle\\)

Work in the single-particle basis

\\[
(|u\rangle, |d\rangle, |e\rangle).
\\]

Use the ansatz

\\[
|1\rangle =
0|u\rangle
+
\cos\frac{\theta_J}{2}|d\rangle
+
e^{-i\phi_J}\sin\frac{\theta_J}{2}|e\rangle.
\\]

Here \\(\theta_J,\phi_J\\) should be found by comparing this ansatz to the simulation state at the current timestep. Use the existing logic for computing average active-manifold angles, i.e. the logic based on normalized \\(s_i\\)-like quantities rather than raw \\(J_i\\), since different \\(N_J\\) sectors have different spin lengths and raw \\(J_i\\) values are not directly comparable across sectors.

Reuse the existing `active_manifold_angles` logic where appropriate.

---

### 2. Define the instantaneous mean direction \\(|c\rangle\\)

Use the ansatz

\\[
|c\rangle
=
\cos\frac{\theta_S}{2}|u\rangle
+
e^{-i\phi_S}\sin\frac{\theta_S}{2}|1\rangle.
\\]

Here \\(\theta_S,\phi_S\\) should be found by comparing this ansatz to the instantaneous effective \\(S\\)-Bloch vector of the simulation state. The angles \\(\theta_J,\phi_J\\) appearing inside \\(|1\rangle\\) are the angles found in step 1.

In the explicit \\((u,d,e)\\) basis,

\\[
|c\rangle =
\begin{pmatrix}
\cos(\theta_S/2) \\
e^{-i\phi_S}\sin(\theta_S/2)\cos(\theta_J/2) \\
e^{-i(\phi_S+\phi_J)}\sin(\theta_S/2)\sin(\theta_J/2)
\end{pmatrix}.
\\]

---

### 3. Define the \\(J\\)-fluctuation direction \\(|j\rangle\\)

Using the \\(J\\)-angles found in step 1, construct the state orthogonal to \\(|1\rangle\\) inside the \\((d,e)\\) manifold:

\\[
|j\rangle =
0|u\rangle
-
\sin\frac{\theta_J}{2}|d\rangle
+
e^{-i\phi_J}\cos\frac{\theta_J}{2}|e\rangle.
\\]

In vector form:

\\[
|j\rangle =
\begin{pmatrix}
0 \\
-\sin(\theta_J/2) \\
e^{-i\phi_J}\cos(\theta_J/2)
\end{pmatrix}.
\\]

---

### 4. Define the \\(S\\)-fluctuation direction \\(|s\rangle\\)

Using the \\(S\\)-angles from step 2 and the \\(J\\)-angles from step 1, construct

\\[
|s\rangle =
-\sin\frac{\theta_S}{2}|u\rangle
+
e^{-i\phi_S}\cos\frac{\theta_S}{2}|1\rangle.
\\]

In the \\((u,d,e)\\) basis,

\\[
|s\rangle =
\begin{pmatrix}
-\sin(\theta_S/2) \\
e^{-i\phi_S}\cos(\theta_S/2)\cos(\theta_J/2) \\
e^{-i(\phi_S+\phi_J)}\cos(\theta_S/2)\sin(\theta_J/2)
\end{pmatrix}.
\\]

The three states \\(|c\rangle,|j\rangle,|s\rangle\\) should form an orthonormal single-particle basis up to numerical precision.

---

### 5. Construct the four local fluctuation operators

Construct the four \\(3\times3\\) single-particle operators:

\\[
o_1 =
\frac{|c\rangle\langle j|+|j\rangle\langle c|}{2},
\\]

\\[
o_2 =
\frac{|c\rangle\langle j|-|j\rangle\langle c|}{2i},
\\]

\\[
o_3 =
\frac{|c\rangle\langle s|+|s\rangle\langle c|}{2},
\\]

\\[
o_4 =
\frac{|c\rangle\langle s|-|s\rangle\langle c|}{2i}.
\\]

Then construct the corresponding collective operators

\\[
O_a = \sum_i o_a^{(i)},
\qquad a=1,2,3,4.
\\]

Do not construct full \\(3^N\\)-dimensional tensor-product operators unless absolutely necessary. Use the existing reduced \\((N_J,n_e)\\) basis and exploit symmetry. Each collective operator should be represented as a sparse operator acting on the reduced simulation basis.

If possible, construct \\(O_a\\) from precomputed collective one-body transition operators

\\[
A_{\mu\nu}=\sum_i |\mu_i\rangle\langle \nu_i|,
\qquad
\mu,\nu\in\{u,d,e\}.
\\]

Then

\\[
O_a =
\sum_{\mu,\nu}
(o_a)_{\mu\nu} A_{\mu\nu}.
\\]

This avoids building large tensor-product matrices.

---

### 6. Construct the covariance matrix

For the current state \\(|\psi(t)\rangle\\), compute

\\[
\mu_a = \langle O_a\rangle.
\\]

Then construct the \\(4\times4\\) covariance matrix

\\[
C_{ab}
=
\frac{1}{2}
\langle O_aO_b+O_bO_a\rangle
-
\mu_a\mu_b.
\\]

For efficiency, avoid explicitly constructing \\(O_aO_b\\). Instead use

\\[
C_{ab}
=
\mathrm{Re}\left[
\langle O_a\psi | O_b\psi\rangle
\right]
-
\mu_a\mu_b.
\\]

---

### 7. Minimum fluctuation direction

Diagonalize the \\(4\times4\\) covariance matrix and take

\\[
\lambda_{\min}(C).
\\]

This is the minimum generalized transverse variance.

---

### 8. Generalized squeezing parameter

Compute

\\[
\xi^2_{\rm gen}(t)
=
\frac{N\lambda_{\min}(C)}
{\langle N_c/2\rangle^2},
\\]

where

\\[
N_c =
\sum_i |c_i\rangle\langle c_i|.
\\]

Compute \\(\langle N_c\rangle\\) using the same collective-operator construction:

\\[
N_c =
\sum_{\mu,\nu}
(|c\rangle\langle c|)_{\mu\nu} A_{\mu\nu}.
\\]

If the state is well polarized along \\(|c\rangle\\), then \\(\langle N_c/2\rangle\approx N/2\\), but the code should compute it explicitly rather than assuming this.

---

## Data that may need to be saved per timestep

Save enough data to reconstruct the squeezing calculation after the simulation:

- the full state vector \\(|\psi(t)\rangle\\) in the reduced \\((N_J,n_e)\\) basis, or equivalent data sufficient to reconstruct it;
- the timestep values \\(t\\);
- \\(N\\), \\(\Gamma\\), \\(\Omega\\), \\(\delta\\), and any sector truncation metadata;
- the basis/index mapping for \\((N_J,n_e)\\);
- if not recomputed during post-processing, save \\(\theta_J(t),\phi_J(t),\theta_S(t),\phi_S(t)\\).

Prefer recomputing angles from the state during post-processing if this avoids storing redundant data.
