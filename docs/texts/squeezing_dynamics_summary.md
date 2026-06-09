# Squeezing Dynamics Summary

We define the generalized squeezing parameter from the $4\times 4$ covariance matrix $C$ as

```latex
\[
\xi^2(t)=\frac{N\,\lambda_{\min}(C)}{\langle N_c/2\rangle^2},
\qquad
C_{ab}
=
\frac{1}{2}\langle O_a O_b + O_b O_a\rangle
-\langle O_a\rangle\langle O_b\rangle.
\]
```

The four fluctuation operators are

```latex
\[
O_1=\sum_i \frac{|c_i\rangle\langle j_i|+|j_i\rangle\langle c_i|}{2},
\qquad
O_2=\sum_i \frac{|c_i\rangle\langle j_i|-|j_i\rangle\langle c_i|}{2i},
\]
\[
O_3=\sum_i \frac{|c_i\rangle\langle s_i|+|s_i\rangle\langle c_i|}{2},
\qquad
O_4=\sum_i \frac{|c_i\rangle\langle s_i|-|s_i\rangle\langle c_i|}{2i}.
\]
```

In the single-particle basis $(|u\rangle,|d\rangle,|e\rangle)$, the three dressed directions are

```latex
\[
|1\rangle
=
0|u\rangle
+
\cos\frac{\Theta_J}{2}|d\rangle
+
e^{-i\phi_J}\sin\frac{\Theta_J}{2}|e\rangle,
\]
\[
|j\rangle
=
0|u\rangle
-
\sin\frac{\Theta_J}{2}|d\rangle
+
e^{-i\phi_J}\cos\frac{\Theta_J}{2}|e\rangle,
\]
\[
|c\rangle
=
\cos\frac{\Theta_S}{2}|u\rangle
+
e^{-i\phi_S}\sin\frac{\Theta_S}{2}|1\rangle,
\]
\[
|s\rangle
=
-\sin\frac{\Theta_S}{2}|u\rangle
+
e^{-i\phi_S}\cos\frac{\Theta_S}{2}|1\rangle.
\]
```

Here $|c\rangle$ is the instantaneous single-particle mean direction at time $t$, $|j\rangle$ is the orthogonal $J$-fluctuation direction, and $|s\rangle$ is the orthogonal $S$-fluctuation direction; the minimum eigenvalue $\lambda_{\min}(C)$ gives the minimum-uncertainty direction and its variance.

The $2\times 2$ squeezing-dynamics grid shown here uses $N=1000$, $dN=80$, $100$ trajectories, a binomial initial $N_J$ distribution, $\delta=1$, $\Gamma=1$, and a phase-1 drive chosen such that

```latex
\[
\Omega = 0.2165\,\Omega_c,
\qquad
\cos\Theta_J = 0.5
\quad \text{during phase 1.}
\]
```

The upper-left panel plots $10\log_{10}(\xi^2)$. We observe squeezing already in phase 1 and continuing into phase 2; unlike the original paper discussion, we do not combine these two phases here. The best value is close to $-1.2\,\mathrm{dB}$. Stronger drive should improve the squeezing further. Once the drive is turned off, $\xi^2$ becomes slightly worse.

The upper-right panel plots the four ordered eigenvalues of $C$. All eigenvalues start at

```latex
\[
\lambda_a = \frac{N}{4},
\]
```

which corresponds to an initially circular noise distribution. During phases 1 and 2 the eigenvalues split, showing the development of anisotropic fluctuations.

The lower-left panel plots the mean active population

```latex
\[
N_c=\sum_i |c_i\rangle\langle c_i|.
\]
```

For a coherent state polarized along $|c\rangle$, one expects $\langle N_c\rangle \approx N$. Since the $dN=80$ binomial state is already close to a coherent state, this is indeed what is seen initially.

The lower-right panel plots the active-manifold excitation fraction

```latex
\[
\langle e\rangle
=
\frac{\langle N_e\rangle}{\langle N_J\rangle},
\qquad
N_J=N_d+N_e.
\]
```

This measures the excitation fraction inside the active $(|d\rangle,|e\rangle)$ manifold. For $\cos\Theta_J=0.5$ in phase 1, the steady-state expectation is

```latex
\[
\langle e\rangle
=
\sin^2\!\left(\frac{\Theta_J}{2}\right)
=
\frac{1-\cos\Theta_J}{2}
=
0.25.
\]
```

As expected, the excitation fraction drops back toward zero once the drive is turned off.

## Figure caption

```latex
\caption{
Generalized squeezing dynamics for $N=1000$ atoms with sector width $dN=80$, using a binomial initial $N_J$ distribution and averaging over $100$ MCWF trajectories. The protocol uses $\delta=1$, $\Gamma=1$, and a phase-1 drive $\Omega=0.2165\,\Omega_c$, corresponding to $\cos\Theta_J=0.5$ in phase 1. Upper left: generalized squeezing parameter $10\log_{10}(\xi^2)$, showing squeezing through phases 1 and 2 with a minimum close to $-1.2\,\mathrm{dB}$. Upper right: ordered covariance-matrix eigenvalues, which start from the isotropic value $\lambda_a=N/4$ and split as squeezing develops. Lower left: mean active population $\langle N_c\rangle$, which starts near $\langle N_c\rangle=N$ as expected for a nearly coherent initial state. Lower right: active-manifold excitation fraction $\langle e\rangle=\langle N_e\rangle/\langle N_J\rangle$, which approaches the phase-1 steady-state value $\sin^2(\Theta_J/2)=(1-\cos\Theta_J)/2=0.25$ and drops back toward zero after the drive is turned off.
}
```
