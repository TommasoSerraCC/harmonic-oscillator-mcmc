# Physical background

## From the partition function to a path integral

The quantum partition function of a system with Hamiltonian $H$ is

$$Z = \mathrm{Tr}\left(e^{-\beta H}\right) = \int \mathrm{d}x \,
\langle x | e^{-\beta H} | x \rangle .$$

Writing the matrix element as a path integral after the Wick rotation
$t \to -i\tau$ turns the oscillating weight $e^{iS}$ of real-time
quantum mechanics into a real, positive Boltzmann weight. The only
structural difference with respect to the real-time case is that the
trace forces **periodic boundary conditions** in imaginary time:

$$Z \propto \int_{x(0)=x(\beta)} \mathcal{D}x(\tau)\, e^{-S_E[x]}, \qquad
S_E[x] = \int_0^\beta \mathrm{d}\tau
\left[\frac{m}{2}\left(\frac{\mathrm{d}x}{\mathrm{d}\tau}\right)^2
+ V(x)\right].$$

Thermal averages follow the same construction,

$$\langle O \rangle =
\frac{\int \mathcal{D}x \, e^{-S_E[x]}\, O[x]}
     {\int \mathcal{D}x \, e^{-S_E[x]}},$$

so an average over quantum states becomes an average over periodic paths
weighted by $e^{-S_E}$. That is an ordinary classical statistical
mechanics problem, and it can be sampled with a Markov chain.

## Discretisation

Splitting the imaginary-time interval into $N_t$ slices of width
$a = \beta / N_t$ and introducing dimensionless variables

$$\eta = a\omega, \qquad y_i = \frac{x_i}{l}, \qquad
l = \sqrt{\frac{1}{m\omega}},$$

the Euclidean action of the harmonic oscillator becomes

$$S_E = \sum_{i=0}^{N_t-1}
\left[ y_i^2\left(\frac{\eta}{2} + \frac{1}{\eta}\right)
- \frac{1}{\eta}\, y_i y_{i+1} \right],
\qquad y_{N_t} \equiv y_0 .$$

This is a one-dimensional chain of $N_t$ coupled Gaussian variables.
Two limits matter and are probed independently:

* $\eta \to 0$ at fixed $\beta\hbar\omega$ is the **continuum limit**.
  Observables approach their physical value with $O(\eta^2)$
  corrections, so they are extrapolated with a linear fit in $\eta^2$.
* $\beta\hbar\omega \to \infty$ is the **zero-temperature limit**, where
  the system is frozen into the ground state.

Both must be respected at once to reconstruct ground-state properties:
a low temperature with a coarse lattice, or a fine lattice at a high
temperature, both fail. This is visible in the ground-state histograms
of the report.

## Observables

The primary observables are averaged over the sites of each sampled
path and then over the chain:

$$\overline{y^k} = \frac{1}{N_t}\sum_{i=0}^{N_t-1} y_i^k, \qquad k=1,2,3.$$

By the symmetry $y \to -y$ of the potential, $\langle y\rangle$ and
$\langle y^3\rangle$ vanish; measuring them is a useful check that the
chain is sampling correctly.

The energy uses the lattice virial-like estimator

$$\overline{H} = \frac{1}{2\eta}
+ \frac{1}{2N_t}\sum_i y_i^2
- \frac{1}{2N_t}\sum_i \frac{(y_{i+1}-y_i)^2}{\eta^2},$$

whose exact value in these units coincides with $\langle y^2\rangle$:

$$\langle H \rangle = \langle y^2 \rangle
= \frac{1}{2}\,\frac{1+e^{-\beta\hbar\omega}}{1-e^{-\beta\hbar\omega}}
= \frac{1}{2}\coth\!\left(\frac{\beta\hbar\omega}{2}\right).$$

## Energy gaps from connected correlators

The connected correlator of an operator $O$,

$$C_O(n) = \langle O(0) O(n) \rangle - \langle O \rangle^2,$$

decays at large separation as

$$C_O(n) \sim e^{-n\eta\,\Delta E},$$

where $\Delta E = E_{\bar m} - E_0$ is the gap to the *lowest* excited
state $\bar m$ with a non-vanishing overlap $\langle 0 | O | \bar m
\rangle$. Rather than fitting the exponential directly, the code forms
the **effective gap**

$$\Delta E(n) = \frac{1}{\eta}\,\log\frac{C(n)}{C(n+1)},$$

which develops a plateau at the physical gap once the excited states
above $\bar m$ have died out.

Choosing the interpolating operator therefore selects which gap is
measured:

* $y$ has $\langle 0 | y | n \rangle = \tfrac{1}{\sqrt2}\,\delta_{n,1}$,
  so $C_y$ isolates $E_1-E_0$ cleanly.
* $y^2$ overlaps only with states 0 and 2, giving $E_2-E_0$.
* $y^3$ overlaps with both 1 and 3, so its plateau converges slowly and
  is reached only where the correlator is already noisy.
* $A = y^3 - \tfrac{3}{2}y$ is built to be proportional to the Hermite
  polynomial $H_3$, so it overlaps only with state 3 and gives
  $E_3-E_0$. Note that constructing it requires knowing the analytic
  solution, which is possible here but not in a general problem.

## Ground-state wave function

The histogram of all sampled positions estimates the diagonal density
matrix $\rho(y,y) = \sum_n P_n |\psi_n(y)|^2$. At low enough temperature
only the ground state survives, and the histogram should reproduce

$$\lvert\psi_0(y)\rvert^2 = \frac{1}{\sqrt{\pi}}\, e^{-y^2}.$$
