# The Markov chain

The chain combines two local updates, both with acceptance probability
exactly 1, so no proposal is ever rejected and no tuning of a step size
is required.

## Heat bath

Because the discrete action is quadratic, the conditional distribution
of a single site given all the others is a Gaussian that can be sampled
directly. Isolating the terms of $S_E$ that contain $y_i$,

$$S_E \supset \alpha\, y_i^2 - \gamma\, y_i, \qquad
\alpha = \frac{\eta}{2} + \frac{1}{\eta}, \qquad
\gamma = \frac{y_{i-1} + y_{i+1}}{\eta},$$

so that

$$P(y_i \mid y_{j\neq i}) \propto
\exp\!\left[-\alpha\left(y_i - \mu\right)^2\right], \qquad
\mu = \frac{\gamma}{2\alpha}, \qquad
\sigma^2 = \frac{1}{2\alpha}.$$

The new value is drawn from that Gaussian and is completely independent
of the old one, which is why the acceptance is unity. The normal
deviates come from the Box-Muller transform applied to the `ran2`
uniform generator.

## Over-relaxation (microcanonical update)

The over-relaxation update reflects each site about the conditional mean,

$$y_i \to 2\mu - y_i .$$

Since the local action is quadratic and symmetric about $\mu$, this
leaves $S_E$ **exactly** invariant, so the move is always accepted. It
is deterministic and therefore cannot be used on its own — it does not
change the action and is not ergodic — but interleaved with the heat
bath it moves the configuration a long way across the constant-action
surface at negligible cost, which sharply reduces autocorrelation.

The invariance is not merely an argument: the test suite measures the
relative drift of the action across twenty over-relaxation sweeps and
requires it to stay at the level of floating-point round-off.

## Composition and measurement frequency

One recorded measurement corresponds to

$$10 \times \left(1 \text{ heat-bath sweep} + 5 \text{ over-relaxation sweeps}\right)$$

over the whole lattice. Because ten heat-bath sweeps separate successive
records, the stored series already has integrated autocorrelation times
of order one, which keeps the blocking analysis short and the storage
manageable.

Combining the two updates also guarantees the formal requirements of the
chain: the heat bath provides irreducibility and aperiodicity, while
over-relaxation alone would leave the action fixed.

## Sweep ordering

The Fortran implementation sweeps the lattice sequentially and updates
in place, so a site sees the already-updated value of its left neighbour.
This is a valid Markov chain: each individual site update preserves the
conditional distribution given the current state of the rest.

The NumPy re-implementation in {mod}`ui.core.ground_state_sim` instead
uses a **checkerboard** decomposition, updating all even sites and then
all odd sites in two vectorised operations. With nearest-neighbour
coupling the even sites are conditionally independent given the odd
ones, so this samples the same distribution while being fast enough to
run interactively inside the GUI.
