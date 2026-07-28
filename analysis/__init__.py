"""Reusable numerics for the harmonic-oscillator MCMC analysis.

The package holds the parts of the analysis that are worth testing on
their own, separated from the command line scripts and from the GUI:

``analysis.statistics``
    Estimators for correlated Markov chain data: blocking, jackknife and
    autocorrelation times.
``analysis.pipeline``
    The full "read one raw chain, write the result tables" procedure.
``analysis.formatting``
    Presentation helpers shared by the GUI and the scripts.
"""

from analysis import formatting, pipeline, statistics

__all__ = ["formatting", "pipeline", "statistics"]
