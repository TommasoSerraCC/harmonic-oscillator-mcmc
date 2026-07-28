"""Presentation helpers shared by the GUI and the plotting scripts."""

import math

__all__ = ["format_value_with_uncertainty"]


def format_value_with_uncertainty(value, error, n_sig=2):
    """Round a value and its uncertainty to a common number of decimals.

    The number of decimals is chosen so that the uncertainty carries
    ``n_sig`` significant figures, which is the usual convention when
    quoting an experimental result.

    Parameters
    ----------
    value : float
        Central value.
    error : float
        Uncertainty on ``value``.
    n_sig : int, optional
        Significant figures to keep on the uncertainty. Default is 2.

    Returns
    -------
    tuple of str
        ``(value_string, error_string)``, both rounded to the same
        number of decimals. Non-finite inputs and a zero uncertainty are
        passed through unrounded, since no sensible precision exists.

    Examples
    --------
    >>> format_value_with_uncertainty(0.499612, 0.000341)
    ('0.49961', '0.00034')
    >>> format_value_with_uncertainty(1234.0, 56.0)
    ('1234', '56')
    """
    if error == 0 or not math.isfinite(error) or not math.isfinite(value):
        return f'{value}', f'{error}'
    mag = math.floor(math.log10(abs(error)))
    decimals = -int(mag) + n_sig - 1
    rounded_err = round(error, decimals)
    rounded_val = round(value, decimals)
    if decimals > 0:
        return f'{rounded_val:.{decimals}f}', f'{rounded_err:.{decimals}f}'
    return f'{int(rounded_val)}', f'{int(rounded_err)}'
