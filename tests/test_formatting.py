"""Tests for :mod:`analysis.formatting`."""

import doctest
import math

import pytest

from analysis import formatting
from analysis.formatting import format_value_with_uncertainty as fmt


@pytest.mark.parametrize('value, error, expected', [
    (0.499612, 0.000341, ('0.49961', '0.00034')),
    (0.4996, 0.0003, ('0.49960', '0.00030')),
    (2.0005, 0.0008, ('2.00050', '0.00080')),
    (1234.0, 56.0, ('1234', '56')),
    (0.9999, 0.0003, ('0.99990', '0.00030')),
])
def test_rounds_to_two_significant_figures_of_the_error(value, error, expected):
    assert fmt(value, error) == expected


def test_value_and_error_share_the_same_precision():
    value_str, error_str = fmt(0.123456, 0.00789)
    assert len(value_str.split('.')[1]) == len(error_str.split('.')[1])


def test_significant_figures_are_configurable():
    assert fmt(0.499612, 0.000341, n_sig=1) == ('0.4996', '0.0003')
    assert fmt(0.499612, 0.000341, n_sig=3) == ('0.499612', '0.000341')


@pytest.mark.parametrize('error', [0.0, math.inf, math.nan])
def test_a_meaningless_error_is_passed_through(error):
    value_str, error_str = fmt(1.5, error)
    assert value_str == '1.5'
    assert error_str == str(error)


def test_a_non_finite_value_is_passed_through():
    value_str, _ = fmt(math.nan, 0.1)
    assert value_str == 'nan'


def test_a_negative_value_keeps_its_sign():
    assert fmt(-0.0574, 0.0025) == ('-0.0574', '0.0025')


def test_docstring_examples_are_correct():
    results = doctest.testmod(formatting, verbose=False)
    assert results.failed == 0
