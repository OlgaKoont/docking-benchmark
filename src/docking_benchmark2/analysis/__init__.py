"""Analysis module - statistics and plots."""

from .statistics import calculate_protein_statistics, compare_methods
from .plots import generate_protein_plots

__all__ = ['calculate_protein_statistics', 'compare_methods', 'generate_protein_plots']

