"""dock-check: does docking work on your target, before you trust it?

A virtual screen ranks compounds. Nothing in the output tells you
whether the ranking means anything. The cheapest way to find out is to
dock a ligand the target is already known to bind, against decoys
matched to it in size, and see where it lands. If the known binder does
not beat most random molecules, the ranking is not usable and no
compound in it has been shown anything.

This runs that check and gives a verdict.
"""
__version__ = "0.1.0"

from .control import run_control
from .report import verdict, write_report

__all__ = ["run_control", "verdict", "write_report", "__version__"]
