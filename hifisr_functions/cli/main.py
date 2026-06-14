"""Minimal CLI scaffold for the new architecture."""

from __future__ import annotations

import argparse

from hifisr_functions import __version__


FUNCTION_PURITY = {
    "build_parser": "pure",
    "main": "impure",
}


def build_parser():
    parser = argparse.ArgumentParser(prog="hifisr")
    parser.add_argument(
        "--version",
        action="store_true",
        help="Print the hifisr_functions version and exit.",
    )
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.version:
        print(__version__)
        return 0
    parser.print_help()
    return 0


__all__ = ["build_parser", "main"]
