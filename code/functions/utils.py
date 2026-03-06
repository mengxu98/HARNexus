from __future__ import annotations

import os
from typing import Dict

from functions.log_message import log_message

__all__ = [
    "log_message",
    "COLOR_CELLTYPES",
    "COLOR_STAGES",
    "DEFAULT_NODE_COLOR",
    "check_dir",
]

COLOR_CELLTYPES: Dict[str, str] = {
    "Radial glia": "#8076A3",
    "Neuroblasts": "#ED5736",
    "Excitatory neurons": "#0AA344",
    "Inhibitory neurons": "#2177B8",
    "Astrocytes": "#D70440",
    "Oligodendrocyte progenitor cells": "#F9BD10",
    "Oligodendrocytes": "#B14B28",
    "Microglia": "#006D87",
    "Endothelial cells": "#5E7987",
}


def _hex_to_rgb_plot(h: str):
    h = h.lstrip("#")
    return tuple(int(h[i : i + 2], 16) for i in (0, 2, 4))


def _rgb_to_hex_plot(rgb) -> str:
    return "#{:02X}{:02X}{:02X}".format(*[int(round(x)) for x in rgb])


def _build_color_stages() -> Dict[str, str]:
    c1 = _hex_to_rgb_plot("#0AA344")
    c2 = _hex_to_rgb_plot("#006D87")
    part1 = [
        _rgb_to_hex_plot(tuple(c1[i] + (c2[i] - c1[i]) * (k / 6.0) for i in range(3)))
        for k in range(7)
    ]
    c3 = _hex_to_rgb_plot("#2B73AF")
    c4 = _hex_to_rgb_plot("#003D74")
    part2 = [
        _rgb_to_hex_plot(tuple(c3[i] + (c4[i] - c3[i]) * (k / 7.0) for i in range(3)))
        for k in range(8)
    ]
    return {f"S{i + 1}": c for i, c in enumerate(part1 + part2)}


COLOR_STAGES = _build_color_stages()
DEFAULT_NODE_COLOR = "rgba(200, 200, 200, 0.8)"


def check_dir(path: str) -> str:
    """Create a directory if needed and return the normalized path."""
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)
    return path
