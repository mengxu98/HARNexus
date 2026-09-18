from __future__ import annotations

import importlib.util
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, Optional

_log_cache: dict[str, object] = {}


def _resolve_thisutils_log_message():
    """Dynamically resolve and load log_message from thisutils."""
    if "fn" in _log_cache:
        return _log_cache["fn"]

    explicit = os.environ.get("LOG_MESSAGE_PY")
    if explicit and Path(explicit).is_file():
        script_path = Path(explicit).resolve()
    else:
        script_path = None
        try:
            expr = (
                'p <- system.file("scripts/log_message.py", package = "thisutils"); '
                'if (!nzchar(p)) p <- system.file("python/log_message.py", package = "thisutils"); '
                'cat(p)'
            )
            out = subprocess.run(
                ["Rscript", "--vanilla", "-e", expr],
                capture_output=True,
                text=True,
                timeout=30,
            )
            candidate = out.stdout.strip()
            if out.returncode == 0 and candidate and Path(candidate).is_file():
                script_path = Path(candidate)
        except (OSError, subprocess.SubprocessError):
            pass

    if script_path and script_path.is_file():
        spec = importlib.util.spec_from_file_location("_thisutils_logger", script_path)
        if spec and spec.loader:
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            _log_cache["fn"] = module.log_message
            return module.log_message

    def fallback(message, *args, **kwargs):
        print(message, file=sys.stderr)

    _log_cache["fn"] = fallback
    return fallback


def log_message(message, *args, **kwargs):
    """Log a message using thisutils.log_message."""
    return _resolve_thisutils_log_message()(message, *args, **kwargs)

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
