"""
Python implementation of log_message function
Provides formatted logging with timestamps, colors, and message types
"""

import sys
import os
import re
import inspect
from datetime import datetime
from typing import Optional, List


class LogMessage:
    """Log message formatter with color and style support"""

    # ANSI color codes
    COLORS = {
        "black": "\033[30m",
        "red": "\033[31m",
        "green": "\033[32m",
        "yellow": "\033[33m",
        "blue": "\033[34m",
        "magenta": "\033[35m",
        "cyan": "\033[36m",
        "white": "\033[37m",
        "grey": "\033[90m",
        "orange": "\033[38;2;255;165;0m",
        "br_red": "\033[91m",
        "br_green": "\033[92m",
        "br_yellow": "\033[93m",
        "br_blue": "\033[94m",
        "br_magenta": "\033[95m",
        "br_cyan": "\033[96m",
        "br_white": "\033[97m",
        "none": "\033[0m",
    }

    # Background colors
    BG_COLORS = {
        "black": "\033[40m",
        "red": "\033[41m",
        "green": "\033[42m",
        "yellow": "\033[43m",
        "blue": "\033[44m",
        "magenta": "\033[45m",
        "cyan": "\033[46m",
        "white": "\033[47m",
        "none": "\033[0m",
    }

    # Text styles
    STYLES = {
        "bold": "\033[1m",
        "dim": "\033[2m",
        "italic": "\033[3m",
        "underline": "\033[4m",
        "strikethrough": "\033[9m",
        "inverse": "\033[7m",
    }

    # Message type symbols and colors
    MESSAGE_TYPES = {
        "info": {"symbol": "ℹ", "color": "blue"},
        "success": {"symbol": "✓", "color": "green"},
        "warning": {"symbol": "!", "color": "yellow"},
        "error": {"symbol": "✗", "color": "red"},
        "running": {"symbol": "◌", "color": "orange"},
    }

    RESET = "\033[0m"

    # Inline format styles (cli-style)
    INLINE_FORMATS = {
        ".pkg": {"color": "blue", "style": []},
        ".code": {"color": "grey", "style": []},
        ".val": {"color": "blue", "style": []},
        ".arg": {"color": "none", "style": []},
        ".fun": {"color": "none", "style": []},
        ".file": {"color": "blue", "style": []},
        ".path": {"color": "blue", "style": []},
        ".field": {"color": "blue", "style": []},
        ".emph": {"color": "none", "style": ["italic"]},
        ".strong": {"color": "none", "style": ["bold"]},
    }

    def __init__(self):
        self._check_color_support()

    def _check_color_support(self):
        """Check if colors should be enabled (prefer enabled unless NO_COLOR is set)."""
        # disable colors (e.g., R's reticulate output), unless explicitly set NO_COLOR
        self.color_support = os.environ.get("NO_COLOR") is None

    def _hex_to_rgb(self, hex_color: str) -> tuple:
        """Convert hex color to RGB tuple"""
        hex_color = hex_color.lstrip("#")
        if len(hex_color) != 6:
            raise ValueError(f"Invalid hex color: #{hex_color}")
        return tuple(int(hex_color[i : i + 2], 16) for i in (0, 2, 4))

    def _rgb_to_ansi(self, rgb: tuple, bg: bool = False) -> str:
        """Convert RGB to ANSI color code"""
        r, g, b = rgb
        if bg:
            return f"\033[48;2;{r};{g};{b}m"
        else:
            return f"\033[38;2;{r};{g};{b}m"

    def _apply_color(self, text: str, color: Optional[str], bg: bool = False) -> str:
        """Apply color to text"""
        if not self.color_support or not color:
            return text

        # Handle hex colors
        if color.startswith("#"):
            try:
                rgb = self._hex_to_rgb(color)
                ansi_color = self._rgb_to_ansi(rgb, bg)
                return f"{ansi_color}{text}{self.RESET}"
            except ValueError:
                return text

        # Handle named colors
        color_map = self.BG_COLORS if bg else self.COLORS
        if color in color_map:
            return f"{color_map[color]}{text}{self.RESET}"

        return text

    def _apply_style(self, text: str, styles: Optional[List[str]]) -> str:
        """Apply text styles"""
        if not self.color_support or not styles:
            return text

        style_codes = []
        for style in styles:
            if style in self.STYLES:
                style_codes.append(self.STYLES[style])

        if style_codes:
            return f"{''.join(style_codes)}{text}{self.RESET}"

        return text

    def _get_indent(self, level: int, symbol: str) -> str:
        """Generate indentation string"""
        if symbol != "  ":
            return symbol * level + " "
        elif level > 1:
            return "  " * (level - 1)
        else:
            return ""

    def _format_message(
        self,
        message: str,
        message_type: str = "info",
        timestamp: bool = True,
        timestamp_format: str = "%Y-%m-%d %H:%M:%S",
        level: int = 1,
        symbol: str = "  ",
        text_color: Optional[str] = None,
        back_color: Optional[str] = None,
        text_style: Optional[List[str]] = None,
        multiline_indent: bool = False,
        timestamp_style: bool = True,
    ) -> str:
        """Format the complete message"""

        # Get message type info
        msg_info = self.MESSAGE_TYPES.get(message_type, self.MESSAGE_TYPES["info"])
        msg_symbol = msg_info["symbol"]
        msg_color = msg_info["color"]

        # Build timestamp
        timestamp_str = ""
        if timestamp:
            timestamp_str = f"[{datetime.now().strftime(timestamp_format)}] "

        # Build indentation
        indent = self._get_indent(level, symbol)

        # Prepare message-type symbol (with color)
        symbol_colored = (
            self._apply_color(msg_symbol, msg_color)
            if self.color_support
            else msg_symbol
        )
        symbol_part = f"{symbol_colored} "

        # Handle multiline messages
        if "\n" in message:
            lines = message.split("\n")
            formatted_lines = []

            for i, line in enumerate(lines):
                if i == 0 or multiline_indent:
                    # First line or multiline_indent=True: full formatting (symbol before timestamp)
                    prefix = symbol_part + timestamp_str + indent
                else:
                    # Subsequent lines: alignment spaces + indent
                    alignment_spaces = (
                        (" " * (len(symbol_part) + len(timestamp_str)))
                        if timestamp
                        else (" " * len(symbol_part))
                    )
                    prefix = alignment_spaces + indent

                # Apply formatting to line
                formatted_line = self._apply_formatting(
                    line, text_color, back_color, text_style, timestamp_style, prefix
                )
                formatted_lines.append(formatted_line)

            return "\n".join(formatted_lines)

        # Single line message (symbol before timestamp)
        prefix = symbol_part + timestamp_str + indent

        # Apply formatting
        formatted_msg = self._apply_formatting(
            message, text_color, back_color, text_style, timestamp_style, prefix
        )

        # Already added colored symbol in prefix
        return formatted_msg

    def _apply_formatting(
        self,
        text: str,
        text_color: Optional[str],
        back_color: Optional[str],
        text_style: Optional[List[str]],
        timestamp_style: bool,
        prefix: str,
    ) -> str:
        """Apply all formatting to text"""

        # Apply styles first
        if text_style:
            text = self._apply_style(text, text_style)

        # Apply text color
        if text_color:
            text = self._apply_color(text, text_color)

        # Apply background color
        if back_color:
            text = self._apply_color(text, back_color, bg=True)

        return prefix + text

    def _parse_inline_expressions(self, message: str, caller_frame=None) -> str:
        """
        Parse inline expressions with cli-style formatting

        Supports:
        - {.pkg package_name} / {pkg package_name}
        - {.pkg {variable}} / {pkg {variable}}
        - {.code some_code} / {code some_code}
        - {.val {expression}} / {val {expression}}
        - {expression}  # bare expression evaluation, no formatting
        etc.
        """
        max_iterations = 15
        iteration = 0

        allowed_tags = set(t.lstrip(".") for t in self.INLINE_FORMATS.keys())

        def eval_expression(expr: str) -> str:
            expr = expr.strip()
            if not expr:
                return ""
            if caller_frame is None:
                return expr
            try:
                value = eval(expr, caller_frame.f_globals, caller_frame.f_locals)
                return str(value)
            except Exception:
                return expr

        while iteration < max_iterations:
            iteration += 1

            # 1) first parse the inline format: {.tag content} or {tag content}
            # content can be non-curly brace text or single layer {expr}
            format_match = re.search(
                r"\{\.?([A-Za-z_][A-Za-z0-9_]*)\s+(\{[^{}]*\}|[^{}]+)\}", message
            )

            if format_match:
                tag = format_match.group(1)
                content = format_match.group(2)

                # only process allowed tags
                if tag in allowed_tags:
                    if content.startswith("{") and content.endswith("}"):
                        evaluated = eval_expression(content[1:-1])
                    else:
                        evaluated = content

                    formatted = self._apply_inline_format(evaluated, f".{tag}")
                    message = (
                        message[: format_match.start()]
                        + formatted
                        + message[format_match.end() :]
                    )
                    # continue to next iteration
                    continue
                else:
                    # non-supported tag, skip to bare expression stage
                    pass

            # 2) then parse the bare expression: {expr} (not starting with .tag or tag)
            # to avoid conflict with format syntax, here exclude {.xxx ...} and {xxx ...}
            bare_match = re.search(r"\{([^{}]+)\}", message)
            if bare_match:
                inner = bare_match.group(1).strip()
                # if it looks like a format prefix (.tag or tag followed by space), skip this iteration
                if re.match(r"^\.?[A-Za-z_][A-Za-z0-9_]*\s+", inner):
                    # no bare expression to process, end loop
                    break
                evaluated = eval_expression(inner)
                message = (
                    message[: bare_match.start()]
                    + evaluated
                    + message[bare_match.end() :]
                )
                continue

            # no content to parse, end loop
            break

        return message

    def _apply_inline_format(self, text: str, format_type: str) -> str:
        """Apply formatting to inline content"""
        if format_type not in self.INLINE_FORMATS:
            return text

        fmt = self.INLINE_FORMATS[format_type]

        # apply color
        if fmt["color"] and fmt["color"] != "none":
            text = self._apply_color(text, fmt["color"])

        # apply style
        if fmt["style"]:
            text = self._apply_style(text, fmt["style"])

        return text

    def _format_traceback(self, depth: int = 1, skip_frames: int = 3) -> str:
        """
        Format traceback information for error messages

        Parameters:
        - depth: number of stack frames to show
        - skip_frames: number of internal frames to skip

        Returns:
        - formatted stack information string
        """
        stack = inspect.stack()

        # skip log_message internal call frames
        # skip_frames: _format_traceback, log_message, _format_message, etc.
        start_idx = skip_frames
        end_idx = min(start_idx + depth, len(stack))

        lines = []
        for i in range(start_idx, end_idx):
            frame_info = stack[i]
            # format file location information
            location = f'  File "{os.path.basename(frame_info.filename)}", line {frame_info.lineno}, in {frame_info.function}'
            lines.append(location)

            # show code line if available
            if frame_info.code_context:
                code_line = frame_info.code_context[0].strip()
                lines.append(f"    {code_line}")

        return "\n".join(lines) if lines else ""


# Global instance
_logger = LogMessage()


def log_message(
    *args,
    verbose: bool = True,
    message_type: str = "info",
    timestamp: bool = True,
    timestamp_format: str = "%Y-%m-%d %H:%M:%S",
    level: int = 1,
    symbol: str = "  ",
    text_color: Optional[str] = None,
    back_color: Optional[str] = None,
    text_style: Optional[List[str]] = None,
    multiline_indent: bool = False,
    timestamp_style: bool = True,
    show_traceback: bool = True,
    traceback_depth: int = 1,
) -> None:
    """
    Print formatted message with timestamp, colors, styling, and inline expressions

    Parameters:
    -----------
    *args : str
        Message parts to concatenate. Supports cli-style inline expressions:
        - {.pkg package_name} - Package names (cyan + bold)
        - {.code code_snippet} - Code snippets (grey)
        - {.val variable_name} - Variable values (blue)
        - {.arg parameter_name} - Function parameters (green)
        - {.fun function_name} - Function names (magenta)
        - {.file file_path} - File paths (underline)
        - {.path directory_path} - Directory paths (underline)
        - {.field field_name} - Object fields (cyan)
        - {.emph text} - Emphasized text (italic)
        - {.strong text} - Strong text (bold)

        Expressions can be nested: {.pkg {package_name}}
    verbose : bool, default True
        Whether to print the message
    message_type : str, default "info"
        Type of message: "info", "success", "warning", "error", "running"
    timestamp : bool, default True
        Whether to show timestamp
    timestamp_format : str, default "%Y-%m-%d %H:%M:%S"
        Timestamp format string
    level : int, default 1
        Indentation level
    symbol : str, default "  "
        Symbol used for indentation
    text_color : str, optional
        Text color (named color or hex code)
    back_color : str, optional
        Background color (named color or hex code)
    text_style : list, optional
        Text styles: ["bold", "italic", "underline", "dim", "strikethrough", "inverse"]
    multiline_indent : bool, default False
        Whether to apply formatting to each line in multiline messages
    timestamp_style : bool, default True
        Whether to apply styling to timestamp
    show_traceback : bool, default True
        Whether to show traceback for error messages
    traceback_depth : int, default 1
        Number of stack frames to show in traceback

    Returns:
    --------
    None
    """

    if not verbose:
        return

    # Validate message_type
    valid_types = ["info", "success", "warning", "error", "running"]
    if message_type not in valid_types:
        message_type = "info"

    # Get caller frame for inline expression evaluation and traceback
    # Need to handle both direct calls and calls through convenience functions
    current_frame = inspect.currentframe()
    caller_frame = current_frame.f_back

    # If called through convenience function (log_info, log_error, etc.), skip one more frame
    if caller_frame and caller_frame.f_code.co_name in [
        "log_info",
        "log_success",
        "log_warning",
        "log_error",
        "log_running",
    ]:
        caller_frame = caller_frame.f_back

    # Build message from args
    if not args:
        message = ""
    else:
        message = "".join(str(arg) for arg in args)

    # Parse inline expressions
    message = _logger._parse_inline_expressions(message, caller_frame)

    # Format and print message
    formatted_message = _logger._format_message(
        message=message,
        message_type=message_type,
        timestamp=timestamp,
        timestamp_format=timestamp_format,
        level=level,
        symbol=symbol,
        text_color=text_color,
        back_color=back_color,
        text_style=text_style,
        multiline_indent=multiline_indent,
        timestamp_style=timestamp_style,
    )

    # Add traceback for error messages
    if message_type == "error" and show_traceback:
        traceback_info = _logger._format_traceback(depth=traceback_depth)
        if traceback_info:
            formatted_message += "\n" + traceback_info

    print(formatted_message)


# Convenience functions for common message types
def log_info(*args, **kwargs):
    """Log info message"""
    log_message(*args, message_type="info", **kwargs)


def log_success(*args, **kwargs):
    """Log success message"""
    log_message(*args, message_type="success", **kwargs)


def log_warning(*args, **kwargs):
    """Log warning message"""
    log_message(*args, message_type="warning", **kwargs)


def log_error(*args, **kwargs):
    """Log error message"""
    log_message(*args, message_type="error", **kwargs)


def log_running(*args, **kwargs):
    """Log running message"""
    log_message(*args, message_type="running", **kwargs)
