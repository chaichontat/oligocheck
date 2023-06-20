import logging
from typing import Any, Literal

from rich.logging import RichHandler
from rich.traceback import install


def setup_logging() -> None:
    """Setup logging."""
    install()
    if not logging.root.handlers:
        logging.basicConfig(
            level="INFO",
            format="%(message)s",
            datefmt="[%X]",
            handlers=[RichHandler(rich_tracebacks=True)],
        )


def log(*args: Any, type_: Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"] = "INFO") -> None:
    """Log a message."""
    if not logging.getLogger("rich").hasHandlers():
        setup_logging()
    logging.log(getattr(logging, type_), " ".join(str(a) for a in args))
