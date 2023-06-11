import json
from typing import Any, Callable, Literal, ParamSpec, Type, TypeVar, cast

from rich.console import Console
from rich.syntax import Syntax

console = Console()

P = ParamSpec("P")
T = TypeVar("T")


# https://stackoverflow.com/a/71968448
def copy_signature(kwargs_call: Callable[P, Any]) -> Callable[[Callable[..., T]], Callable[P, T]]:
    """Decorator does nothing but returning the casted original function"""

    def return_func(func: Callable[..., T]) -> Callable[P, T]:
        return cast(Callable[P, T], func)

    return return_func


def jprint(d: Any, **kwargs: Any):
    console.print(Syntax(json.dumps(d, indent=2), "json", **kwargs))
