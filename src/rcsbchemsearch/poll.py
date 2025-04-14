import errno
import os
import time
from collections.abc import Callable
from datetime import timedelta
from typing import Any, TypeVar

T = TypeVar("T")


def poll_until(
    resource: Callable[[], T | None],
    *,
    timeout: float | timedelta = 600,
    poll_every: float | timedelta = 60,
    predicate: Callable[[T], bool] = bool,
    sleep: Callable[[float], Any] = time.sleep,
) -> T:
    """
    Polls `t := resource()` until `predicate(t)` returns `True`, raising `TimeoutError` if a timeout is reached.

    :param resource: A function that takes no args and returns a resource.
    :param timeout: The timeout, as a `timedelta` or in seconds.
           This function is guaranteed to exit before the timeout, subject to how long `resource()` takes.
           If `resource()` is called `n` seconds before timeout `T` but takes `n + x` seconds for some positive `x`,
           then this function will exit at `T + x` seconds.
    :param poll_every: The duration (as a `timedelta` or in seconds) between consecutive `resource()` calls.
           This is the time since the last call, **not** the time since the result was received.
           If `resource()` takes longer than `poll_every`, then there will be no wait period.
    :param predicate: A boolean function called on each `resource()` result to decide whether it was successful.
           The default returns `False` for any falsy `T` value such as `""`, `0`, and `None`.
           If you want to consider any non-`None` value successful, use `lambda t: t is not None`.
           (Technically, will return `v` for any truthy `predicate(v)`, not just `predicate(v) is True`.)
    :param sleep: A function that waits, normally `time.sleep`.
           When calling from an `async` context, pass `sleep=asyncio.sleep` to free resources during sleep.

    :returns: The value returned from `resource()`.
    :raises TimeoutError: If the timeout is reached.
    """
    if isinstance(timeout, timedelta):
        timeout = timeout.total_seconds()
    if isinstance(poll_every, timedelta):
        poll_every = poll_every.total_seconds()
    t_start = time.monotonic()
    t_end = t_start + timeout
    while (t0 := time.monotonic()) < t_end:
        if predicate(t := resource()):
            return t
        # Sleep for `poll_every` minus the time that `resource()` took, but cap to the timeout.
        # `0.0` is needed because `time.sleep(-1)` errors.
        t1 = time.monotonic()
        sleep(min(0.0, poll_every - t1 + t0, t_end - t1))
    raise TimeoutError(errno.ETIMEDOUT, os.strerror(errno.ETIMEDOUT))
