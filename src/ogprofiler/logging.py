"""Contextual console and file logging."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

_FORMAT = (
    "%(asctime)s %(levelname)s stage=%(stage)s task=%(task)s component=%(component)s %(message)s"
)


class ContextLogger(logging.LoggerAdapter[logging.Logger]):
    def process(self, msg: object, kwargs: Any) -> tuple[object, Any]:
        context: dict[str, object] = {"stage": "-", "task": "-", "component": "-"}
        context.update(self.extra or {})
        supplied = kwargs.get("extra", {})
        context.update(supplied)
        kwargs["extra"] = context
        return msg, kwargs

    def bind(self, **context: object) -> ContextLogger:
        merged = dict(self.extra or {})
        merged.update(context)
        return ContextLogger(self.logger, merged)


def configure_logging(log_file: Path, level: str = "INFO") -> ContextLogger:
    logger = logging.getLogger("ogprofiler")
    logger.setLevel(level.upper())
    logger.propagate = False
    logger.handlers.clear()
    formatter = logging.Formatter(_FORMAT)
    console = logging.StreamHandler()
    console.setFormatter(formatter)
    logger.addHandler(console)
    log_file.parent.mkdir(parents=True, exist_ok=True)
    file_handler = logging.FileHandler(log_file, encoding="utf-8")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    return ContextLogger(logger, {})
