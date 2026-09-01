from __future__ import annotations

from pathlib import Path
import importlib.util


SCRIPT = Path(__file__).parents[1] / "scripts/statistics/build_stage2a_outputs.py"
SPEC = importlib.util.spec_from_file_location("stage2a", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader
SPEC.loader.exec_module(MODULE)


def test_holm_is_monotone_in_sorted_p_order():
    raw = [0.04, 0.001, 0.02, 0.5]
    adjusted = MODULE.holm(raw)
    ordered = sorted(zip(raw, adjusted))
    assert all(ordered[i][1] <= ordered[i + 1][1] for i in range(len(ordered) - 1))
    assert all(a >= p for p, a in zip(raw, adjusted))


def test_rank_biserial_direction_and_wall_clock():
    import numpy as np
    assert MODULE.rank_biserial(np.array([1.0, 2.0, 3.0])) == 1.0
    assert MODULE.rank_biserial(np.array([-1.0, -2.0, -3.0])) == -1.0
    assert MODULE.wall_seconds("1:36:14") == 5774
    assert MODULE.wall_seconds("47:50.31") == 2870.31
