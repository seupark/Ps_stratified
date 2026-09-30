
import sys
import json
import pandas as pd

sys.path.insert(0, snakemake.params.modules_dir)
import modules


def _load(paths):
    if not paths:
        return pd.DataFrame()
    return pd.concat(
        [pd.read_csv(p, sep="\t", dtype={"group": str}) for p in paths],
        ignore_index=True,
    )


# Pool every TSS set so the limits are global; per-TSS figures then render on a
# shared y-axis (one range per figure type).
combined = _load(snakemake.input.combined)
fc = _load(snakemake.input.fc)

limits = modules.compute_axis_limits(combined, fc)

with open(snakemake.output.limits, "w") as f:
    json.dump(limits, f, indent=2)
