
import sys
import json
import pandas as pd

sys.path.insert(0, snakemake.params.modules_dir)
import modules  #


def _load(paths):
    return pd.concat(
        [pd.read_csv(p, sep="\t", dtype={"group": str}) for p in paths],
        ignore_index=True,
    )


combined = _load(snakemake.input.combined)
fc = _load(snakemake.input.fc)

with open(snakemake.input.limits) as f:
    limits = json.load(f)

out = snakemake.output.svg
out_base = out[:-4] if out.endswith(".svg") else out
modules.plot_per_condition(combined, fc, out_base,
                           main_ylim=limits.get("per_condition_main"),
                           fc_ylim=limits.get("per_condition_fc"))
