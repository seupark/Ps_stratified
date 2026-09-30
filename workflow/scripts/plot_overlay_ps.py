
import sys
import json
import pandas as pd

sys.path.insert(0, snakemake.params.modules_dir)
import modules

combined = pd.concat(
    [pd.read_csv(p, sep="\t", dtype={"group": str}) for p in snakemake.input.combined],
    ignore_index=True,
)

with open(snakemake.input.limits) as f:
    limits = json.load(f)

out = snakemake.output.svg
out_base = out[:-4] if out.endswith(".svg") else out
modules.plot_overlay_ps(combined, out_base, ylim=limits.get("overlay_ps"))
