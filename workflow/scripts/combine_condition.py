
import sys
import pandas as pd

sys.path.insert(0, snakemake.params.modules_dir)  
import modules  

condition = snakemake.wildcards.condition  
ref_group = str(snakemake.params.ref_group)  

rep_dfs = [
    pd.read_csv(p, sep="\t", dtype={"group": str})
    for p in snakemake.input.curves  
]

combined = modules.combine_replicates(rep_dfs, condition)
fc = modules.fold_change_replicates(rep_dfs, condition, ref_group=ref_group)

combined.to_csv(snakemake.output.combined, sep="\t", index=False)  
fc.to_csv(snakemake.output.fc, sep="\t", index=False)              
