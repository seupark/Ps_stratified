
import sys

sys.path.insert(0, snakemake.params.modules_dir)  
import modules  

range_dict = modules.load_regions(snakemake.input.regions)
_max_group = snakemake.params.max_group
_max_group = None if _max_group is None else int(_max_group)
df = modules.generate_df_from_cooler(
    snakemake.input.mcool,
    int(snakemake.wildcards.res),
    range_dict,
    sample_name=snakemake.wildcards.sample,
    max_group=_max_group,
)
df.to_csv(snakemake.output.curve, sep="\t", index=False)
