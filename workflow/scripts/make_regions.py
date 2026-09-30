
import sys
from pathlib import Path

sys.path.insert(0, snakemake.params.modules_dir)  
import modules  

bed = snakemake.input.tss_bed          
extend = int(snakemake.params.extend_len)  
min_regions = int(snakemake.params.min_regions)  
out_regions = snakemake.output.regions  

if extend > 0:
    ext_path = str(Path(out_regions).with_name(f"tss_extended_clamped_{extend}.bed"))
    bed = modules.extend_tss_clamped(bed, extend, output_file=ext_path)

range_dict = modules.create_range_dict(bed, min_regions=min_regions)
modules.save_regions(range_dict, out_regions)
