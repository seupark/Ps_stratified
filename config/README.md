# Configuration

Edit [`config.yaml`](config.yaml) to point the pipeline at the data and tune the analysis.
Paths are resolved **relative to the directory Snakemake is run from** (the repository root).

## Fields

### Paths
| Key | Meaning |
|-----|---------|
| `mcool_dir` | Directory holding the balanced multi-resolution `.mcool` files. Default `resources/mcool_files`. |
| `results_dir` | Where outputs are written. Default `results`. |
| `tss_beds` | Map of `{tag: bed_path}`; each tag namespaces a full result set, so several TSS BEDs are processed in one run. BEDs live under `resources/bed_files/`. A single `tss_bed` key is still accepted. |

### Analysis parameters
| Key | Meaning |
|-----|---------|
| `resolutions` | One or more bin sizes to process; each **must exist** as a zoom level in the `.mcool` files. |
| `extend` | `0` = raw TSS intervals; `>0` = clamped, overlap-free extension (bp). |
| `min_regions` | Drop chromosomes with fewer inter-TSS regions than this. |
| `max_group` | TSS crossings above this collapse into `">N"`; `null` = no collapse. |
| `ref_group` | Fold-change reference group (`"0"` = within-region contacts). |

### Experimental design
`conditions` maps each condition to its replicates:

```yaml
conditions:
  <condition>:
    <replicate_sample>: <mcool_filename>   # filename is relative to mcool_dir
```

Filenames are listed explicitly because the source `.mcool` files may **not be named uniformly**.
