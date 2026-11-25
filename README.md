# sta425_subsample

Omnibenchmark **subsample** module for the STA425 replicate benchmark.

Draws `n` replicates per condition (WT and snf2) without replacement, reproducibly
from `seed`. All genes are kept, so every downstream DE result covers the same
gene set.

## Interface

```
run.R --data.counts counts.tsv.gz --data.coldata coldata.tsv \
      --n N --seed S --output_dir DIR --name ID
```

| Output | Format |
|---|---|
| `counts.tsv.gz` | same as the input count table, restricted to the drawn samples |
| `coldata.tsv` | same as the input sample sheet, restricted to the drawn samples |
| `subsample.json` | `n`, `seed`, RNG settings and the drawn sample ids per condition |

## Reproducibility

- The RNG is fixed to Mersenne-Twister with `sample.kind = "Rejection"` (the R ≥ 3.6
  default), so a seed gives the same draw on every R version.
- Sample ids are sorted before drawing, so the result does not depend on the
  column order of the input.
- The draws for WT and snf2 come from one RNG stream (WT first).

## Environment

R with `optparse` and `jsonlite` (`envs/r_base.yaml` in the benchmark repository).
