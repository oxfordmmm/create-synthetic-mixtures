# mix_pipeline

Standalone synthetic read-mix generator. Takes one challenge strain and one or
more background strains, and produces mixed FASTQs/BAMs at a range of
proportions — without running the full reference-panel/marker/detection
pipeline in `genome_analysis/`.

Each background is mixed against the challenge as its own independent
pairwise experiment (`challenge_vs_<background>`), not one N-way combined mix.

## Setup

Requires [conda](https://docs.conda.io) (or mamba) and Nextflow.

```bash
conda env create -f environment.yaml
conda activate mix_pipeline
```

This only installs Nextflow. Each process resolves its own tool environment
on the fly via `conda.enabled = true` (see `nextflow.config` and
`../envs/mapping.yml`, `../envs/mix.yml`) — you don't need to separately
install minimap2, samtools, seqtk, etc.

## Configure

Edit `params.yaml`:

```yaml
challenge:
  name:  barcode03
  fasta: /path/to/challenge_reference.fasta
  reads: /path/to/challenge_reads.fastq.gz

backgrounds:
  - name:  barcode01
    fasta: /path/to/background01_reference.fasta   # optional — omit to skip
    reads: /path/to/background01_reads.fastq.gz     # pre-filtering these reads
  - name:  barcode10
    reads: /path/to/background10_reads.fastq.gz

outdir: /path/to/output/

mix_total_reads: 50000
mix_proportions: [0, 0.0001, 0.001, 0.01, 0.1, 0.5, 1.0]
mix_replicates: 3
```

- `fasta` on a background is optional — if supplied, that background's reads
  are pre-filtered to reads mapping to it (removes off-target contamination)
  before mixing; if omitted, its raw reads are used as-is.
- Each `mix_proportions` value `p` produces a mix with proportion `p` challenge
  reads and `1-p` background reads, repeated `mix_replicates` times.

## Run

```bash
bash run.sh
```

Pass `-resume` to reuse cached results, or override any param on the command
line:

```bash
bash run.sh -resume --mix_total_reads 10000
```

To run against a different params file (e.g. the bundled smoke test):

```bash
bash run.sh -params-file params_test.yaml
```

### Output

Per-process outputs live under `${outdir}` (they use `storeDir`, so they
double as the resume cache — nothing is copied out from a working directory):

- `filtered_reads/<name>/` — reads filtered to their own reference
- `index/<name>/` — minimap2 index of the challenge reference
- `mixes/<challenge>_vs_<background>/reads/` — mixed FASTQs
- `mixes/<challenge>_vs_<background>/bam/` — mixed reads mapped to the challenge reference
- `logs/` — Nextflow timeline/report/trace HTML+text reports

Nextflow's own run log is written separately to `./logs/nextflow_<timestamp>.log`
(next to `run.sh`, not under `outdir`).

## Quick test run

`params_test.yaml` points at small synthetic data in `test_data/` (a 20kb
random challenge genome + background genome, 200 reads each) with a tiny
mix sweep, so it finishes in seconds and has no dependency on real assembly
data:

```bash
bash run.sh -params-file params_test.yaml
```
