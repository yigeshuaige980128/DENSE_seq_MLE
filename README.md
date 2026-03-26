# DENSE_seq_MLE (Perl rewrite of MAP_SD)

This directory is a Perl implementation of the MAP_SD repository functionality.

## Files

- `map_sd.pl`: Main estimator (single fit, window scan, bootstrap).
- `misc_scripts/simulate_reads.pl`: Simulate distorted reads from depth profiles (stdin input).
- `misc_scripts/trim_reads.pl`: Trim/subsample FASTQ read lengths to match another FASTQ length distribution.
- `test/example_sd5e-2.txt`: Example input data.

## Input format

Tab-delimited file without header:

1. Chromosome ID
2. Position in bp
3. Map position in cM
4. Somatic count for parental allele A
5. Somatic count for parental allele a
6. Germline count for parental allele A
7. Germline count for parental allele a

## Usage

Basic:

```bash
perl map_sd.pl -d test/example_sd5e-2.txt
```

Window scan:

```bash
perl map_sd.pl -d test/example_sd5e-2.txt -w 100000
```

Bootstrap:

```bash
perl map_sd.pl -d test/example_sd5e-2.txt -b 100
```

Options (compatible with MAP_SD):

- `-d <file>` input data file
- `-r <file>` optional chromosome rate file (kept for compatibility)
- `-m <int>` minimum size (kept for compatibility)
- `-w <int>` window size
- `-b <int>` bootstrap iterations
- `-e <float>` per-read error rate
- `-t <float>` optimization tolerance

## Note

Perl runtime is required in PATH (`perl`).
