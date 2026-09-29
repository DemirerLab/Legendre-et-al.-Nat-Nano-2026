# FLEX 5′ junction analysis

Python scripts for classifying and quantifying Cre-*lox* recombination outcomes from Nanopore amplicon reads of a FLEX (flip-excision) inversion reporter. The pipeline reads demultiplexed `.fastq` files, distinguishes lox2272-mediated from loxP-mediated inversions at the 5′ junction, and generates the publication figures and a Word summary document.

## What it does

Each amplicon spans a shared forward primer, a constant 58 bp anchor, and a variable junction whose sequence records which *lox* site drove inversion. For every read the pipeline normalizes strand orientation by locating the forward primer, confirms the molecule by fuzzy-matching the anchor, then classifies the post-anchor junction.

Classification is two-step and tolerant of Nanopore error.

1. **lox2272 detection.** A 7 bp motif unique to the lox2272-flip junction (`CGCTGAC`) is searched in a window around the junction start, allowing one mismatch and a few base pairs of anchor-offset slack.
2. **loxP detection.** If the motif is absent, the first 30 bp of the junction are compared to the loxP reference head by edit distance, allowing indels.

Reads that pass the anchor check but match neither reference are counted as unclassified. Length gates separate unflipped molecules (~928 bp) from flipped molecules (~579–584 bp) before classification.

## Repository contents

| Script | Purpose |
|--------|---------|
| `analyze_flex.py` | Core read classifier. Prints per-condition counts of unflipped, loxP-flip, lox2272-flip, and unclassified reads, plus the percentage composition. |
| `plot_flex.py` | Figures 1–5. Total flipped-read counts, per-position junction identity, per-replicate composition, and lox2272 / loxP flip percentages. |
| `plot_alignment.py` | Figure 6. Nucleotide alignment of the most common lox2272-flip junction per condition against the reference. |
| `plot_unique_reads.py` | Figure 7. Per-sample panels of the most abundant unique junction sequences with mismatch highlighting. |
| `plot_flip_pct.py` | Standalone dot-and-bar plot of lox2272-flip percentage from summarized counts. |
| `make_doc_flex.py` | Assembles the figures and tables into a `.docx` analysis summary. |
| `*_no_ubq10.py` | Parallel versions of the plotting scripts that exclude the pUBQ10 condition. Outputs are written to `svg_no_ubq10/`. |

## Requirements

Python 3.9+ with the following packages:

```
matplotlib
numpy
python-docx
```

Install them with:

```bash
python -m venv .venv
source .venv/bin/activate     # Windows PowerShell: .venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

## Input data

The scripts read demultiplexed per-sample `.fastq` files. Point them at your data directory with the `FLEX_FASTQ_DIR` environment variable, otherwise they default to a `fastq/` folder next to the scripts.

```bash
export FLEX_FASTQ_DIR=/path/to/your/fastq        # macOS / Linux
$env:FLEX_FASTQ_DIR = "C:\path\to\your\fastq"    # Windows PowerShell
```

The condition-to-file mapping is defined at the top of each script (the `PAIRS`, `FLIP_FILES`, and `SAMPLES` tables). Edit those entries to match your own sample sheet and file names. Raw `.fastq` files are excluded from version control by `.gitignore`.

## Usage

Run the scripts from the repository directory after setting `FLEX_FASTQ_DIR`.

```bash
# 1. Print the classification tables
python analyze_flex.py

# 2. Generate Figures 1–5
python plot_flex.py

# 3. Generate Figure 6 (junction alignment)
python plot_alignment.py

# 4. Generate Figure 7 (unique junction sequences)
python plot_unique_reads.py

# 5. Assemble the Word summary
python make_doc_flex.py
```

The pUBQ10-excluded variants are run the same way, for example `python plot_flex_no_ubq10.py`.

## Outputs

Each plotting script writes a high-resolution `.png` to the repository directory and a vector `.svg` copy to `svg/` (or `svg_no_ubq10/` for the excluded-condition variants). `make_doc_flex.py` writes a `.docx` that embeds the figures alongside the count tables and methods notes.

## License

MIT
