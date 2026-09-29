# nuclei-counter_cellPose

Interactive click-to-count tool for manual nuclei counting in Leica `.lif` microscopy files, with automated protoplast counting via [Cellpose](https://github.com/MouseLand/cellpose) running in a background thread.

## What it does

The tool opens an interactive GUI displaying each image series from your `.lif` file(s). You click on nuclei to count them manually, while Cellpose automatically counts protoplasts in a separate channel in the background. Results are exported as CSV and Excel files with per-series counts and nuclei-per-cell ratios.

## Requirements

Python 3.10+ with the following packages:

```
readlif
cellpose
matplotlib
numpy
pandas
Pillow
openpyxl
```

Install all dependencies with:

```bash
pip install -r requirements.txt
```

> **Note:** This tool is compatible with Cellpose v4 (CPSAM). The background threading architecture was designed specifically to avoid SEGFAULT issues that arise from spawning multiple PyTorch processes on macOS/Apple Silicon when using the CPSAM model.

## Usage

1. Open `manual_nuclei_counter.ipynb` in Jupyter
2. Edit the **CONFIG cell** (Cell 1) to point to your `.lif` files and set your channel indices
3. Run all cells (`Kernel → Restart & Run All`) — the counter window opens automatically

### Configuration options

| Parameter | Description |
|-----------|-------------|
| `LIF_FOLDER` | Path to a folder — all `.lif` files inside will be loaded |
| `LIF_PATHS` | Explicit list of `.lif` files (used when `LIF_FOLDER = None`) |
| `OUTPUT_DIR` | Directory where results will be saved |
| `CHANNEL_NUCLEI` | Channel index (0-based) to display for manual counting |
| `CHANNEL_CYTO` | Channel index (0-based) for automated Cellpose protoplast counting |
| `USE_MAX_PROJECTION` | If `True`, max-projects z-stacks before display |
| `RUN_CELLPOSE` | Set `False` to skip automated counting entirely |
| `SAVE_QC` | If `True`, saves side-by-side raw/mask PNG overlays per series |

### GUI controls

| Action | How |
|--------|-----|
| Count a nucleus | Left-click on it |
| Remove last click | `Undo` button |
| Clear current series | `Reset` button |
| Go to previous/next series | `< Prev` / `Next >` buttons |
| Export CSV + Excel | `Save` button (can click multiple times) |
| Adjust brightness/contrast | Drag the **Black pt** / **White pt** sliders |

The title bar shows live counts: manual nuclei, Cellpose protoplasts, and nuclei-per-cell ratio. Cellpose runs in a background thread and updates the title bar when it finishes — no need to wait before clicking.

## Output files

Results are saved to `OUTPUT_DIR` with the following structure:

```
OUTPUT_DIR/
├── <lif_stem>/
│   ├── manual_counts.csv       # per-series counts for this LIF file
│   ├── manual_counts.xlsx      # same data + summary stats sheet
│   └── cellpose_qc/            # side-by-side raw/mask PNGs (if SAVE_QC = True)
├── all_counts.csv              # combined results across all LIF files
└── all_counts.xlsx             # combined results + per-LIF summary sheet
```

Each row in the output contains: LIF filename, series name, automated protoplast count, manual nuclei count, nuclei-per-cell ratio, and raw click coordinates.

## License

MIT
