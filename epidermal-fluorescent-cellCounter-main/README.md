# lif-microscopy-analysis

Interactive Jupyter workflow for analyzing Leica confocal `.lif` microscopy files. It combines manual nuclear click-counting with parallelized automated cytosol segmentation via [Cellpose](https://github.com/MouseLand/cellpose), and exports publication-ready QC images and a final CSV.

## What it does

The notebook opens an interactive viewer that steps through every image series in one `.lif` file or a whole folder of them. You click on nuclei to count them by hand while automated cytosol segmentation runs in the background on worker threads, so the interface never blocks. Segmentation results appear as overlays when they finish, QC images are written to a structured output tree, and every navigation event is checkpointed to JSON so a crashed or closed session can be resumed.

## Analysis modes

| Mode | Nuclear counting | Cytosol counting |
|------|-----------------|-----------------|
| **Protoplast** | Manual click (yellow) | Automated — Cellpose `cyto3` |
| **Epidermal** | Manual click (yellow) | Automated — Cellpose `cyto3` (large diameter) + watershed fallback |
| **Manual** | Manual click (yellow) | Manual click (magenta) |

Pseudocolor is fixed throughout the viewer and the saved QC images: nucleus channels render **yellow**, cytosol channels render **magenta**.

## Requirements

Python 3.10+ with the following packages:

```
readlif
cellpose
numpy
scipy
scikit-image
opencv-python-headless
matplotlib
ipympl
ipywidgets
pandas
Pillow
tqdm
jupyterlab
```

Install all dependencies with:

```bash
pip install -r requirements.txt
```

A dedicated conda environment is recommended:

```bash
conda create -n lifanalysis python=3.10
conda activate lifanalysis
pip install -r requirements.txt
```

> **Note:** The first time Cellpose runs a given model (e.g. `cyto3`) it downloads ~280 MB of model weights and caches them in `~/.cellpose/models/`, so the first run needs internet connectivity. The background-threading architecture (a `ThreadPoolExecutor`, not multiple processes) was chosen deliberately, since Cellpose and scikit-image both release the GIL and process pools cannot serialize matplotlib figures or open `LifFile` handles across process boundaries in notebooks.

### GPU acceleration (optional)

Install PyTorch with CUDA support before Cellpose, then set `"gpu": True` in the Cellpose config blocks inside `CONFIG`:

```bash
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu118
pip install cellpose
```

## Usage

1. Open `lif_microscopy_analysis.ipynb` in Jupyter
2. Edit the **CONFIG cell** to set `input_path`, `output_dir`, and `default_mode`
3. Run all cells top to bottom (`Kernel → Restart & Run All`) — the interactive viewer launches automatically in the **Launch Workflow** cell
4. Click on nuclei to count them; automated cytosol jobs run in the background
5. When every image is counted and all background jobs finish, run the **Final CSV Export** cell

### Key configuration options

| Parameter | Description |
|-----------|-------------|
| `input_path` | Path to a single `.lif` file **or** a folder of `.lif` files |
| `output_dir` | Directory where results, QC images, masks, and logs are written |
| `default_mode` | `"protoplast"`, `"epidermal"`, or `"manual"` |
| `default_nucleus_channel` / `default_cytosol_channel` | 0-based channel indices (overrideable per image in the GUI) |
| `z_projection` | Projection for z-stacks — `"max"`, `"mean"`, or `"first"` |
| `cellpose_protoplast` / `cellpose_epidermal` | Cellpose parameters (diameter, thresholds, GPU, min size) per mode |
| `watershed_epidermal` | Fallback watershed parameters used when Cellpose is unavailable or returns 0 cells |
| `max_workers` | Number of background segmentation threads |
| `save_masks` / `save_overlays` | Toggle saving of label masks and QC overlays |

## Notebook structure

The notebook is organized as a small class-based application:

- `LifLoader` — discovers `.lif` files, enumerates series, extracts dimensions and scale
- `LIFWorkflowApp` — main widget controller that steps through items and loads channels on demand
- `BackgroundJobQueue` — `ThreadPoolExecutor`-backed queue for automated segmentation
- `run_cellpose` / `run_epidermal_segmentation` — the two automated counting paths
- `QCExporter` — saves pseudocolored PNGs with count overlays
- `CheckpointManager` — atomic JSON checkpoint written after every navigation event, enabling resume

Utility cells at the end provide a background job-status monitor, a resume-from-checkpoint path, a no-GUI smoke test, and a dataset overview table.

## Output files

Results are written to `output_dir` in a structured tree, including per-series counts, saved label masks, pseudocolored QC overlays, run logs, and a final aggregated CSV combining all counted series. Raw `.lif` files and generated results are excluded from version control via `.gitignore`.

## Limitations

- Only the first timepoint (T=0) of a series is processed; multi-timepoint time-lapse analysis is not supported.
- The viewer displays one channel at a time; composite RGB views are not implemented.
- Cellpose GPU acceleration is off by default.
- Channels must be assigned per-item when they vary across files.

## License

MIT
