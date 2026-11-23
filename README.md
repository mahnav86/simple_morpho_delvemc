# simple_morpho_delvemc

Search for ultra-faint Milky Way satellite (dwarf galaxy) candidates in DELVE-MC DR3, using a
simple binned hotspot search with morphocloud star/galaxy separation.

This is the `simple_adl_morpho` package, derived from `simple_adl` (Sidney Mau, MIT licence).
It is a work in progress.

## What is in this repository

| File | Purpose |
|---|---|
| `query_morpho.py` | Cone search (default radius 1.9 deg) on DELVE-MC DR3, with morphocloud star/galaxy separation. Used by `survey.py`. |
| `query_healpix_morpho.py` | HEALPix version: selects objects by `ring256` for a set of nside=32 pixels, no calibration. |
| `morpho_separation.py` | TAP fetch and morphocloud classifier used by the queries above. |
| `query_dl.py`, `query_healpix.py`, `query_healpix_nei.py` | Earlier Data Lab queries (SQL star/galaxy cuts). |
| `survey.py`, `search.py` | Survey setup and the hotspot search for one position. |
| `parallel_search.py` | Runs `search.py` over many positions in parallel. |
| `make_list.py` | Builds `candidate_list.csv` from the results. |
| `plot_hotspot.py`, `parallel_plot_hotspot.py` | Diagnostic plots for candidates. |
| `isochrone.py`, `projector.py`, `coordinate_tools.py`, `associations.py` | Supporting code. |
| `config.yaml` | Survey, catalogue and output settings. |
| `search_coordinates.csv`, `search_coordinates_healpix.csv` | Positions to search. |

## Installation

1. Clone the repository and put its **parent folder** on `PYTHONPATH`, because the code imports
   itself as `simple_adl_morpho.<module>`:

   ```bash
   git clone https://github.com/mahnav86/simple_morpho_delvemc.git
   export PYTHONPATH=$PWD:$PYTHONPATH
   ```

   The folder that contains `simple_adl_morpho/` must be the one on `PYTHONPATH`.

2. Create a Python environment with:
   `numpy`, `scipy`, `pandas`, `astropy`, `healpy`, `pyyaml`, `fitsio`, `matplotlib`, `pyvo`,
   and `morphocloud` (the star/galaxy classifier).

3. The Data Lab TAP service (https://datalab.noirlab.edu/tap) must be reachable. The queries
   retry three times if the server returns an error.

## Configuration

`config.yaml` sets:

- `catalog.mag_max`: the magnitude limit. It is passed to the r-band cut in the queries and to
  the star/galaxy filters in `survey.py`.
- `isochrone`: age and metallicity of the fitted isochrone.
- `output`: where candidate lists, logs, results and plots are written.

## Single search
Activate the conda environment with

```bash
 conda activate morphocloud
```

Search one position by coordinates, or by HEALPix pixel (nside 32):

```bash
python simple_adl_morpho/search.py --config simple_adl_morpho/config.yaml --outfile out.csv --ra 56.36 --dec -60.45
python simple_adl_morpho/search.py --config simple_adl_morpho/config.yaml --outfile out.csv --nside 32 --ipix 11460
```

### Choosing the query mode

`search.py` and `plot_hotspot.py` take `--query cone|healpix`:

- `cone` (default): queries a 1.9 deg circle around the search position by coordinates
  (`query_morpho`).
- `healpix`: queries the `ring256` pixels covering the nside=32 search pixel
  (`query_healpix_morpho`). Use with `--nside 32 --ipix <pixel>`.

The default is `cone`, because `config.yaml` has no `query` key. `--query` overrides it for that
run only. Use the same mode for a search and its plot.

```bash
# cone (default)
python simple_adl_morpho/search.py --config simple_adl_morpho/config.yaml --outfile out.csv --nside 32 --ipix 11812
# healpix
python simple_adl_morpho/search.py --config simple_adl_morpho/config.yaml --outfile out.csv --nside 32 --ipix 11812 --query healpix
```

The two modes return different catalogues, so the peaks differ. State the mode when you compare
or report results.

Plot one candidate:

```bash
python simple_adl_morpho/plot_hotspot.py --ra 56.36 --dec -60.45 --mod 19.5
```

## Parallel search procedure

1. Add target coordinates to `search_coordinates.csv`.
2. In `parallel_search.py`, set the number of threads and the file paths for your machine.
   On Linux, check the number of threads per core with `lscpu`. Try a short test run with fewer
   threads first.
3. Run `python simple_adl_morpho/parallel_search.py`. On a cluster, use `nohup` so the job keeps
   running after you disconnect, then follow it with `tail -f nohup.out`.
4. Check per-search logs in `log_dir/`.
5. Run `python simple_adl_morpho/make_list.py` to create `candidate_list.csv`, and review it.
6. In `parallel_plot_hotspot.py`, set the thread count and paths, then run it to make plots.

## Checking the output

- `grep -L "Job done" *` in the log directory lists incomplete runs.
- If a log has no matching `.csv` in the results, delete any partial output and rerun that
  search without `nohup`, so you can watch it.
- If grep finds "Job done" but the log looks incomplete, open the log in a different editor, since
  some editors do not load the whole file.

## Generated files

Logs, results, plots, candidate lists and `nohup.out` are generated and are listed in
`.gitignore`. They are not part of the repository.

## Known issues

- `parallel_search.py` still has a hard-coded morphocloud Python path and log paths for one
  machine. Check them before a run.
- Some known objects are missed by the search. See the findings notes in the project folder.

## Licence

To be decided. The original `simple_adl` code is MIT-licensed (Sidney Mau). Add a licence file
before making the repository public.
