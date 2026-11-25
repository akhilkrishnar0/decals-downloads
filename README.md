# DECaLS-download

Download Legacy Surveys g/r/z images, inverse-variance maps and PSFs.

## Install

```bash
pip install --upgrade git+https://github.com/akhilkrishnar0/decals-pipeline.git
```

## Example

```python
from decals_pipeline import DecalsPipeline

pipe = DecalsPipeline()
pipe.run(['NGC6902'], size=3000)
```

Outputs are saved in `decals_data/`: original downloads, split science FITS, PSFs and diagnostic plots. Size is in pixels (maximum 3000; server limits may differ). Existing downloads are reused when valid; use `force_download=True` to refresh.

## Manual coordinates (if SIMBAD is unavailable)

If SIMBAD fails, the pipeline automatically asks for RA and Dec in decimal degrees. Alternatively, bypass SIMBAD entirely for a single target:

```python
pipe.run(['NGC6902'], size=3000, ra=306.117375, dec=-43.653681)
```

Coordinates are ICRS decimal degrees (RA: 0–360, Dec: −90–90). Verify the centre before downloading. For multiple galaxies, omit ra/dec to resolve or enter each centre separately.

## Plots

- `plots/NGC6902/grz.png`: science-only review with galaxy centre, independent of inverse variance.
- `plots/NGC6902.png`: 12-panel diagnostic when the native inverse-variance file has a supported single-brick layout.
- `plots/NGC6902/native_invvar/`: separate plots of explicitly labelled native inverse-variance HDUs when the full diagnostic is unavailable.

Science FITS are split before plotting. Native maps retain their own WCS and are not mosaicked or treated as aligned uncertainties. Unsupported diagnostic layouts issue warnings; the science preview and split images remain available. Set `do_plot=False` to skip previews. Set `include_invvar=True` only for strict WCS-checked inverse-variance splitting; mismatches raise errors.
