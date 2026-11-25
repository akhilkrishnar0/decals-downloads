# DECaLS Downloads

Download and split Legacy Surveys g/r/z science images, native inverse-variance maps and PSFs.

## Installation

```bash
pip install --upgrade git+https://github.com/akhilkrishnar0/decals-downloads.git
```

## Example: NGC6902

```python
from decals_pipeline_class import DecalsPipeline

galaxies_new = ["NGC6902"]
GALAXY_CENTERS = {"NGC6902": (306.117, -43.653)}  # RA, Dec in degrees

pipe = DecalsPipeline(layer="ls-dr9", galaxy_centers=GALAXY_CENTERS)
pipe.run(
    galaxies_new, size=3000, skip_download=True,
    do_plot=False, do_psf_split=False,
    force_download=False, do_save_split=True
)
```

Coordinates are supplied by the user. If omitted, SIMBAD is queried; on failure, enter RA and Dec manually in degrees.

Set `do_plot=True` for the 12-panel diagnostic. Outputs are saved in `decals_data/`.

Native inverse-variance maps use the original HDU 2/4/6 selection and may have a different WCS and footprint from the science images. Verify HDU band and brick identity before scientific use; these are not aligned uncertainty maps.
