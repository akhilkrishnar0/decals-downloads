# DECaLS Downloads

Download and split Legacy Surveys g/r/z science images, native inverse-variance maps and PSFs.

## Installation

```bash
pip install --upgrade git+https://github.com/akhilkrishnar0/decals-downloads.git
```

## Example: UGC9024

```python
from decals_pipeline_class import DecalsPipeline

galaxies_new = ["UGC9024"]
GALAXY_CENTERS = {"UGC9024": (211.669,  22.070)}  # RA, Dec in degrees

pipe = DecalsPipeline(layer="ls-dr9", galaxy_centers=GALAXY_CENTERS)
pipe.run(
    galaxies_new, size=1500, skip_download=True,
    do_plot=False, do_psf_split=False,
    force_download=False, do_save_split=True
)
```

Coordinates are supplied by the user. If omitted, SIMBAD is queried; on failure, enter RA and Dec manually in degrees.

Set `do_plot=True` for the 12-panel diagnostic. Outputs are saved in `decals_data/`.

Native inverse-variance maps use the original HDU 2/4/6 selection and may have a different WCS and footprint from the science images. Verify HDU band and brick identity before scientific use; these are not aligned uncertainty maps.
## Example output — NGC6902

![UGC9024 DECaLS diagnostic plot](UGC9024(1).png)
