# DECaLS Downloads

Download and split Legacy Surveys g/r/z science images, native inverse-variance maps and PSFs.

Data acquisition and diagnostic visualisation only.

[Installation](#installation) · [Image download](#example-ugc9024) · [Brick depth](#optional-dr9-brick-depth--ngc1566)

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

## Example output — UGC9024

![UGC9024 DECaLS diagnostic plot](figures/UGC9024.png)

## Why use brick-depth maps?

In some cases, the inverse-variance maps downloaded from DECaLS do not match the dimensions or spatial coverage of the corresponding science images because the requested cutout can extend across multiple survey bricks. This prevents a direct pixel-by-pixel estimation of photometric uncertainties. To address this, we retrieve the DR9 brick-level depth maps covering the entire science image. The Legacy Surveys divide the sky into bricks of approximately 0.25° × 0.25°, each containing coadded science images, inverse-variance maps, PSF information and depth maps. The depth map represents the spatially varying inverse variance of the PSF-fitted point-source flux, accounting for variations in observational coverage, background noise and image quality. By sampling the relevant brick-depth maps onto the science image's WCS, we obtain spatially consistent sensitivity maps in the g, r and z bands. The corresponding point-source flux uncertainty is given by $\sigma_{\rm PSF}=1/\sqrt{\mathrm{depth}}$, in nanomaggies. However, these depth maps are not equivalent to pixel-level inverse-variance maps, and additional assumptions are required to estimate uncertainties for extended emission or aperture-based photometry.

**References**

- [Legacy Surveys DR9 — Data Release Description](https://www.legacysurvey.org/dr9/description/)
- [Legacy Surveys DR9 — Data Products and Depth Map Definitions](https://www.legacysurvey.org/dr9/files/)
- [DR9 Brick Data Repository — NERSC](https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr9/)
- [Discussion on Bricks, Inverse-Variance Maps and Cutout Shapes](https://discuss.legacysurvey.org/t/question-about-bricks-inverse-variance-maps-and-cutout-shapes/10502)

The last reference also discusses the `&invvar` cutout option at the native pixel scale of 0.262 arcsec/pixel, which is worth investigating as a direct solution for obtaining matching pixel-level inverse-variance maps.

## Example output — NGC1566

![NGC1566 DECaLS diagnostic plot](figures/NGC1566.png)

## Optional: DR9 brick depth — NGC1566

If native inverse-variance maps have a different footprint, download DR9
PSF-depth bricks and sample them on the science image grid:

```python
from decals_pipeline_class import DecalsPipeline

pipe = DecalsPipeline(
    layer="ls-dr9",
    galaxy_centers={"NGC1566": (65.0042, -54.9378)}
)
pipe.run(["NGC1566"], size=3000, skip_download=True,
         do_plot=True, do_psf_split=False,
         do_save_split=True, do_depth=True, depth_region="south")
```

For previously downloaded science images only:

```python
pipe.depth_maps("NGC1566", region="south")
```

- `outputs/depth/NGC1566/`: g/r/z depth and PSF flux sigma FITS on the science WCS and shape; brick coverage CSV.
- `plots/NGC1566_depth.png`: science, depth and PSF flux sigma (3 × 3).
- `brick_maps/`: cached original DR9 FITS.

All bricks intersecting the science pixel centres are considered. Each centre
uses its unique survey brick, sampled by nearest neighbour. Missing or
nonpositive values remain NaN. Choose the matching DR9 north/south region
explicitly; NGC1566 uses south.

**Depth is PSF-fitted point-source flux inverse variance**, not science-pixel
inverse variance. Sigma = 1/sqrt(depth), in nanomaggies. No aperture-noise
conversion or variance resampling is performed. Original inverse-variance
files are unchanged.
### Example depth output — NGC1566

![NGC1566 DR9 brick-depth diagnostic plot](figures/NGC1566_depth.png)

## Repository structure

```text
src/        Python modules
examples/   Runnable download and brick-depth examples
figures/    Example diagnostic plots
```
