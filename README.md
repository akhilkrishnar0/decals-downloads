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

Outputs are saved in `decals_data/`: original downloads, split science FITS, PSFs and `plots/NGC6902.png`. Size is in pixels (maximum 3000; server limits may differ). Existing downloads are reused when valid; use `force_download=True` to refresh.

## Diagnostic plot

Four rows show g/r/z science images, native inverse variance, sigma (1/sqrt(invvar)), and PSFs. Each map uses its own WCS. Native inverse variance can cover a different sky region; it is not an aligned uncertainty map. The full plot requires a single-brick seven-HDU inverse-variance product. For other layouts, use `do_plot=False` and inspect the native products separately.

```python
pipe.run(['NGC6902'], size=3000, do_plot=False, include_invvar=False)
```

Set `include_invvar=True` only to export inverse variance after strict image-grid checks. Invalid weights become NaN in sigma plots.
