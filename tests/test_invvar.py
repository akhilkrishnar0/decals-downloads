import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from decals_pipeline_class import DecalsPipeline


def test_cutout_cube(tmp_path):
    w = WCS(naxis=2)
    w.wcs.crpix = [3, 3]
    w.wcs.crval = [65, -55]
    w.wcs.cdelt = [-0.262/3600, 0.262/3600]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    h = w.to_header()
    h["BANDS"] = "griz"
    for i, b in enumerate("griz"):
        h[f"BAND{i}"] = b
    cube = np.stack([np.full((5, 5), i, dtype=np.float32) for i in (1, 2, 3, 4)])
    h["IMAGETYP"] = "IMAGE"
    image = fits.PrimaryHDU(cube, h)
    h["IMAGETYP"] = "INVVAR"
    path = tmp_path / "cutout.fits"
    fits.HDUList([image, fits.ImageHDU(cube*10, h)]).writeto(path)
    images = DecalsPipeline._read_bands(path, "IMAGE")
    inv = DecalsPipeline._read_bands(path, "INVVAR")
    assert [x[0][0,0] for x in images] == [1,2,4]
    assert [x[0][0,0] for x in inv] == [10,20,40]
    DecalsPipeline._check_grid(images, inv)
