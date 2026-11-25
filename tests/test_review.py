import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from decals_pipeline import DecalsPipeline

class ReviewTest(unittest.TestCase):
    def test_multihdu_diagnostic_does_not_block_science(self):
        with tempfile.TemporaryDirectory() as directory:
            pipe = DecalsPipeline(base_dir=directory)
            g = "NGC6902"
            for kind in ("fits", "invvar"):
                (pipe.download_dir / kind).mkdir()
            w = WCS(naxis=2)
            w.wcs.crpix = [4, 4]
            w.wcs.cdelt = [-0.0001, 0.0001]
            w.wcs.crval = [306.117375, -43.653681]
            w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
            header = w.to_header()
            fits.writeto(pipe.download_dir / "fits" / f"{g}.fits",
                         np.ones((3, 8, 8)), header)
            h = fits.ImageHDU(np.ones((8, 8)), header=header, name="INVVAR")
            fits.HDUList([fits.PrimaryHDU(), h]).writeto(
                pipe.download_dir / "invvar" / f"{g}_invvar.fits")
            with patch.object(pipe, "download", return_value=SkyCoord(306.117375, -43.653681, unit="deg")):
                with self.assertWarns(RuntimeWarning):
                    pipe.run([g], do_psf_split=False)
            self.assertTrue((pipe.plot_dir / g / "grz.png").is_file())
            self.assertTrue((pipe.plot_dir / g / "native_invvar" / "hdu_1.png").is_file())
            for band in "grz":
                self.assertTrue((pipe.base_dir / "outputs" / "images" / f"{g}_{band}.fits").is_file())

if __name__ == "__main__":
    unittest.main()
