import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u
from decals_query import DECaLSQuery
class DecalsPipeline:
    def __init__(self, base_dir="decals_data", layer="ls-dr9", galaxy_centers=None):
        self.base_dir = Path(base_dir).resolve()
        self.layer = layer
        self.galaxy_centers = {self.clean_name(k): v for k, v in (galaxy_centers or {}).items()}

        # directories
        self.download_dir = self.base_dir / "downloads"
        self.plot_dir = self.base_dir / "plots"

        self.download_dir.mkdir(parents=True, exist_ok=True)
        self.plot_dir.mkdir(parents=True, exist_ok=True)

    def clean_name(self, g):
        return g.replace(" ", "")

    # ---------- COORD ----------

    
    def get_coord(self, target):
        name = self.clean_name(target)
        if name in self.galaxy_centers:
            ra, dec = self.galaxy_centers[name]
            return self._degree_coord(ra, dec)
        try:
            r = Simbad.query_object(target)
            if r is None or len(r) == 0:
                raise ValueError("Not found")
            if "ra" in r.colnames:
                return self._degree_coord(r["ra"][0], r["dec"][0])
            return SkyCoord(str(r["RA"][0]) + " " + str(r["DEC"][0]),
                            unit=(u.hourangle, u.deg), frame="icrs")
        except Exception as exc:
            print(f"SIMBAD lookup failed for {target}: {exc}")
            while True:
                try:
                    ra = input(f"{target} RA (deg): ")
                    dec = input(f"{target} Dec (deg): ")
                except EOFError as err:
                    raise RuntimeError("Pass galaxy_centers to DecalsPipeline.") from err
                try:
                    return self._degree_coord(ra, dec)
                except ValueError as err:
                    print(err)

    @staticmethod
    def _degree_coord(ra, dec):
        ra, dec = float(ra), float(dec)
        if not (np.isfinite(ra) and np.isfinite(dec) and 0 <= ra < 360 and -90 <= dec <= 90):
            raise ValueError("Enter RA in [0, 360) and Dec in [-90, 90] degrees.")
        return SkyCoord(ra, dec, unit="deg", frame="icrs")

    # ---------- CHECK ----------
    def already_downloaded(self, target):
        g = self.clean_name(target)
        f = self.download_dir / "fits" / f"{g}.fits"
        return f.exists()

    # ---------- DOWNLOAD ----------
    def download(self, target, size=1200):
        g = self.clean_name(target)

        coord = self.get_coord(target)
        query = DECaLSQuery(output_dir=str(self.download_dir), layer=self.layer)

        print(f"⬇ Downloading {g}")
        query.query_region(coord, size=size, download_type='all', galid=g)

    # ---------- LOAD ----------
    def load_data(self, target):
        g = self.clean_name(target)

        fits_file   = self.download_dir / "fits" / f"{g}.fits"
        invvar_file = self.download_dir / "invvar" / f"{g}_invvar.fits"
        psf_file    = self.download_dir / "psf" / f"{g}_psf.fits"

        if not fits_file.exists():
            raise FileNotFoundError(f"{fits_file} not found")

        with fits.open(fits_file) as hdul:
            img = hdul[0].data
            images = [img[i, :, :].copy() for i in range(3)]
            wcs = WCS(hdul[0].header).celestial
        with fits.open(invvar_file) as hdul:
            invvars = [hdul[i].data.copy() for i in (2, 4, 6)]
        with fits.open(psf_file) as hdul:
            planes = [h.data.copy() for h in hdul if h.data is not None and h.data.ndim == 2]
            if len(planes) == 3:
                psfs = planes
            else:
                cube = hdul[0].data
                if cube is not None and cube.ndim == 3 and cube.shape[-1] == 3 and cube.shape[0] != 3:
                    psfs = [cube[:, :, i].copy() for i in range(3)]
                elif cube is not None and cube.ndim == 3 and cube.shape[0] == 3 and cube.shape[-1] != 3:
                    psfs = [cube[i].copy() for i in range(3)]
                else:
                    raise ValueError("Unrecognized PSF layout")
        sigmas = []
        for inv in invvars:
            s = np.full(inv.shape, np.nan, dtype=float)
            m = np.isfinite(inv) & (inv > 0)
            s[m] = 1.0 / np.sqrt(inv[m])
            sigmas.append(s)
        return images, invvars, sigmas, psfs, wcs

    # ---------- PLOT ----------
    def plot(self, target):
        g = self.clean_name(target)
        bands = ['g','r','z']

        images, invvars, sigmas, psfs, wcs = self.load_data(target)

        with fits.open(self.download_dir / "invvar" / f"{g}_invvar.fits") as hdul:
            invwcs = [WCS(hdul[i].header).celestial for i in (2, 4, 6)]
        fig = plt.figure(figsize=(12,16))

        for row in range(4):
            for col in range(3):
                ax = plt.subplot(4,3,row*3+col+1,
                                 projection=wcs if row == 0 else invwcs[col] if row < 3 else None)

                if row == 0:
                    data = images[col]
                    vmin,vmax = np.percentile(data[np.isfinite(data)],[21,97]) if np.isfinite(data).any() else (0, 1)
                    ax.imshow(data, vmin=vmin, vmax=vmax)
                    ax.set_title(f"{bands[col]} image")

                elif row == 1:
                    data = invvars[col]
                    vmin,vmax = np.percentile(data[np.isfinite(data)],[1,99]) if np.isfinite(data).any() else (0, 1)
                    ax.imshow(data, vmin=vmin, vmax=vmax)
                    ax.set_title(f"{bands[col]} invvar")

                elif row == 2:
                    data = sigmas[col]
                    vmin,vmax = np.percentile(data[np.isfinite(data)],[1,99]) if np.isfinite(data).any() else (0, 1)
                    ax.imshow(data, vmin=vmin, vmax=vmax)
                    ax.set_title(f"{bands[col]} sigma")

                else:
                    ax.imshow(psfs[col])
                    ax.set_title(f"{bands[col]} PSF")
                    ax.axis("off")
                    continue

                ax.coords.grid(color='white', ls='dotted', alpha=0.5)
                ax.set_xlabel("RA")
                ax.set_ylabel("Dec")

        out = self.plot_dir / f"{g}.png"
        plt.tight_layout()
        plt.savefig(out, dpi=200)
        plt.close()

        print(f"✔ Saved plot {out}")

    # ---------- PSF SPLIT ----------
    def split_psf(self, target):
        g = self.clean_name(target)
        bands = ["g","r","z"]

        psf_file = self.download_dir / "psf" / f"{g}_psf.fits"

        with fits.open(psf_file) as hdul:
            for i, band in enumerate(bands):
                data = hdul[i].data
                header = hdul[i].header

                out = self.download_dir / "psf" / f"{g}_PSF_{band}.fits"

                h = header.copy()
                h["BAND"] = band

                fits.writeto(out, data, header=h, overwrite=True)

        print(f"✔ Split PSF for {g}")

    def save_split_data(self, target):
        g = self.clean_name(target)
        bands = ["g","r","z"]
    
        # create output dirs
        out_img = self.base_dir / "outputs" / "images"
        out_inv = self.base_dir / "outputs" / "invvar"
        out_img.mkdir(parents=True, exist_ok=True)
        out_inv.mkdir(parents=True, exist_ok=True)
    
        fits_file   = self.download_dir / "fits" / f"{g}.fits"
        invvar_file = self.download_dir / "invvar" / f"{g}_invvar.fits"
    
        # ---------- IMAGE ----------
        with fits.open(fits_file) as hdul:
            img = hdul[0].data
            header = hdul[0].header
    
            for i, band in enumerate(bands):
                data = img[i,:,:]
    
                h = WCS(header).celestial.to_header()
                for key in ('BUNIT', 'TELESCOP', 'INSTRUME', 'SURVEY', 'MAGZERO'):
                    if key in header:
                        h[key] = header[key]
                h["BAND"] = band
    
                out = out_img / f"{g}_{band}.fits"
                fits.writeto(out, data, header=h, overwrite=True)
    
        # ---------- INVVAR ----------
        with fits.open(invvar_file) as hdul:
            invs = [hdul[2].data, hdul[4].data, hdul[6].data]
    
            for i, band in enumerate(bands):
                data = invs[i]
    
                h = hdul[(2, 4, 6)[i]].header.copy()
                h["BAND"] = band
    
                out = out_inv / f"{g}_{band}_invvar.fits"
                fits.writeto(out, data, header=h, overwrite=True)
    
        print(f"✔ Saved split IMAGE + INVVAR for {g}")

    

    def depth_maps(self, target, region="south"):
        """Get DR9 PSF depth and PSF flux sigma on the downloaded science grid."""
        if self.layer not in ("ls-dr9", "ls-dr9-south", "ls-dr9-north"):
            raise ValueError("Brick depth option supports DR9 only")
        if self.layer.endswith("-north") and region != "north" or self.layer.endswith("-south") and region != "south":
            raise ValueError("Depth region must match the science layer")
        from decals_depth import make_depth_maps
        g = self.clean_name(target)
        return make_depth_maps(
            self.download_dir / "fits" / f"{g}.fits",
            self.base_dir / "brick_maps",
            self.base_dir / "outputs" / "depth" / g,
            self.plot_dir / f"{g}_depth.png", region=region)

    # ---------- RUN ----------
    def run(self, targets, size=1200,
            skip_download=True,
            force_download=False,
            do_plot=True,
            do_psf_split=True,
            do_save_split=True, do_depth=False, depth_region='south'):

        if isinstance(targets, str):
            targets = [targets]

        for t in targets:
            print(f"\n=== {t} ===")

            # DOWNLOAD CONTROL (fixed)
            if force_download:
                print(f"🔁 Force downloading {t}")
                self.download(t, size=size)

            elif skip_download and self.already_downloaded(t):
                print(f"⏭ Skipping download for {t}")

            else:
                self.download(t, size=size)

            # PLOT
            if do_plot:
                try:
                    self.plot(t)
                except Exception as e:
                    print(f"❌ Plot failed: {e}")
            # SAVE SPLIT DATA
            if do_save_split:
                try:
                    self.save_split_data(t)
                except Exception as e:
                    print(f"❌ Save split failed: {e}")

            # PSF SPLIT
            if do_psf_split:
                try:
                    self.split_psf(t)
                except Exception as e:
                    print(f"❌ PSF split failed: {e}")

            if do_depth:
                self.depth_maps(t, region=depth_region)
