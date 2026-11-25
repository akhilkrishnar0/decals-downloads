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

    @staticmethod
    def _read_bands(path, kind):
        """Select g/r/z by FITS metadata, not HDU position."""
        wanted = ("g", "r", "z")
        found = {}
        with fits.open(path, memmap=False) as hdul:
            for hdu in hdul:
                if hdu.data is None:
                    continue
                h = hdu.header
                typ = str(h.get("IMAGETYP", h.get("EXTNAME", ""))).strip().upper()
                if kind == "INVVAR" and typ not in ("INVVAR", "INVERSE_VARIANCE"):
                    continue
                if kind == "IMAGE" and typ not in ("", "IMAGE", "PRIMARY"):
                    continue
                data = hdu.data
                if data.ndim == 3:
                    labels = [str(h.get(f"BAND{i}", "")).strip().lower()
                              for i in range(data.shape[0])]
                    if not all(labels):
                        labels = list(str(h.get("BANDS", "")).strip().lower())
                    if len(labels) != data.shape[0]:
                        raise ValueError(f"{path}: cannot determine cube band order")
                    pairs = zip(labels, data)
                elif data.ndim == 2:
                    band = str(h.get("BAND", h.get("FILTER", ""))).strip().lower()
                    pairs = [(band, data)]
                else:
                    continue
                for band, plane in pairs:
                    if band not in wanted:
                        continue
                    if band in found:
                        raise ValueError(f"{path}: multiple {kind} planes for {band}; "
                                         "redownload using the aligned cutout endpoint")
                    found[band] = (plane.copy(), h.copy())
        missing = set(wanted) - set(found)
        if missing:
            raise ValueError(f"{path}: missing {kind} bands {sorted(missing)}; "
                             "inspect FITS headers or redownload with force_download=True")
        return [found[b] for b in wanted]

    @staticmethod
    def _check_grid(images, invvars):
        """Reject mismatched shapes or sky grids; never resample silently."""
        for band, (image, ih), (inv, vh) in zip("grz", images, invvars):
            if image.shape != inv.shape:
                raise ValueError(f"{band}: science {image.shape} vs invvar {inv.shape}; redownload both")
            ny, nx = image.shape
            xx = np.array([0, nx-1, 0, nx-1, (nx-1)/2])
            yy = np.array([0, 0, ny-1, ny-1, (ny-1)/2])
            iw, vw = WCS(ih).celestial, WCS(vh).celestial
            a, d = iw.pixel_to_world_values(xx, yy)
            x, y = vw.world_to_pixel_values(a, d)
            if not (np.allclose(x, xx, atol=1e-3, rtol=0) and
                    np.allclose(y, yy, atol=1e-3, rtol=0)):
                raise ValueError(f"{band}: science and invvar WCS differ; redownload both")

    @staticmethod
    def _split_header(header, band):
        h = WCS(header).celestial.to_header()
        for key in ("BUNIT", "TELESCOP", "INSTRUME", "SURVEY", "VERSION",
                    "MAGZERO", "IMAGETYP"):
            if key in header:
                h[key] = header[key]
        h["BAND"] = band
        return h

    # ---------- LOAD ----------
    def load_data(self, target):
        g = self.clean_name(target)
        images = self._read_bands(self.download_dir / "fits" / f"{g}.fits", "IMAGE")
        invvars = self._read_bands(self.download_dir / "invvar" / f"{g}_invvar.fits", "INVVAR")
        self._check_grid(images, invvars)
        psfs = self._read_bands(self.download_dir / "psf" / f"{g}_psf.fits", "PSF")
        sigmas = []
        for inv, _ in invvars:
            s = np.full(inv.shape, np.nan, dtype=float)
            m = np.isfinite(inv) & (inv > 0)
            s[m] = 1.0 / np.sqrt(inv[m])
            sigmas.append(s)
        wcs = WCS(images[0][1]).celestial
        return [x[0] for x in images], [x[0] for x in invvars], sigmas, [x[0] for x in psfs], wcs

    # ---------- PLOT ----------
    def plot(self, target):
        g = self.clean_name(target)
        bands = ['g','r','z']

        images, invvars, sigmas, psfs, wcs = self.load_data(target)

        invwcs = [wcs] * 3  # grids validated by load_data
        fig = plt.figure(figsize=(12,16))

        for row in range(4):
            for col in range(3):
                ax = plt.subplot(4,3,row*3+col+1,
                                 projection=wcs if row == 0 else invwcs[col] if row < 3 else None)

                if row == 0:
                    data = images[col]
                    vmin,vmax = np.percentile(data[np.isfinite(data)],[21,97]) if np.isfinite(data).any() else (0, 1)
                    ax.imshow(data, origin='lower', vmin=vmin, vmax=vmax)
                    ax.set_title(f"{bands[col]} image")

                elif row == 1:
                    data = invvars[col]
                    vmin,vmax = np.percentile(data[np.isfinite(data)],[1,99]) if np.isfinite(data).any() else (0, 1)
                    ax.imshow(data, origin='lower', vmin=vmin, vmax=vmax)
                    ax.set_title(f"{bands[col]} invvar")

                elif row == 2:
                    data = sigmas[col]
                    vmin,vmax = np.percentile(data[np.isfinite(data)],[1,99]) if np.isfinite(data).any() else (0, 1)
                    ax.imshow(data, origin='lower', vmin=vmin, vmax=vmax)
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
        folder = self.download_dir / "psf"
        planes = self._read_bands(folder / f"{g}_psf.fits", "PSF")
        for band, (data, header) in zip("grz", planes):
            h = header.copy()
            h["BAND"] = band
            fits.writeto(folder / f"{g}_PSF_{band}.fits", data, header=h, overwrite=True)
        print(f"Saved split PSF for {g}")

    def save_split_data(self, target):
        g = self.clean_name(target)
        out_img = self.base_dir / "outputs" / "images"
        out_inv = self.base_dir / "outputs" / "invvar"
        out_img.mkdir(parents=True, exist_ok=True)
        out_inv.mkdir(parents=True, exist_ok=True)
        images = self._read_bands(self.download_dir / "fits" / f"{g}.fits", "IMAGE")
        for band, (data, header) in zip("grz", images):
            fits.writeto(out_img / f"{g}_{band}.fits", data,
                         self._split_header(header, band), overwrite=True)
        invvars = self._read_bands(self.download_dir / "invvar" / f"{g}_invvar.fits", "INVVAR")
        self._check_grid(images, invvars)
        for band, (data, header) in zip("grz", invvars):
            fits.writeto(out_inv / f"{g}_{band}_invvar.fits", data,
                         self._split_header(header, band), overwrite=True)
        print(f"Saved split IMAGE + INVVAR for {g}")

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
