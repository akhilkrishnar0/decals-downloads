import os
import warnings
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u



class DecalsPipeline:
    def __init__(self, base_dir="decals_data", layer="ls-dr9"):
        self.base_dir = Path(base_dir).resolve()
        self.layer = layer

        # directories
        self.download_dir = self.base_dir / "downloads"
        self.plot_dir = self.base_dir / "plots"

        self.download_dir.mkdir(parents=True, exist_ok=True)
        self.plot_dir.mkdir(parents=True, exist_ok=True)

    def clean_name(self, g):
        return g.replace(" ", "")

    # ---------- COORD ----------
    @staticmethod
    def _degree_coord(ra, dec):
        ra, dec = float(ra), float(dec)
        if not (np.isfinite(ra) and np.isfinite(dec) and 0 <= ra < 360 and -90 <= dec <= 90):
            raise ValueError("Use decimal degrees: 0 <= RA < 360 and -90 <= Dec <= 90.")
        return SkyCoord(ra, dec, unit=u.deg, frame="icrs")

    def get_coord(self, target, ra=None, dec=None):
        if ra is not None or dec is not None:
            if ra is None or dec is None:
                raise ValueError("Supply both ra and dec in decimal degrees.")
            return self._degree_coord(ra, dec)
        try:
            result = Simbad.query_object(target)
            if result is None or len(result) == 0:
                raise ValueError(f"{target} not found in SIMBAD")
            if 'ra' in result.colnames and 'dec' in result.colnames:
                return self._degree_coord(result['ra'][0], result['dec'][0])
            return SkyCoord(result['RA'][0] + " " + result['DEC'][0],
                            unit=(u.hourangle, u.deg), frame="icrs")
        except Exception as exc:
            print(f"SIMBAD lookup failed for {target}: {exc}")
            print("Enter ICRS coordinates manually in decimal degrees.")
            while True:
                try:
                    ra = input(f"{target} RA (deg): ")
                    dec = input(f"{target} Dec (deg): ")
                except EOFError as err:
                    raise RuntimeError("No interactive input available. Pass ra= and dec= to run().") from err
                try:
                    return self._degree_coord(ra, dec)
                except ValueError as err:
                    print(f"Invalid coordinates: {err}")

    # ---------- CHECK ----------
    def already_downloaded(self, target):
        g = self.clean_name(target)
        f = self.download_dir / "fits" / f"{g}.fits"
        return f.exists()

    # ---------- DOWNLOAD ----------
    def download(self, target, size=1200, overwrite=True, ra=None, dec=None):
        from .query import DECaLSQuery
        g = self.clean_name(target)

        coord = self.get_coord(target) if ra is None and dec is None else self.get_coord(target, ra=ra, dec=dec)
        query = DECaLSQuery(output_dir=str(self.download_dir), layer=self.layer, overwrite=overwrite)

        print(f"⬇ Downloading {g}")
        query.query_region(coord, size=size, download_type='all', galid=g)
        return coord

    def _check_image_invvar(self, target):
        name = self.clean_name(target)
        with fits.open(self.download_dir / 'fits' / f'{name}.fits') as images, fits.open(self.download_dir / 'invvar' / f'{name}_invvar.fits') as inv:
            cube = images[0].data
            if cube is None or cube.ndim != 3 or cube.shape[0] != 3:
                raise ValueError('Expected primary image cube (3, ny, nx), ordered g/r/z')
            if len(inv) != 7:
                raise ValueError('Expected one image/INVVAR pair per g/r/z band (7 HDUs). Multiple-brick subimages require explicit selection or mosaicking.')
            wcs = WCS(images[0].header).celestial
            yy, xx = np.meshgrid(np.linspace(0, cube.shape[1]-1, 3), np.linspace(0, cube.shape[2]-1, 3), indexing='ij')
            for band, index in zip(('g', 'r', 'z'), (2, 4, 6)):
                hdu = inv[index]
                label = str(hdu.header.get('BAND', hdu.header.get('FILTER', band))).strip().lower()
                if label != band:
                    raise ValueError(f'Unexpected band in inverse-variance HDU {index}: {label}')
                if hdu.data is None or hdu.data.shape != cube.shape[1:]:
                    raise ValueError(f'{band} image/inverse-variance shapes differ; native subimages cannot be used as aligned resampled maps')
                ix, iy = WCS(hdu.header).celestial.world_to_pixel(wcs.pixel_to_world(xx, yy))
                if not (np.allclose(ix, xx, rtol=0, atol=1e-3) and np.allclose(iy, yy, rtol=0, atol=1e-3)):
                    raise ValueError(f'{band} image/inverse-variance WCS mismatch. Use matched native products or a documented variance-resampling procedure.')

    # ---------- LOAD ----------
    def load_data(self, target):
        self._check_image_invvar(target)
        g = self.clean_name(target)

        fits_file   = self.download_dir / "fits" / f"{g}.fits"
        invvar_file = self.download_dir / "invvar" / f"{g}_invvar.fits"
        psf_file    = self.download_dir / "psf" / f"{g}_psf.fits"

        if not fits_file.exists():
            raise FileNotFoundError(f"{fits_file} not found")

        with fits.open(fits_file) as hdul:
            img = hdul[0].data
            if img is None or img.ndim != 3 or img.shape[0] != 3:
                raise ValueError("Expected a primary image cube shaped (3, ny, nx), ordered g/r/z")
            images = [img[i].copy() for i in range(3)]
            wcs = WCS(hdul[0].header).celestial
        with fits.open(invvar_file) as hdul:
            invvars = [hdul[i].data.copy() for i in (2, 4, 6)]
        if any(a.shape != b.shape for a, b in zip(images, invvars)):
            raise ValueError("Image and inverse-variance shapes differ")
        psfs, _ = self._read_psfs(psf_file)

        # SIGMA
        sigmas = []
        for inv in invvars:
            s = np.full(inv.shape, np.nan, dtype=float)
            m = np.isfinite(inv) & (inv > 0)
            s[m] = 1.0 / np.sqrt(inv[m])
            sigmas.append(s)

        return images, invvars, sigmas, psfs, wcs

    # ---------- PLOT ----------
    def review(self, target, ra, dec):
        """Science-only g/r/z preview; independent of inverse variance and PSF."""
        from .review import prepare_decals_review
        coord = self._degree_coord(ra, dec)
        g = self.clean_name(target)
        return prepare_decals_review(self.download_dir, self.plot_dir / g, g,
                                     coord.ra.deg, coord.dec.deg)

    def plot_native_invvar(self, target):
        """Plot explicitly identified native inverse-variance HDUs separately."""
        g = self.clean_name(target)
        outdir = self.plot_dir / g / "native_invvar"
        outdir.mkdir(parents=True, exist_ok=True)
        paths = []
        with fits.open(self.download_dir / "invvar" / f"{g}_invvar.fits") as hdul:
            for idx, h in enumerate(hdul):
                label = " ".join(str(h.header.get(k, "")) for k in
                                 ("EXTNAME", "IMTYPE", "PRODTYPE")).lower()
                if not any(x in label for x in ("invvar", "inverse variance", "inverse_variance")):
                    continue
                if h.data is None or h.data.ndim != 2:
                    continue
                wcs = WCS(h.header).celestial
                if not wcs.has_celestial:
                    continue
                inv = h.data
                sigma = np.full(inv.shape, np.nan)
                good = np.isfinite(inv) & (inv > 0)
                sigma[good] = 1 / np.sqrt(inv[good])
                band = h.header.get("BAND", h.header.get("FILTER", "?"))
                fig = plt.figure(figsize=(10, 5))
                try:
                    for j, (name, data) in enumerate((("invvar", inv), ("sigma", sigma)), 1):
                        ax = fig.add_subplot(1, 2, j, projection=wcs)
                        finite = data[np.isfinite(data)]
                        lo, hi = np.percentile(finite, (1, 99)) if finite.size else (0, 1)
                        ax.imshow(data, origin="lower", cmap="viridis", vmin=lo, vmax=hi)
                        ax.set_title(f"HDU {idx} | {band} {name}")
                        ax.set_xlabel("RA")
                        ax.set_ylabel("Dec")
                    fig.tight_layout()
                    path = outdir / f"hdu_{idx}.png"
                    fig.savefig(path, dpi=150)
                    paths.append(path)
                finally:
                    plt.close(fig)
        if not paths:
            warnings.warn("No explicitly labelled native inverse-variance HDUs found; original FITS retained.", RuntimeWarning)
        return paths

    def plot(self, target):
        g = self.clean_name(target)
        bands = ['g','r','z']

        with fits.open(self.download_dir / "fits" / f"{g}.fits") as hdul:
            cube = hdul[0].data
            if cube is None or cube.ndim != 3 or cube.shape[0] != 3:
                raise ValueError("Expected g/r/z science cube")
            images = [a.copy() for a in cube]
            wcs = WCS(hdul[0].header).celestial
        with fits.open(self.download_dir / "invvar" / f"{g}_invvar.fits") as hdul:
            if len(hdul) != 7:
                raise ValueError("Diagnostic requires a single native brick with seven HDUs")
            invvars = []
            invwcs = []
            for band, i in zip(bands, (2, 4, 6)):
                h = hdul[i]
                label = str(h.header.get("BAND", h.header.get("FILTER", band))).strip().lower()
                if label != band or h.data is None or h.data.ndim != 2:
                    raise ValueError("Unexpected inverse-variance band or layout")
                invvars.append(h.data.copy())
                invwcs.append(WCS(h.header).celestial)
        sigmas = []
        for inv in invvars:
            sigma = np.full(inv.shape, np.nan)
            good = np.isfinite(inv) & (inv > 0)
            sigma[good] = 1 / np.sqrt(inv[good])
            sigmas.append(sigma)
        psfs, _ = self._read_psfs(self.download_dir / "psf" / f"{g}_psf.fits")
        fig = plt.figure(figsize=(12, 16))
        try:
            for row, (label, arrays) in enumerate(zip(
                    ("image", "invvar", "sigma", "PSF"),
                    (images, invvars, sigmas, psfs))):
                for col, data in enumerate(arrays):
                    projection = wcs if row == 0 else invwcs[col] if row < 3 else None
                    ax = fig.add_subplot(4, 3, row * 3 + col + 1, projection=projection)
                    finite = data[np.isfinite(data)]
                    limits = np.percentile(finite, (21, 97) if row == 0 else (1, 99)) if finite.size else (0, 1)
                    ax.imshow(data, origin="lower", cmap="viridis", vmin=limits[0], vmax=limits[1])
                    ax.set_title(f"{bands[col]} {label}")
                    if row == 3:
                        ax.axis("off")
                    else:
                        ax.coords.grid(color="white", ls="dotted", alpha=0.5)
                        ax.set_xlabel("RA")
                        ax.set_ylabel("Dec")
            fig.tight_layout()
            out = self.plot_dir / f"{g}.png"
            fig.savefig(out, dpi=200)
        finally:
            plt.close(fig)
        print(f"Saved plot: {out}")
        return out

    # ---------- PSF SPLIT ----------
    def split_psf(self, target):
        g = self.clean_name(target)
        bands = ["g","r","z"]

        psf_file = self.download_dir / "psf" / f"{g}_psf.fits"

        psfs, headers = self._read_psfs(psf_file)
        for band, data, header in zip(bands, psfs, headers):
            header["BAND"] = band
            fits.writeto(self.download_dir / "psf" / f"{g}_PSF_{band}.fits",
                         data, header=header, overwrite=True)

        print(f"✔ Split PSF for {g}")

    @staticmethod
    def _read_psfs(path):
        # Retain the uploaded split_psf convention: one 2-D HDU per g/r/z band.
        # Reject ambiguous cubes rather than guessing the band axis.
        with fits.open(path) as hdul:
            planes = [h for h in hdul if h.data is not None and h.data.ndim == 2]
            if len(planes) != 3:
                raise ValueError("Expected three 2-D PSF HDUs in g/r/z order; inspect the FITS headers to identify other layouts")
            return [h.data.copy() for h in planes], [h.header.copy() for h in planes]

    @staticmethod
    def _image_header(header, band):
        result = WCS(header).celestial.to_header()
        for key in ('BUNIT', 'TELESCOP', 'INSTRUME', 'SURVEY', 'MAGZERO'):
            if key in header:
                result[key] = header[key]
        result['BAND'] = band
        result.add_history('Split from original DECaLS product without resampling.')
        return result

    def save_split_data(self, target, include_invvar=True):
        """Split science independently; inverse variance requires verified matching grids."""
        if include_invvar:
            self._check_image_invvar(target)
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
            if img is None or img.ndim != 3 or img.shape[0] != 3:
                raise ValueError("Expected image cube (3, ny, nx) in g/r/z order")
    
            for i, band in enumerate(bands):
                data = img[i,:,:]
    
                h = self._image_header(header, band)
    
                out = out_img / f"{g}_{band}.fits"
                fits.writeto(out, data, header=h, overwrite=True)
    
        if not include_invvar:
            print(f"Saved split science images for {g}; inverse variance not processed. Existing inverse-variance outputs, if any, are unchanged and not validated by this run.")
            return

        # ---------- INVVAR ----------
        with fits.open(invvar_file) as hdul:
            invs = [hdul[2].data, hdul[4].data, hdul[6].data]
    
            for i, band in enumerate(bands):
                data = invs[i]
    
                h = self._image_header(hdul[(2, 4, 6)[i]].header, band)
    
                out = out_inv / f"{g}_{band}_invvar.fits"
                fits.writeto(out, data, header=h, overwrite=True)
    
        print(f"✔ Saved split IMAGE + INVVAR for {g}")

    

    # ---------- RUN ----------
    def run(self, targets, size=1200,
            skip_download=True,
            force_download=False,
            do_plot=True,
            do_psf_split=True,
            do_save_split=True, include_invvar=False, ra=None, dec=None):

        if isinstance(targets, str):
            targets = [targets]

        targets = list(targets)
        if ra is not None or dec is not None:
            if len(targets) != 1 or ra is None or dec is None:
                raise ValueError('Explicit ra and dec require exactly one target and both coordinates.')
            self._degree_coord(ra, dec)

        for t in targets:
            print(f"\n=== {t} ===")

            # DOWNLOAD CONTROL (fixed)
            if force_download:
                print(f"🔁 Force downloading {t}")
                coord = self.download(t, size=size, ra=ra, dec=dec)

            elif skip_download and self.already_downloaded(t):
                coord = self.download(t, size=size, overwrite=False, ra=ra, dec=dec)

            else:
                coord = self.download(t, size=size, overwrite=not skip_download, ra=ra, dec=dec)

            if do_save_split:
                self.save_split_data(t, include_invvar=False)
            if do_plot:
                self.review(t, coord.ra.deg, coord.dec.deg)
                try:
                    self.plot(t)
                except (ValueError, OSError) as exc:
                    warnings.warn(f"{t}: full diagnostic unavailable ({exc}); g/r/z preview saved.", RuntimeWarning)
                    try:
                        self.plot_native_invvar(t)
                    except (ValueError, OSError) as native_exc:
                        warnings.warn(f"{t}: native diagnostic unavailable ({native_exc}).", RuntimeWarning)
            if do_psf_split:
                self.split_psf(t)
            if do_save_split and include_invvar:
                self.save_split_data(t, include_invvar=True)

