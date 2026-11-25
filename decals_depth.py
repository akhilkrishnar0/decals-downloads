"""DR9 PSF-depth maps sampled onto the science image WCS."""
from pathlib import Path
import csv
import time
import numpy as np
import requests
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt

BASE = "https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr9"
BANDS = ("g", "r", "z")


def download(url, path, tries=3):
    path = Path(path)
    if path.is_file() and path.stat().st_size:
        try:
            with fits.open(path, memmap=False) as h:
                h.verify("exception")
                if any(x.data is not None for x in h):
                    return path
        except (OSError, ValueError):
            pass
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name("part-" + path.name)
    for attempt in range(tries):
        try:
            with requests.get(url, stream=True, timeout=(20, 180)) as response:
                response.raise_for_status()
                with tmp.open("wb") as f:
                    for block in response.iter_content(chunk_size=1024*1024):
                        if block:
                            f.write(block)
            with fits.open(tmp, memmap=False) as h:
                h.verify("exception")
                if not any(x.data is not None for x in h):
                    raise ValueError("Empty FITS")
            tmp.replace(path)
            return path
        except (requests.RequestException, OSError, ValueError) as exc:
            tmp.unlink(missing_ok=True)
            if isinstance(exc, requests.HTTPError) and exc.response is not None and exc.response.status_code == 404:
                raise
            if attempt == tries - 1:
                raise
            time.sleep(2*(attempt+1))


def get_bricks(cache):
    path = download(f"{BASE}/survey-bricks.fits.gz", cache/"survey-bricks.fits.gz")
    with fits.open(path, memmap=False) as h:
        table = h[1].data
        return [(str(n).strip() if not isinstance(n, bytes) else n.decode().strip(),
                 float(a), float(b), float(c), float(d))
                for n, a, b, c, d in zip(table["BRICKNAME"], table["RA1"], table["RA2"],
                                        table["DEC1"], table["DEC2"])]


def assign_bricks(wcs, shape, bricks):
    """Assign each science pixel centre to its unique survey brick bounds."""
    ny, nx = shape
    # Filter by full pixel-centre Dec/RA bounds; RA is unwrapped around CRVAL.
    ref = float(wcs.wcs.crval[0]) % 360
    def unwrap(ra):
        return (ra-ref+180) % 360-180
    # Compute world coordinates in strips to avoid a full sky-coordinate mesh.
    bounds = [np.inf, -np.inf, np.inf, -np.inf]
    for y0 in range(0, ny, 128):
        yy, xx = np.mgrid[y0:min(y0+128, ny), :nx]
        ra, dec = wcs.pixel_to_world_values(xx, yy)
        a = unwrap(ra)
        good = np.isfinite(a) & np.isfinite(dec)
        if good.any():
            bounds = [min(bounds[0], a[good].min()), max(bounds[1], a[good].max()),
                      min(bounds[2], dec[good].min()), max(bounds[3], dec[good].max())]
    candidates = []
    for name, a, b, c, d in bricks:
        if d < bounds[2] or c > bounds[3]:
            continue
        lo = float(unwrap(a))
        width = (b-a) % 360
        if width == 0 and b != a:
            width = 360
        if any(lo+s <= bounds[1] and lo+width+s >= bounds[0] for s in (-360, 0, 360)):
            candidates.append((name, a, b, c, d))
    ids = np.full(shape, -1, dtype=np.int32)
    for y0 in range(0, ny, 128):
        yy, xx = np.mgrid[y0:min(y0+128, ny), :nx]
        ra, dec = wcs.pixel_to_world_values(xx, yy)
        block = ids[y0:y0+len(yy)]
        for i, (_, a, b, c, d) in enumerate(candidates):
            width = (b-a) % 360
            if width == 0 and b != a:
                width = 360
            mask = ((ra-a)%360 < width) & (dec >= c) & (dec < d)
            block[mask & (block == -1)] = i
    return ids, candidates


def sample_depth(path, ra, dec):
    with fits.open(path, memmap=False) as h:
        ext = next(x for x in h if x.data is not None and x.data.ndim == 2)
        w = WCS(ext.header).celestial
        x, y = w.world_to_pixel_values(ra, dec)
        data = ext.data
        good = (np.isfinite(x) & np.isfinite(y) & (x >= -.5) & (y >= -.5)
                & (x < data.shape[1]-.5) & (y < data.shape[0]-.5))
        result = np.full(ra.shape, np.nan, dtype=np.float32)
        xi = np.floor(x[good]+.5).astype(int)
        yi = np.floor(y[good]+.5).astype(int)
        values = np.asarray(data[yi, xi], dtype=np.float32)
        values[~np.isfinite(values) | (values <= 0)] = np.nan
        result[good] = values
        return result


def make_depth_maps(science_file, cache, output, plot_file, region="south"):
    """Download all intersecting DR9 depth bricks for the selected survey region.

    Nearest-neighbour sampling at each science pixel centre; no summing of
    overlapping bricks. Missing or nonpositive depth is NaN.
    """
    if region not in ("south", "north"):
        raise ValueError("region must be 'south' or 'north'")
    cache, output, plot_file = Path(cache), Path(output), Path(plot_file)
    output.mkdir(parents=True, exist_ok=True)
    plot_file.parent.mkdir(parents=True, exist_ok=True)
    with fits.open(science_file) as h:
        cube = h[0].data
        if cube is None or cube.ndim != 3 or cube.shape[0] != 3:
            raise ValueError("Expected g/r/z science cube (3, ny, nx)")
        images = cube.copy()
        wcs = WCS(h[0].header).celestial
    shape = images.shape[1:]
    ids, bricks = assign_bricks(wcs, shape, get_bricks(cache))
    used = np.unique(ids[ids >= 0])
    if not len(used):
        raise ValueError("No survey bricks intersect the science grid")
    print("DR9 bricks:", ", ".join(bricks[i][0] for i in used))
    coverage = []
    results = {}
    fig = plt.figure(figsize=(12, 12))
    try:
        for col, band in enumerate(BANDS):
            depth = np.full(shape, np.nan, dtype=np.float32)
            for i in used:
                brick = bricks[i][0]
                name = f"legacysurvey-{brick}-depth-{band}.fits.fz"
                url = f"{BASE}/{region}/coadd/{brick[:3]}/{brick}/{name}"
                yy, xx = np.nonzero(ids == i)
                try:
                    path = download(url, cache/region/brick/name)
                except requests.HTTPError as exc:
                    if exc.response is None or exc.response.status_code != 404:
                        raise
                    print(f"No {region} depth: {brick} {band}")
                    coverage.append((brick, band, len(xx), 0, url))
                    continue
                ra, dec = wcs.pixel_to_world_values(xx, yy)
                values = sample_depth(path, ra, dec)
                depth[yy, xx] = values
                coverage.append((brick, band, len(xx), int(np.isfinite(values).sum()), url))
            sigma = np.full(shape, np.nan, dtype=np.float32)
            good = np.isfinite(depth) & (depth > 0)
            sigma[good] = 1 / np.sqrt(depth[good])
            for name, array, unit in (("depth", depth, "1/nanomaggy2"),
                                      ("sigma_psf", sigma, "nanomaggy")):
                header = wcs.to_header()
                header["BAND"] = band
                header["BUNIT"] = unit
                header["REGION"] = region
                header.add_history("DR9 PSF depth, nearest-neighbour sampled at science pixel centres.")
                header.add_history("Not per-pixel inverse variance or aperture-sum uncertainty.")
                path = output / f"{band}_{name}.fits"
                fits.writeto(path, array, header, overwrite=True)
                results[f"{band}_{name}"] = path
            print(f"{band}: {good.mean():.1%} positive depth coverage")
            for row, (name, array) in enumerate((("image", images[col]),
                                                 ("PSF depth", depth),
                                                 ("PSF flux sigma", sigma))):
                ax = fig.add_subplot(3, 3, row*3+col+1, projection=wcs)
                vals = array[np.isfinite(array)]
                lo, hi = np.percentile(vals, (21, 97) if row == 0 else (1, 99)) if vals.size else (0, 1)
                if hi <= lo:
                    hi = lo + max(abs(lo)*1e-6, 1e-12)
                ax.imshow(array, origin="lower", cmap="viridis", vmin=lo, vmax=hi)
                ax.set_title(f"{band} {name}")
                ax.coords.grid(color="white", ls="dotted", alpha=.5)
                ax.set_xlabel("RA")
                ax.set_ylabel("Dec")
        fig.tight_layout()
        fig.savefig(plot_file, dpi=180)
    finally:
        plt.close(fig)
    with (output/"coverage.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(("brick", "band", "assigned_pixels", "valid_pixels", "source_url"))
        writer.writerows(coverage)
    results["plot"] = plot_file
    print(f"Saved depth diagnostic: {plot_file}")
    return results
