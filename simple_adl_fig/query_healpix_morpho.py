"""HEALPix version of query_morpho: DELVE-MC DR3 data for a set of nside=32
pixels, with morphocloud star/galaxy separation and NO calibration offsets.

Author: Mahdieh Navabi

query(profile, pixels, rmax, stars=True, galaxies=False)
    pixels : nside=32 HEALPix pixel numbers (e.g. the search pixel plus its
             neighbours, from hp.get_all_neighbours).

How it works
------------
1. Each nside=32 pixel is expanded to the nside=256 pixels that cover it
   (hp.query_polygon on its boundary, inclusive=True).
2. One TAP query on delve_dr3.mc_object selects
       ring256 IN (...)  AND rmag < rmax  AND brickuniq
3. Objects are then kept only if hp.ang2pix(32, ra, dec) is one of the
   requested pixels, so the result is exact pixel membership (the ring256
   set can slightly over-cover the pixel edges).
4. The same morphocloud classifier as query_morpho (morpho_separation.classify)
   splits stars and galaxies.

No calibration is applied: magnitudes are only dereddened (E(B-V)).
"""
import time
from functools import lru_cache

import healpy as hp
import numpy as np
import pyvo

from simple_adl_fig.morpho_separation import (
    TAP_URL, MAG_SENT, R_g, R_r, R_i, R_z,
    classify,
)

NSIDE_32 = 32
NSIDE_256 = 256


def nside32toring256(pixels):
    """All nside=256 (RING) pixels covering the given nside=32 pixels.

    Parameters
    ----------
    pixels : int or iterable of int, nside=32 pixel numbers

    Returns
    -------
    numpy int64 array of unique nside=256 pixel numbers
    """
    out = []
    for pixel_32 in np.atleast_1d(np.asarray(pixels, dtype=np.int64)):
        # Polygon from the pixel corners, as in the original nside32toring256
        vertices = np.transpose(hp.boundaries(NSIDE_32, int(pixel_32), nest=False))
        out.append(hp.query_polygon(NSIDE_256, vertices, inclusive=True, nest=False))
        # Disc around the centre covers the curved edges the 4-corner polygon misses
        vec = hp.pix2vec(NSIDE_32, int(pixel_32), nest=False)
        radius = hp.max_pixrad(NSIDE_32)
        out.append(hp.query_disc(NSIDE_256, vec, radius, inclusive=True, nest=False))
    return np.unique(np.concatenate(out)).astype(np.int64)


def _fetch_ring256(ring256_list, rmax):
    """TAP query of delve_dr3.mc_object for the given nside=256 pixels."""
    in_list = ",".join(str(int(p)) for p in ring256_list)
    adql = f"""
    SELECT objid, brick, brickuniq, ring256,
           ra, dec,
           gmag, gerr, grms, ndetg,
           rmag, rerr, rrms, ndetr,
           imag, ierr, irms, ndeti,
           zmag, zerr, zrms, ndetz,
           ebv, ndet, nalfdetiter,
           chi, sharp, prob, mag_auto, fwhm, ellipticity, asemi, bsemi
    FROM delve_dr3.mc_object
    WHERE ring256 IN ({in_list})
      AND rmag < {rmax}
    """
    svc = pyvo.dal.TAPService(TAP_URL)
    for attempt in range(1, 4):
        try:
            df = svc.run_async(adql, maxrec=10_000_000).to_table().to_pandas()
            break
        except Exception as exc:
            if attempt == 3:
                raise
            print(f"TAP query failed ({type(exc).__name__}), retrying ({attempt}/3)...")
            time.sleep(10 * attempt)

    # Same dereddening convention as morpho_separation.fetch (no calibration)
    df["ebv"] = df["ebv"].astype(float)
    for lo, R in (("g", R_g), ("r", R_r), ("i", R_i), ("z", R_z)):
        mag = df[f"{lo}mag"].astype(float)
        df[f"{lo}mag0"] = (mag - R * df["ebv"]).where(mag < MAG_SENT)
    return df


def _select_pixels(dl, pixels):
    """Keep rows whose position falls in one of the nside=32 pixels."""
    pix = hp.ang2pix(NSIDE_32, dl["ra"].to_numpy(dtype=float),
                     dl["dec"].to_numpy(dtype=float), lonlat=True)
    keep = np.isin(pix, np.asarray(pixels, dtype=np.int64))
    return dl[keep]


@lru_cache(maxsize=8)
def _labelled(pixels, rmax):
    """Fetch + pixel-select + classify once per pixel set; stars and galaxies share it."""
    ring256_list = nside32toring256(pixels)
    dl = _fetch_ring256(ring256_list, rmax)
    dl = dl[dl["brickuniq"].map({"t": True, "f": False, True: True, False: False}).astype(bool)]
    dl = _select_pixels(dl, pixels).reset_index(drop=True)
    if len(dl) == 0:
        return None
    return classify(dl)


def query(profile, pixels, rmax, stars=True, galaxies=False):
    """Return morphocloud-separated DELVE-MC objects inside the given HEALPix pixels.

    `profile` is accepted for compatibility with query_dl / query_morpho and is
    not used (the public Data Lab TAP service needs no profile).
    """
    key = tuple(int(p) for p in np.atleast_1d(pixels))
    labelled = _labelled(key, float(rmax))
    if labelled is None:
        empty = np.array([])
        return (empty, empty) if (stars and galaxies) else empty

    star_tab = labelled[labelled["STAR_LEAK1"]].reset_index(drop=True).to_records(index=False)
    gal_tab = labelled[~labelled["STAR_LEAK1"]].reset_index(drop=True).to_records(index=False)

    if stars and galaxies:
        return star_tab, gal_tab
    if stars:
        return star_tab
    if galaxies:
        return gal_tab
    return None
